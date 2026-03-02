# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Builds the electricity demand for base regions based on population and GDP.
"""

import logging

import country_converter as coco
import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import rasterio as rio
import xarray as xr
from rasterstats import zonal_stats

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()


def normed(s: pd.Series) -> pd.Series:
    return s / s.sum()


def redistribute_attribute(
    orig: gpd.GeoDataFrame, dest: gpd.GeoDataFrame, attr: str
) -> pd.Series:
    """
    Redistributes an attribute from origin GeoDataFrame to destination GeoDataFrame based on overlapping areas.

    Computes the intersection between geometries and proportionally assigns the attribute values
    according to the share of overlapping area.

    Parameters
    ----------
    orig : gpd.GeoDataFrame
        Source GeoDataFrame containing the attribute to redistribute
    dest : gpd.GeoDataFrame
        Target GeoDataFrame to receive redistributed values
    attr : str
        Name of the attribute column to redistribute

    Returns
    -------
    pd.Series
        Redistributed attribute values indexed by destination's 'name' column
    """
    if orig.crs != 3035:
        orig = orig.to_crs(epsg=3035)
    if dest.crs != 3035:
        dest = dest.to_crs(epsg=3035)
    orig["area_orig"] = orig.area
    overlay = gpd.overlay(dest, orig, keep_geom_type=False)
    overlay["share"] = overlay.area / overlay.area_orig
    overlay["attr"] = overlay[attr].mul(overlay.share)
    return overlay.dissolve("name", aggfunc="sum")["attr"]


def energy_atlas_distribution_keys(
    raster_fn: str, regions: gpd.GeoDataFrame
) -> pd.Series:
    """
    Calculate distribution keys for regions based on energy atlas raster data.

    Parameters
    ----------
    fn : str
        File path to the raster data (GeoTIFF format).
    regions : gpd.GeoDataFrame
        GeoDataFrame containing the regions with a 'country' column.

    Returns
    -------
    pd.Series
        Series of distribution keys indexed by region names.
        Sum of keys per country equals 1 (unless omitted islands in regions).
    """

    raster = rio.open(raster_fn)
    band = raster.read(1)

    distribution_keys = []
    for country, group in regions.groupby("country"):
        if country not in cc.EU27as("ISO2").ISO2.to_list():
            continue
        weights = zonal_stats(
            group, band, affine=raster.transform, nodata=-1, stats="sum"
        )
        group["weights"] = [w["sum"] for w in weights]
        distribution_keys.append(group["weights"] / group["weights"].sum())
    distribution_keys = pd.concat(distribution_keys).reindex(regions.index)
    return distribution_keys


def gb_distribution_keys(
    excel_fn: str, geojson_fn: str, regions: gpd.GeoDataFrame
) -> pd.Series:
    """
    Calculate distribution keys for Great Britain regions based on electricity consumption statistics.

    Parameters
    ----------
    excel_fn : str
        File path to the Excel file containing electricity consumption by local authority districts.
    geojson_fn : str
        File path to the GeoJSON file containing geometries of local authority districts.
    regions : gpd.GeoDataFrame
        GeoDataFrame containing the network regions.

    Returns
    -------
    pd.Series
        Series of distribution keys indexed by region names.
    """

    df = pd.read_excel(excel_fn, skiprows=4, sheet_name="2019")
    df = df.loc[~df["Local authority"].isin(["All local authorities", "Unallocated"])]
    gdf = gpd.read_file(geojson_fn).to_crs(epsg=3035)
    gdf = gdf.rename(columns={"LAD24CD": "Code"}).merge(df, on="Code")

    attr = "Total consumption\n(GWh):\nAll meters"
    redistributed = redistribute_attribute(gdf, regions.reset_index(drop=True), attr)
    distribution_keys = normed(redistributed)
    return distribution_keys


def nuts3_distribution_keys(
    nuts3_fn: str, distribution_key: dict[str, float], regions: gpd.GeoDataFrame
) -> pd.Series:
    """
    Calculate distribution keys for regions based on NUTS3 data (GDP and population).

    Parameters
    ----------
    nuts3_fn : str
        File path to the NUTS3 GeoJSON file containing geometries and attributes.
    distribution_key : dict[str, float]
        Weights for GDP and population in the distribution key calculation.
        Example: {"gdp": 0.6, "pop": 0.4}
    regions : gpd.GeoDataFrame
        GeoDataFrame containing the regions with a 'country' column.

    Returns
    -------
    pd.Series
        Series of distribution keys indexed by region names.
    """

    gdp_weight = distribution_key.get("gdp", 0.6)
    pop_weight = distribution_key.get("pop", 0.4)

    nuts3 = gpd.read_file(nuts3_fn).to_crs(epsg=3035)
    nuts3.rename(columns={"name": "nuts3_name"}, inplace=True)

    regions["pop"] = redistribute_attribute(
        nuts3, regions.reset_index(drop=True), "pop"
    )
    regions["gdp"] = redistribute_attribute(
        nuts3, regions.reset_index(drop=True), "gdp"
    )

    nuts3_keys = []
    for country, group in regions.groupby("country"):
        factors = normed(
            gdp_weight * normed(group["gdp"]) + pop_weight * normed(group["pop"])
        )
        nuts3_keys.append(factors)
    return pd.concat(nuts3_keys).reindex(regions.index)


def upsample_load(
    n: pypsa.Network,
    regions_fn: str,
    load_fn: str,
    raster_fn: str,
    gb_excel_fn: str,
    gb_geojson_fn: str,
    nuts3_fn: str,
    distribution_key: dict[str, float],
    substation_only: bool = False,
) -> xr.DataArray:
    regions = gpd.read_file(regions_fn).set_index("name", drop=False).to_crs(epsg=3035)

    if substation_only:
        substation_lv_i = n.buses.index[n.buses["substation_lv"]]
        regions = regions.reindex(substation_lv_i)
    load = pd.read_csv(load_fn, index_col=0, parse_dates=True)

    ea_keys = energy_atlas_distribution_keys(raster_fn, regions)
    gb_keys = gb_distribution_keys(gb_excel_fn, gb_geojson_fn, regions)
    nuts3_keys = nuts3_distribution_keys(nuts3_fn, distribution_key, regions)

    factors = ea_keys.combine_first(gb_keys).combine_first(nuts3_keys)

    # sanitize: need to renormalize since `gb_keys` only cover Great Britain
    # and Northern Ireland is taken from `nuts3_keys`
    if "GB" in regions.country:
        uk_regions_i = regions.query("country == 'GB'").index
        uk_weights = factors.loc[uk_regions_i].sum()
        factors.loc[uk_regions_i] /= uk_weights

    data_arrays = []

    for cntry, group in regions.geometry.groupby(regions.country):
        if cntry not in load.columns:
            logger.warning(f"Cannot upsample load for {cntry}: no load data defined")
            continue

        load_ct = load[cntry]
        factors_ct = factors.loc[group.index]

        data_arrays.append(
            xr.DataArray(
                factors_ct.values * load_ct.values[:, np.newaxis],
                dims=["time", "bus"],
                coords={"time": load_ct.index.values, "bus": factors_ct.index.values},
            )
        )

    return xr.concat(data_arrays, dim="bus")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_electricity_demand_base")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params

    n = pypsa.Network(snakemake.input.base_network)

    load = upsample_load(
        n,
        regions_fn=snakemake.input.regions,
        load_fn=snakemake.input.load,
        raster_fn=snakemake.input.raster,
        gb_excel_fn=snakemake.input.gb_excel,
        gb_geojson_fn=snakemake.input.gb_geojson,
        nuts3_fn=snakemake.input.nuts3,
        distribution_key=params.distribution_key,
        substation_only=params.substation_only,
    )

    load.name = "electricity demand (MW)"
    comp = dict(zlib=True, complevel=9, least_significant_digit=5)
    load.to_netcdf(snakemake.output[0], encoding={load.name: comp})
