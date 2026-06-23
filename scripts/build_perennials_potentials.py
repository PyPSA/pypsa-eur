# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Reproject NUTS2-level 1st-generation (1G) biofuel and perennial crop yields
(from ``build_perennials_crop_yields_nuts2.py``) onto the clustered network
regions, for use by ``add_perennials()`` in ``prepare_sector_network.py``.

NUTS2 yields (MWh/ha/y for 1G biofuels, t/ha/y for perennials) are mapped to
clustered regions via an area-weighted overlay (NUTS2 geometries intersected
with cluster region geometries, weighted by intersection area); NUTS2
regions not covered by Eurostat data (non-EU countries, small islands,
city-states) are first filled from the nearest valid NUTS2 centroid.
Resulting columns are then grouped into biomass classes
(``resolve_biomass_classes()``) to match the class structure used by
``build_biomass_potentials.py``.

``add_perennials()`` later combines this script's clustered-region 1G yields
with the ENSPRESO 1G biomass potential (MWh/y) to back out the land area
available for conversion to perennial grasses, and converts that area into a
CO2 sequestration potential via ``perennials.potential_co2`` (tCO2/ha/y).

Outputs a single CSV with one column per crop class (cereals, sugar beet,
rapeseed, perennials) indexed by clustered region name.
"""

import logging

import geopandas as gpd
import pandas as pd
from _helpers import configure_logging, resolve_biomass_classes, set_scenario_config

logger = logging.getLogger(__name__)


def build_nuts2_shapes():
    """
    - load NUTS2 geometries
    - add RS, AL, BA country shapes (not covered in NUTS 2013)
    - consistently name ME, MK
    """
    nuts2 = gpd.GeoDataFrame(
        gpd.read_file(snakemake.input.nuts2).set_index("NUTS_ID").geometry
    )

    countries = gpd.read_file(snakemake.input.country_shapes).set_index("name")
    missing_iso2 = countries.index.intersection(["AL", "RS", "XK", "BA"])
    missing = countries.loc[missing_iso2]

    nuts2.rename(index={"ME00": "ME", "MK00": "MK"}, inplace=True)

    return pd.concat([nuts2, missing])


def area(gdf):
    return gdf.to_crs(epsg=3035).area.div(1e6)


def impute_missing_values(df_nuts2, missing_shapes, yield_cols):

    # Keep only rows that have valid yields (drop NaN rows!)
    df_valid = df_nuts2.dropna(subset=yield_cols).copy()

    # Project to 3035
    valid_proj = df_valid.to_crs(3035)
    missing_proj = missing_shapes.to_crs(3035)

    valid_centroids = valid_proj.centroid
    missing_centroids = missing_proj.centroid

    imputed_rows = []

    for missing_id, c_geom in missing_centroids.items():

        # Distance to only VALID NUTS2 rows
        dists = valid_centroids.distance(c_geom)
        nearest_nuts2 = dists.idxmin()

        # Copy yields
        attrs = df_valid.loc[nearest_nuts2, yield_cols]

        # Create new row
        new_row = missing_shapes.copy().loc[[missing_id]]
        for col in yield_cols:
            new_row[col] = attrs[col]

        imputed_rows.append(new_row)

    return pd.concat(imputed_rows)


def convert_nuts2_to_regions_yields(df_nuts2, regions, yield_cols=None):
    """
    Convert NUTS2-level yields (intensive) to PyPSA regions using:

        y_n = Σ_i (y_i * A_i∩n) / Σ_i A_i∩n

    Only NUTS2 rows with non-NaN yields are used. Regions with no
    overlapping valid NUTS2 get NaN.
    """

    nuts = df_nuts2.copy()
    regs = regions.copy()

    # Ensure ID column for NUTS2
    nuts["NUTS_ID"] = nuts.index

    # Identify yield columns if not given
    if yield_cols is None:
        yield_cols = nuts.columns.difference(["geometry", "NUTS_ID"])

    # 1) Keep only NUTS2 rows that have at least one non-NaN yield
    nuts_valid = nuts.dropna(subset=yield_cols, how="all")

    # 2) Reproject to equal-area CRS for areas
    nuts_valid = nuts_valid.to_crs(3035)
    regs = regs.to_crs(3035)

    # 3) Overlay: intersection between regions and valid NUTS2
    overlay = gpd.overlay(regs, nuts_valid, keep_geom_type=True)

    # Return empty if nothing overlaps
    if overlay.empty:
        return pd.DataFrame(index=regs["name"].values, columns=yield_cols, dtype=float)

    # 4) Area of intersections (m² → km²)
    overlay["area_intersection"] = overlay.geometry.area

    # Optional: drop absurdly tiny slivers
    overlay = overlay[overlay["area_intersection"] > 0]

    # 5) Numerators: y_i * A_i∩n
    numerators = overlay[yield_cols].multiply(overlay["area_intersection"], axis=0)

    # 6) Sum per region
    numer_by_region = numerators.groupby(overlay["name"]).sum(min_count=1)
    denom_by_region = overlay.groupby("name")["area_intersection"].sum()

    # 7) Weighted average
    yields_regions = numer_by_region.div(denom_by_region, axis=0)

    # 8) Ensure one row per region (regions with no valid overlaps → NaN)
    yields_regions = yields_regions.reindex(regs["name"].values)

    return yields_regions


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_perennials_potentials",
            clusters="39",
            planning_horizons=2050,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    regions = gpd.read_file(snakemake.input.regions_onshore)
    nuts2 = build_nuts2_shapes()

    yields = pd.read_csv(snakemake.input.perennials_yields_1G_biofuels, index_col=0)

    df_nuts2 = gpd.GeoDataFrame(nuts2.geometry).join(yields)

    # Address the missign countries (df_nuts2 contains all NUTS2 + missing shapes (with NaNs))
    missing_countries = ["AL", "RS", "BA", "XK"]
    missing_shapes = df_nuts2.loc[missing_countries, ["geometry"]]

    # Impute yields for missing shapes basest of nearest valid Nuts2
    yield_cols = [
        "Bioethanol barley, wheat, grain maize, oats, other cereals and rye",
        "Sugar from sugar beet",
        "Rape seed",
        "perennials",
    ]
    imputed_missing = impute_missing_values(df_nuts2, missing_shapes, yield_cols)
    df_nuts2 = pd.concat([
        df_nuts2.drop(index=missing_countries),
        imputed_missing
    ])

    # convert nuts2 yields to regions
    df = convert_nuts2_to_regions_yields(df_nuts2, regions)

    params = snakemake.params.biomass

    classes = resolve_biomass_classes(
        params["classes"], snakemake.config["sector"].get("perennials", False)
    )
    grouper = {v: k for k, vv in classes.items() for v in vv}

    # keep the column 'perennials' which is otherwise dropped
    for col in df.columns:
        if col not in grouper:
            grouper[col] = col

    df = df.T.groupby(grouper).sum().T

    df.index.name = "name"
    df.to_csv(snakemake.output.csv_file)
