# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build spatial distribution of industries from Hotmaps database.

Inputs
-------

- ``resources/regions_onshore_elec_s{simpl}_{clusters}.geojson``
- ``resources/pop_layout_elec_s{simpl}_{clusters}.csv``

Outputs
-------

- ``resources/industrial_distribution_key_elec_s{simpl}_{clusters}.csv``

Description
-------

This rule uses the `Hotmaps database <https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database>`. After removing entries without valid locations, it assigns each industrial site to a bus region based on its location.
Then, it calculates the nodal distribution key for each sector based on the emissions of the industrial sites in each region. This leads to a distribution key of 1 if there is only one bus per country and <1 if there are multiple buses per country. The sum over buses of one country is 1.

The following subcategories of industry are considered:
- Iron and steel
- Cement
- Refineries
- Paper and printing
- Chemical industry
- Glass
- Non-ferrous metals
- Non-metallic mineral products
- Other non-classified
Furthermore, the population distribution is added
- Population
"""

import logging
import uuid
from itertools import product

import country_converter as coco
import geopandas as gpd
import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()


def locate_missing_industrial_sites(df):
    """
    Locate industrial sites without valid locations based on city and
    countries.

    Should only be used if the model's spatial resolution is coarser
    than individual cities.
    """
    try:
        from geopy.extra.rate_limiter import RateLimiter
        from geopy.geocoders import Nominatim
    except ImportError:
        raise ModuleNotFoundError(
            "Optional dependency 'geopy' not found."
            "Install via 'conda install -c conda-forge geopy'"
            "or set 'industry: hotmaps_locate_missing: false'."
        )

    locator = Nominatim(user_agent=str(uuid.uuid4()))
    geocode = RateLimiter(locator.geocode, min_delay_seconds=2)

    def locate_missing(s):
        if pd.isna(s.City) or s.City == "CONFIDENTIAL":
            return None

        loc = geocode([s.City, s.Country], geometry="wkt")
        if loc is not None:
            logger.debug(f"Found:\t{loc}\nFor:\t{s['City']}, {s['Country']}\n")
            return f"POINT({loc.longitude} {loc.latitude})"
        else:
            return None

    missing = df.index[df.geom.isna()]
    df.loc[missing, "coordinates"] = df.loc[missing].apply(locate_missing, axis=1)

    # report stats
    num_still_missing = df.coordinates.isna().sum()
    num_found = len(missing) - num_still_missing
    share_missing = len(missing) / len(df) * 100
    share_still_missing = num_still_missing / len(df) * 100
    logger.warning(
        f"Found {num_found} missing locations. \nShare of missing locations reduced from {share_missing:.2f}% to {share_still_missing:.2f}%."
    )

    return df


def prepare_hotmaps_database(regions):
    """
    Load hotmaps database of industrial sites and map onto bus regions.
    """
    df = pd.read_csv(snakemake.input.hotmaps_industrial_database, sep=";", index_col=0)

    df[["srid", "coordinates"]] = df.geom.str.split(";", expand=True)

    if snakemake.params.hotmaps_locate_missing:
        df = locate_missing_industrial_sites(df)

    # remove those sites without valid locations
    df.drop(df.index[df.coordinates.isna()], inplace=True)

    df["coordinates"] = gpd.GeoSeries.from_wkt(df["coordinates"])

    gdf = gpd.GeoDataFrame(df, geometry="coordinates", crs="EPSG:4326")

    gdf = gpd.sjoin(gdf, regions, how="inner", predicate="within")

    gdf.rename(columns={"index_right": "bus"}, inplace=True)
    gdf["country"] = gdf.bus.str[:2]

    # the .sjoin can lead to duplicates if a geom is in two overlapping regions
    if gdf.index.duplicated().any():
        # get all duplicated entries
        duplicated_i = gdf.index[gdf.index.duplicated()]
        # convert from raw data country name to iso-2-code
        code = cc.convert(gdf.loc[duplicated_i, "Country"], to="iso2")  # noqa: F841
        # screen out malformed country allocation
        gdf_filtered = gdf.loc[duplicated_i].query("country == @code")
        # concat not duplicated and filtered gdf
        gdf = pd.concat([gdf.drop(duplicated_i), gdf_filtered])

    return gdf


def build_nodal_distribution_key(hotmaps, regions, countries):
    """
    Build nodal distribution keys for each sector.
    """
    sectors = hotmaps.Subsector.unique()

    keys = pd.DataFrame(index=regions.index, columns=sectors, dtype=float)

    pop = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    pop["country"] = pop.index.str[:2]
    ct_total = pop.total.groupby(pop["country"]).sum()
    keys["population"] = pop.total / pop.country.map(ct_total)

    for sector, country in product(sectors, countries):
        regions_ct = regions.index[regions.index.str.contains(country)]

        facilities = hotmaps.query("country == @country and Subsector == @sector")

        if not facilities.empty:
            emissions = facilities["Emissions_ETS_2014"].fillna(
                hotmaps["Emissions_EPRTR_2014"].dropna()
            )
            if emissions.sum() == 0:
                key = pd.Series(1 / len(facilities), facilities.index)
            else:
                # BEWARE: this is a strong assumption
                emissions = emissions.fillna(emissions.mean())
                key = emissions / emissions.sum()
            key = key.groupby(facilities.bus).sum().reindex(regions_ct, fill_value=0.0)
        else:
            key = keys.loc[regions_ct, "population"]

        keys.loc[regions_ct, sector] = key

    return keys


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_distribution_key",
            simpl="",
            clusters=128,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    countries = snakemake.params.countries

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    hotmaps = prepare_hotmaps_database(regions)

    keys = build_nodal_distribution_key(hotmaps, regions, countries)

    keys.to_csv(snakemake.output.industrial_distribution_key)
