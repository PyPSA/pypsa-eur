# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2021-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Compute biogas and solid biomass potentials for each clustered model region
using data from JRC ENSPRESO.
"""

import logging

logger = logging.getLogger(__name__)
import geopandas as gpd
import numpy as np
import pandas as pd

AVAILABLE_BIOMASS_YEARS = [2010, 2020, 2030, 2040, 2050]


def build_nuts_population_data(year=2013):
    pop = pd.read_csv(
        snakemake.input.nuts3_population,
        sep=r"\,| \t|\t",
        engine="python",
        na_values=[":"],
        index_col=1,
    )[str(year)]

    # only countries
    pop.drop("EU28", inplace=True)

    # mapping from Cantons to NUTS3
    cantons = pd.read_csv(snakemake.input.swiss_cantons)
    cantons = cantons.set_index(cantons.HASC.str[3:]).NUTS
    cantons = cantons.str.pad(5, side="right", fillchar="0")

    # get population by NUTS3
    swiss = pd.read_excel(
        snakemake.input.swiss_population, skiprows=3, index_col=0
    ).loc["Residents in 1000"]
    swiss = swiss.rename(cantons).filter(like="CH")

    # aggregate also to higher order NUTS levels
    swiss = [swiss.groupby(swiss.index.str[:i]).sum() for i in range(2, 6)]

    # merge Europe + Switzerland
    pop = pd.concat([pop, pd.concat(swiss)]).to_frame("total")

    # add missing manually
    pop["AL"] = 2893
    pop["BA"] = 3871
    pop["RS"] = 7210

    pop["ct"] = pop.index.str[:2]

    return pop


def enspreso_biomass_potentials(year=2020, scenario="ENS_Low"):
    """
    Loads the JRC ENSPRESO biomass potentials.

    Parameters
    ----------
    year : int
        The year for which potentials are to be taken.
        Can be {2010, 2020, 2030, 2040, 2050}.
    scenario : str
        The scenario. Can be {"ENS_Low", "ENS_Med", "ENS_High"}.

    Returns
    -------
    pd.DataFrame
        Biomass potentials for given year and scenario
        in TWh/a by commodity and NUTS2 region.
    """
    glossary = pd.read_excel(
        str(snakemake.input.enspreso_biomass),
        sheet_name="Glossary",
        usecols="B:D",
        skiprows=1,
        index_col=0,
    )

    df = pd.read_excel(
        str(snakemake.input.enspreso_biomass),
        sheet_name="ENER - NUTS2 BioCom E",
        usecols="A:H",
    )

    df["group"] = df["E-Comm"].map(glossary.group)
    df["commodity"] = df["E-Comm"].map(glossary.description)

    to_rename = {
        "NUTS2 Potential available by Bio Commodity": "potential",
        "NUST2": "NUTS2",
    }
    df.rename(columns=to_rename, inplace=True)

    # fill up with NUTS0 if NUTS2 is not given
    df.NUTS2 = df.apply(lambda x: x.NUTS0 if x.NUTS2 == "-" else x.NUTS2, axis=1)

    # convert PJ to TWh
    df.potential /= 3.6
    df.Unit = "TWh/a"

    dff = df.query("Year == @year and Scenario == @scenario")

    bio = dff.groupby(["NUTS2", "commodity"]).potential.sum().unstack()

    # currently Serbia and Kosovo not split, so aggregate
    bio.loc["RS"] += bio.loc["XK"]
    bio.drop("XK", inplace=True)

    return bio


def disaggregate_nuts0(bio):
    """
    Some commodities are only given on NUTS0 level. These are disaggregated
    here using the NUTS2 population as distribution key.

    Parameters
    ----------
    bio : pd.DataFrame
        from enspreso_biomass_potentials()

    Returns
    -------
    pd.DataFrame
    """
    pop = build_nuts_population_data()

    # get population in nuts2
    pop_nuts2 = pop.loc[pop.index.str.len() == 4]
    by_country = pop_nuts2.total.groupby(pop_nuts2.ct).sum()
    pop_nuts2["fraction"] = pop_nuts2.total / pop_nuts2.ct.map(by_country)

    # distribute nuts0 data to nuts2 by population
    bio_nodal = bio.loc[pop_nuts2.ct]
    bio_nodal.index = pop_nuts2.index
    bio_nodal = bio_nodal.mul(pop_nuts2.fraction, axis=0)

    # update inplace
    bio.update(bio_nodal)

    return bio


def build_nuts2_shapes():
    """
    - load NUTS2 geometries
    - add RS, AL, BA country shapes (not covered in NUTS 2013)
    - consistently name ME, MK
    """
    nuts2 = gpd.GeoDataFrame(
        gpd.read_file(snakemake.input.nuts2).set_index("id").geometry
    )

    countries = gpd.read_file(snakemake.input.country_shapes).set_index("name")
    missing_iso2 = countries.index.intersection(["AL", "RS", "BA"])
    missing = countries.loc[missing_iso2]

    nuts2.rename(index={"ME00": "ME", "MK00": "MK"}, inplace=True)

    return pd.concat([nuts2, missing])


def area(gdf):
    """
    Returns area of GeoDataFrame geometries in square kilometers.
    """
    return gdf.to_crs(epsg=3035).area.div(1e6)


def convert_nuts2_to_regions(bio_nuts2, regions):
    """
    Converts biomass potentials given in NUTS2 to PyPSA-Eur regions based on
    the overlay of both GeoDataFrames in proportion to the area.

    Parameters
    ----------
    bio_nuts2 : gpd.GeoDataFrame
        JRC ENSPRESO biomass potentials indexed by NUTS2 shapes.
    regions : gpd.GeoDataFrame
        PyPSA-Eur clustered onshore regions

    Returns
    -------
    gpd.GeoDataFrame
    """
    # calculate area of nuts2 regions
    bio_nuts2["area_nuts2"] = area(bio_nuts2)

    overlay = gpd.overlay(regions, bio_nuts2, keep_geom_type=True)

    # calculate share of nuts2 area inside region
    overlay["share"] = area(overlay) / overlay["area_nuts2"]

    # multiply all nuts2-level values with share of nuts2 inside region
    adjust_cols = overlay.columns.difference(
        {"name", "area_nuts2", "geometry", "share"}
    )
    overlay[adjust_cols] = overlay[adjust_cols].multiply(overlay["share"], axis=0)

    bio_regions = overlay.dissolve("name", aggfunc="sum")

    bio_regions.drop(["area_nuts2", "share"], axis=1, inplace=True)

    return bio_regions


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_biomass_potentials",
            simpl="",
            clusters="5",
            planning_horizons=2050,
        )

    overnight = snakemake.config["foresight"] == "overnight"
    params = snakemake.params.biomass
    investment_year = int(snakemake.wildcards.planning_horizons)
    year = params["year"] if overnight else investment_year
    scenario = params["scenario"]

    if year > 2050:
        logger.info("No biomass potentials for years after 2050, using 2050.")
        max_year = max(AVAILABLE_BIOMASS_YEARS)
        enspreso = enspreso_biomass_potentials(max_year, scenario)

    elif year not in AVAILABLE_BIOMASS_YEARS:
        before = int(np.floor(year / 10) * 10)
        after = int(np.ceil(year / 10) * 10)
        logger.info(
            f"No biomass potentials for {year}, interpolating linearly between {before} and {after}."
        )

        enspreso_before = enspreso_biomass_potentials(before, scenario)
        enspreso_after = enspreso_biomass_potentials(after, scenario)

        fraction = (year - before) / (after - before)

        enspreso = enspreso_before + fraction * (enspreso_after - enspreso_before)

    else:
        logger.info(f"Using biomass potentials for {year}.")
        enspreso = enspreso_biomass_potentials(year, scenario)

    enspreso = disaggregate_nuts0(enspreso)

    nuts2 = build_nuts2_shapes()

    df_nuts2 = gpd.GeoDataFrame(nuts2.geometry).join(enspreso)

    regions = gpd.read_file(snakemake.input.regions_onshore)

    df = convert_nuts2_to_regions(df_nuts2, regions)

    df.to_csv(snakemake.output.biomass_potentials_all)

    grouper = {v: k for k, vv in params["classes"].items() for v in vv}
    df = df.groupby(grouper, axis=1).sum()

    df *= 1e6  # TWh/a to MWh/a
    df.index.name = "MWh/a"

    df.to_csv(snakemake.output.biomass_potentials)
