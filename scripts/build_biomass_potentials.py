# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2021-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Compute biogas and solid biomass potentials for each clustered model region
using data from JRC ENSPRESO.
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)
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
    pop_nuts2 = pop.loc[pop.index.str.len() == 4].copy()
    by_country = pop_nuts2.total.groupby(pop_nuts2.ct).sum()
    pop_nuts2["fraction"] = pop_nuts2.total / pop_nuts2.ct.map(by_country)

    # distribute nuts0 data to nuts2 by population
    bio_nodal = bio.loc[pop_nuts2.ct]
    bio_nodal.index = pop_nuts2.index
    bio_nodal = bio_nodal.mul(pop_nuts2.fraction, axis=0).astype(float)

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


def build_eurostat(input_eurostat, countries, year, idees_rename):
    """
    Return multi-index for all countries' energy data in TWh/a.
    """
    df = {}
    countries = {idees_rename.get(country, country) for country in countries} - {"CH"}
    for country in countries:
        filename = (
            f"{input_eurostat}/{country}-Energy-balance-sheets-April-2023-edition.xlsb"
        )
        sheet = pd.read_excel(
            filename,
            engine="pyxlsb",
            sheet_name=str(year),
            skiprows=4,
            index_col=list(range(4)),
        )
        df[country] = sheet
    df = pd.concat(df, axis=0)

    # drop columns with all NaNs
    unnamed_cols = df.columns[df.columns.astype(str).str.startswith("Unnamed")]
    df.drop(unnamed_cols, axis=1, inplace=True)
    df.drop(year, axis=1, inplace=True)

    # make numeric values where possible
    df.replace("Z", 0, inplace=True)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.select_dtypes(include=[np.number])

    # write 'International aviation' to the 2nd level of the multiindex
    int_avia = df.index.get_level_values(2) == "International aviation"
    temp = df.loc[int_avia]
    temp.index = pd.MultiIndex.from_frame(
        temp.index.to_frame().fillna("International aviation")
    )
    df = pd.concat([temp, df.loc[~int_avia]])

    # Renaming some indices
    index_rename = {
        "Households": "Residential",
        "Commercial & public services": "Services",
        "Domestic navigation": "Domestic Navigation",
        "International maritime bunkers": "Bunkers",
    }
    columns_rename = {"Total": "Total all products", "UK": "GB"}
    df.rename(index=index_rename, columns=columns_rename, inplace=True)
    df.sort_index(inplace=True)
    df.index.names = [None] * len(df.index.names)

    # convert to MWh/a from ktoe/a
    df *= 11.63 * 1e3

    return df


def add_unsustainable_potentials(df):
    """
    Add unsustainable biomass potentials to the given dataframe.
    The difference between the data of JRC and Eurostat is assumed to be
    unsustainable biomass.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe with sustainable biomass potentials.
    unsustainable_biomass : str
        Path to the file with unsustainable biomass potentials.

    Returns
    -------
    pd.DataFrame
        The dataframe with added unsustainable biomass potentials.
    """
    if "GB" in snakemake.config["countries"]:
        latest_year = 2019
    else:
        latest_year = 2021
    idees_rename = {"GR": "EL", "GB": "UK"}
    df_unsustainable = (
        build_eurostat(
            countries=snakemake.config["countries"],
            year=max(
                min(latest_year, int(snakemake.wildcards.planning_horizons)), 1990
            ),
            input_eurostat=snakemake.input.eurostat,
            idees_rename=idees_rename,
        )
        .xs("Primary production", level=2)
        .droplevel([1, 2, 3])
    )

    df_unsustainable.index = df_unsustainable.index.str.strip()
    df_unsustainable = df_unsustainable.rename(
        {v: k for k, v in idees_rename.items()}, axis=0
    )

    bio_carriers = [
        "Primary solid biofuels",
        "Biogases",
        "Renewable municipal waste",
        "Pure biogasoline",
        "Blended biogasoline",
        "Pure biodiesels",
        "Blended biodiesels",
        "Pure bio jet kerosene",
        "Blended bio jet kerosene",
        "Other liquid biofuels",
    ]

    df_unsustainable = df_unsustainable[bio_carriers]

    df_wo_ch = df.drop(df.filter(regex="CH\d", axis=0).index)

    if snakemake.params["waste_incineration"]:
        df_unsustainable["Primary solid biofuels"] -= df_unsustainable[
            "Renewable municipal waste"
        ]

        df_wo_ch["unsustainable waste"] = (
            df_wo_ch.apply(
                lambda c: c.sum()
                / df_wo_ch.loc[df_wo_ch.index.str[:2] == c.name[:2]].sum().sum()
                * df_unsustainable.loc[c.name[:2], "Renewable municipal waste"],
                axis=1,
            )
            - df_wo_ch["waste"]
        ).clip(lower=0)

    df_wo_ch["unsustainable solid biomass"] = (
        df_wo_ch.apply(
            lambda c: c.sum()
            / df_wo_ch.loc[df_wo_ch.index.str[:2] == c.name[:2]].sum().sum()
            * df_unsustainable.loc[c.name[:2], "Primary solid biofuels"],
            axis=1,
        )
        - df_wo_ch["solid biomass"]
    ).clip(lower=0)

    df_wo_ch["unsustainable biogas"] = (
        df_wo_ch.apply(
            lambda c: c.sum()
            / df_wo_ch.loc[df_wo_ch.index.str[:2] == c.name[:2]].sum().sum()
            * df_unsustainable.loc[c.name[:2], "Biogases"],
            axis=1,
        )
        - df_wo_ch["biogas"]
    ).clip(lower=0)

    df_wo_ch["unsustainable bioliquids"] = df_wo_ch.apply(
        lambda c: c.sum()
        / df_wo_ch.loc[df_wo_ch.index.str[:2] == c.name[:2]].sum().sum()
        * df_unsustainable.filter(regex="gasoline|diesel|kerosene|liquid")
        .sum(axis=1)
        .loc[c.name[:2]],
        axis=1,
    )

    df = df.join(df_wo_ch.filter(like="unsustainable")).fillna(0)
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        import os

        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        snakemake = mock_snakemake(
            "build_biomass_potentials",
            simpl="",
            clusters="37",
            planning_horizons=2030,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

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
    df = df.T.groupby(grouper).sum().T

    df *= 1e6  # TWh/a to MWh/a
    df.index.name = "MWh/a"

    df = add_unsustainable_potentials(df)

    df.to_csv(snakemake.output.biomass_potentials)
