# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build biogas and solid biomass potentials per clustered model region from JRC ENSPRESO.

ENSPRESO gives sustainable biomass potentials per commodity at NUTS2 level for a
chosen scenario and year; years between the available decades are linearly
interpolated. Commodities only reported at country level are disaggregated to
NUTS2 by population. Potentials are mapped to model regions by area overlap and
grouped into the configured biomass classes. Unsustainable potentials are
estimated from the Eurostat primary production of solid biofuels, biogas and
bioliquids, allocated to regions by their share of the sustainable potential,
and scaled by year-dependent shares that phase them out. The potentials limit
the biomass supply in the sector-coupled network.

References
----------
- Ruiz et al. (2019), [ENSPRESO - an open, EU-28 wide, transparent and coherent database of wind, solar and biomass energy potentials](https://doi.org/10.1016/j.esr.2019.100379)
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd

from scripts._helpers import area, configure_logging, set_scenario_config

logger = logging.getLogger(__name__)
AVAILABLE_BIOMASS_YEARS = [2010, 2020, 2030, 2040, 2050]


def _calc_unsustainable_potential(
    df: pd.DataFrame,
    df_unsustainable: pd.DataFrame,
    share_unsus: float,
    resource_type: str,
) -> pd.Series:
    """
    Calculate the unsustainable biomass potential for a given resource type or
    regex.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe with sustainable biomass potentials.
    df_unsustainable : pd.DataFrame
        The dataframe with unsustainable biomass potentials.
    share_unsus : float
        The share of unsustainable biomass potential retained.
    resource_type : str or regex
        The resource type to calculate the unsustainable potential for.

    Returns
    -------
    pd.Series
        The unsustainable biomass potential for the given resource type or regex.
    """

    if "|" in resource_type:
        resource_potential = df_unsustainable.filter(regex=resource_type).sum(axis=1)
    else:
        resource_potential = df_unsustainable[resource_type]

    def _calculate_resource_allocation(c: pd.Series) -> float:
        country = c.name[:2]
        country_total = df.loc[df.index.str[:2] == country].sum().sum()
        if country_total == 0:
            return 0.0
        return c.sum() / country_total * resource_potential.loc[country]

    return (
        df.apply(_calculate_resource_allocation, axis=1).mul(share_unsus).clip(lower=0)
    )


def build_nuts_population_data(
    nuts3_population: str,
    swiss_cantons: str | list,
    swiss_population: str | list,
    year: int = 2013,
) -> pd.DataFrame:
    pop = pd.read_csv(
        nuts3_population,
        sep=r"\,| \t|\t",
        engine="python",
        na_values=[":"],
        index_col=1,
    )[str(year)]
    pop = pop.str.split().str[0].astype(float)

    if swiss_cantons:
        cantons = pd.read_csv(swiss_cantons)
        cantons = cantons.set_index(cantons.HASC.str[3:]).NUTS
        cantons = cantons.str.pad(5, side="right", fillchar="0")

        swiss = pd.read_excel(swiss_population, skiprows=3, index_col=0).loc[
            "Residents in 1000"
        ]
        swiss = swiss.rename(cantons).filter(like="CH")

        swiss = [swiss.groupby(swiss.index.str[:i]).sum() for i in range(2, 6)]
        pop = pd.concat([pop, pd.concat(swiss)])

    pop = pop.to_frame("total")

    # add missing manually
    pop["AL"] = 2778
    pop["BA"] = 3234
    pop["RS"] = 6664
    pop["ME"] = 617
    pop["XK"] = 1587

    pop["ct"] = pop.index.str[:2]

    return pop


def enspreso_biomass_potentials(
    enspreso_biomass: str, year: int = 2020, scenario: str = "ENS_Low"
) -> pd.DataFrame:
    """
    Loads the JRC ENSPRESO biomass potentials.

    Parameters
    ----------
    enspreso_biomass : str
        Path to the JRC ENSPRESO biomass Excel file.
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
        str(enspreso_biomass),
        sheet_name="Glossary",
        usecols="B:D",
        skiprows=1,
        index_col=0,
    )

    df = pd.read_excel(
        str(enspreso_biomass),
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

    return bio


def disaggregate_nuts0(bio: pd.DataFrame, pop: pd.DataFrame) -> pd.DataFrame:
    """
    Some commodities are only given on NUTS0 level. These are disaggregated
    here using the NUTS2 population as distribution key.

    Parameters
    ----------
    bio : pd.DataFrame
        from enspreso_biomass_potentials()
    pop : pd.DataFrame
        from build_nuts_population_data()

    Returns
    -------
    pd.DataFrame
    """

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


def build_nuts2_shapes(nuts2: str, country_shapes: str) -> gpd.GeoDataFrame:
    """
    - load NUTS2 geometries
    - add RS, AL, BA country shapes (not covered in NUTS 2013)
    - consistently name ME, MK
    """
    shapes = gpd.GeoDataFrame(gpd.read_file(nuts2).set_index("NUTS_ID").geometry)

    countries = gpd.read_file(country_shapes).set_index("name")
    missing_iso2 = countries.index.intersection(["AL", "RS", "XK", "BA"])
    missing = countries.loc[missing_iso2]

    shapes.rename(index={"ME00": "ME", "MK00": "MK"}, inplace=True)

    return pd.concat([shapes, missing])


def convert_nuts2_to_regions(
    bio_nuts2: gpd.GeoDataFrame, regions: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
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
        {"name", "area_nuts2", "geometry", "share", "country", "x", "y"}
    )
    overlay[adjust_cols] = overlay[adjust_cols].multiply(overlay["share"], axis=0)

    bio_regions = overlay.dissolve("name", aggfunc="sum")

    bio_regions.drop(["area_nuts2", "share"], axis=1, inplace=True)

    return bio_regions


def upstream_emission_factors(fn: str, year: int) -> pd.DataFrame:
    """
    Upstream emission factors of bioethanol and biodiesel per country
    (tCO2eq/MWh), linearly interpolated to the given year.
    """
    table = pd.read_csv(fn, index_col=0)
    factors = {}
    for fuel in ["bioethanol", "biodiesel"]:
        df = table.filter(regex=rf"^{fuel}_\d{{4}}$")
        years = df.columns.str[-4:].astype(int)
        factors[fuel] = df.apply(lambda row: np.interp(year, years, row), axis=1)
    return pd.DataFrame(factors)


def add_unsustainable_potentials(
    df: pd.DataFrame,
    input_eurostat: str,
    countries: list[str],
    investment_year: int,
    params: dict,
    upstream_emissions: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Add unsustainable biomass potentials to the given dataframe. The difference
    between the data of JRC and Eurostat is assumed to be unsustainable
    biomass.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe with sustainable biomass potentials.
    input_eurostat : str
        Path to the file with Eurostat biomass data.
    countries : list[str]
        Modelled countries.
    investment_year : int
        Planning horizon.
    params : dict
        Biomass configuration.

    Returns
    -------
    pd.DataFrame
        The dataframe with added unsustainable biomass potentials.
    """
    if "GB" in countries:
        latest_year = 2019
    else:
        latest_year = 2021
    idees_rename = {"GR": "EL", "GB": "UK"}
    year = max(min(latest_year, investment_year), 1990)  # noqa: F841
    df_unsustainable = (
        pd.read_csv(input_eurostat)
        .query("year == @year and nrg_bal == 'PPRD'")  # Primary production
        .set_index(["country", "siec"])
        .value.unstack("siec")
    )

    df_unsustainable = df_unsustainable.rename(
        {v: k for k, v in idees_rename.items()}, axis=0
    )

    bio_carriers = [
        "R5110-5150_W6000RI",  # Primary solid biofuels
        "R5300",  # Biogases
        "W6210",  # Renewable municipal waste
        "R5210P",  # Pure biogasoline
        "R5210B",  # Blended biogasoline
        "R5220P",  # Pure biodiesels
        "R5220B",  # Blended biodiesels
        "R5230P",  # Pure bio jet kerosene
        "R5230B",  # Blended bio jet kerosene
        "R5290",  # Other liquid biofuels
    ]

    df_unsustainable = df_unsustainable[bio_carriers]

    # Phase out unsustainable biomass potentials linearly from 2020 to 2035 while phasing in sustainable potentials
    share_unsus = params["share_unsustainable_use_retained"][investment_year]

    df_wo_ch = df.drop(df.filter(regex=r"CH\d*", axis=0).index)

    # Calculate unsustainable solid biomass
    df_wo_ch["unsustainable solid biomass"] = _calc_unsustainable_potential(
        df_wo_ch,
        df_unsustainable,
        share_unsus,
        "R5110-5150_W6000RI",  # Primary solid biofuels
    )

    # Calculate unsustainable biogas
    df_wo_ch["unsustainable biogas"] = _calc_unsustainable_potential(
        df_wo_ch,
        df_unsustainable,
        share_unsus,
        "R5300",  # Biogases
    )

    # Calculate unsustainable bioliquids
    df_wo_ch["unsustainable bioliquids"] = _calc_unsustainable_potential(
        df_wo_ch,
        df_unsustainable,
        share_unsus,
        resource_type="R5210|R5220|R5230|R5290",  # gasoline, diesel, kerosene, liquids
    )

    df_wo_ch["bioliquids phased out"] = _calc_unsustainable_potential(
        df_wo_ch,
        df_unsustainable,
        1 - share_unsus,
        resource_type="R5210|R5220|R5230|R5290",
    )

    if upstream_emissions is not None:
        # weight the factors with the country's shares of biogasoline and other bioliquids
        bioliquids = df_unsustainable.filter(regex="R5210|R5220|R5230|R5290")
        share_bioethanol = (
            bioliquids.filter(regex="R5210").sum(axis=1).div(bioliquids.sum(axis=1))
        )
        factor = upstream_emissions["biodiesel"] + share_bioethanol.fillna(0) * (
            upstream_emissions["bioethanol"] - upstream_emissions["biodiesel"]
        )
        df_wo_ch["unsustainable bioliquids upstream emissions"] = (
            df_wo_ch["unsustainable bioliquids"]
            * factor.reindex(df_wo_ch.index.str[:2]).values
        )

    share_sus = params["share_sustainable_potential_available"][investment_year]
    df.loc[df_wo_ch.index] *= share_sus

    df = (
        df.join(df_wo_ch.filter(regex="unsustainable|phased out"))
        .fillna(0)
        .infer_objects(copy=False)
    )

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_biomass_potentials",
            horizon=2050,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    overnight = snakemake.params.foresight == "overnight"
    params = snakemake.params.biomass
    investment_year = int(snakemake.wildcards.horizon)
    year = params["year"] if overnight else investment_year
    scenario = params["scenario"]
    enspreso_fn = snakemake.input.enspreso_biomass

    if year > 2050:
        logger.info("No biomass potentials for years after 2050, using 2050.")
        max_year = max(AVAILABLE_BIOMASS_YEARS)
        enspreso = enspreso_biomass_potentials(enspreso_fn, max_year, scenario)

    elif year not in AVAILABLE_BIOMASS_YEARS:
        before = int(np.floor(year / 10) * 10)
        after = int(np.ceil(year / 10) * 10)
        logger.info(
            f"No biomass potentials for {year}, interpolating linearly between {before} and {after}."
        )

        enspreso_before = enspreso_biomass_potentials(enspreso_fn, before, scenario)
        enspreso_after = enspreso_biomass_potentials(enspreso_fn, after, scenario)

        fraction = (year - before) / (after - before)

        enspreso = enspreso_before + fraction * (enspreso_after - enspreso_before)

    else:
        logger.info(f"Using biomass potentials for {year}.")
        enspreso = enspreso_biomass_potentials(enspreso_fn, year, scenario)

    pop = build_nuts_population_data(
        snakemake.input.nuts3_population,
        snakemake.input.swiss_cantons,
        snakemake.input.swiss_population,
    )
    enspreso = disaggregate_nuts0(enspreso, pop)

    nuts2 = build_nuts2_shapes(snakemake.input.nuts2, snakemake.input.country_shapes)

    df_nuts2 = gpd.GeoDataFrame(nuts2.geometry).join(enspreso)

    regions = gpd.read_file(snakemake.input.onshore_regions)

    df = convert_nuts2_to_regions(df_nuts2, regions)

    df.to_csv(snakemake.output.biomass_potentials_all)

    grouper = {v: k for k, vv in params["classes"].items() for v in vv}
    df = df.T.groupby(grouper).sum().T

    input_eurostat = snakemake.input.eurostat
    upstream_emissions = None
    if snakemake.input.unsustainable_biomass_upstream:
        upstream_emissions = upstream_emission_factors(
            snakemake.input.unsustainable_biomass_upstream, investment_year
        )

    df = add_unsustainable_potentials(
        df,
        input_eurostat,
        snakemake.params.countries,
        investment_year,
        params,
        upstream_emissions,
    )

    df *= 1e6  # TWh/a to MWh/a
    df.index.name = "MWh/a"

    df.to_csv(snakemake.output.biomass_potentials)
