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

import copy
import logging

import geopandas as gpd
import numpy as np
import pandas as pd

from scripts._helpers import area, configure_logging, load_costs, set_scenario_config
from scripts.build_perennials_yields_eurostat_average import crop_energy_per_tonne

logger = logging.getLogger(__name__)
AVAILABLE_BIOMASS_YEARS = [2010, 2020, 2030, 2040, 2050]
# ENSPRESO year of the 1G biofuel crops behind today's unsustainable
# bioliquids (Eurostat 2019/2021 taken as 2020)
ENSPRESO_1G_YEAR = 2020

# 1G biofuel crop groups and their target biomass class names.
# In the default config these sit in "not included".
# When sector: perennials=True, resolve_biomass_classes() moves them automatically.
ONE_G_BIOFUEL_CLASSES = {
    "Bioethanol barley, wheat, grain maize, oats, other cereals and rye": "biofuels_1G_bioethanol_cereals",
    "Sugar from sugar beet": "biofuels_1G_bioethanol_sugar",
    "Rape seed": "biofuels_1G_biodiesel",
    "Sunflower, soya seed ": "biofuels_1G_biodiesel",
}


def resolve_biomass_classes(classes: dict, perennials_enabled: bool) -> dict:
    """
    Move the 1G biofuel crop groups from 'not included' into their
    biofuels_1G_* classes when perennials are enabled.

    Cases:
      perennials=False  → return classes unchanged.
      perennials=True, all groups in 'not included' → move them (normal case).
      perennials=True, no groups in 'not included' → raise AssertionError.
      perennials=True, some groups elsewhere → move available ones, warn about rest.
    """
    if not perennials_enabled:
        return classes
    classes = copy.deepcopy(classes)
    not_incl = classes.get("not included", [])
    in_ni = [g for g in ONE_G_BIOFUEL_CLASSES if g in not_incl]
    elsewhere = [g for g in ONE_G_BIOFUEL_CLASSES if g not in not_incl]
    if not in_ni:
        raise AssertionError(
            "sector: perennials=true but no 1G-biofuel groups are in "
            "biomass: classes: 'not included'. If all 1G-biofuels are used "
            "there is no land available for perennials."
        )
    if elsewhere:
        logger.warning(
            "perennials: 1G groups already allocated outside 'not included' "
            f"— skipped (running perennials only for available groups): {elsewhere}"
        )
    for g in in_ni:
        classes["not included"].remove(g)
        classes.setdefault(ONE_G_BIOFUEL_CLASSES[g], []).append(g)
    return classes


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
    bio_nuts2: gpd.GeoDataFrame, regions: gpd.GeoDataFrame, extensive: bool = True
) -> gpd.GeoDataFrame | pd.DataFrame:
    """
    Converts values given in NUTS2 to PyPSA-Eur regions based on the overlay of
    both GeoDataFrames.

    Parameters
    ----------
    bio_nuts2 : gpd.GeoDataFrame
        Values indexed by NUTS2 shapes, e.g. JRC ENSPRESO biomass potentials.
    regions : gpd.GeoDataFrame
        PyPSA-Eur clustered onshore regions
    extensive : bool, default True
        If True, the values are totals (e.g. MWh/a) and are split in proportion
        to the share of the NUTS2 area inside each region. If False, the values
        are rates (e.g. yields in MWh/ha/a) and each region gets the average of
        the overlapping NUTS2 values weighted by the overlap area; NUTS2 regions
        without a value are left out of the average.

    Returns
    -------
    gpd.GeoDataFrame or pd.DataFrame
    """
    if not extensive:
        overlay = gpd.overlay(regions, bio_nuts2, keep_geom_type=True)
        overlap = area(overlay)
        cols = overlay.columns.difference({"name", "geometry"})
        values = overlay[cols].astype(float)
        weights = values.notna().mul(overlap, axis=0)
        weighted = values.mul(overlap, axis=0).groupby(overlay["name"]).sum(min_count=1)
        return weighted.div(weights.groupby(overlay["name"]).sum()).reindex(
            regions["name"].values
        )

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


def add_unsustainable_potentials(
    df: pd.DataFrame,
    input_eurostat: str,
    countries: list[str],
    investment_year: int,
    params: dict,
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

    share_sus = params["share_sustainable_potential_available"][investment_year]
    df.loc[df_wo_ch.index] *= share_sus

    df = (
        df.join(df_wo_ch.filter(like="unsustainable"))
        .fillna(0)
        .infer_objects(copy=False)
    )

    return df


def eurostat_primary_production(
    input_eurostat: str, countries: list[str], investment_year: int
) -> tuple[pd.DataFrame, int]:
    """
    Eurostat primary production (TWh) per country and fuel (siec), for the
    latest year used for the unsustainable biomass (2019 with GB, else 2021).

    Returns
    -------
    pd.DataFrame
        Primary production, countries as rows and siec codes as columns.
    int
        The Eurostat year.
    """
    if "GB" in countries:
        latest_year = 2019
    else:
        latest_year = 2021
    idees_rename = {"GR": "EL", "GB": "UK"}
    year = max(min(latest_year, investment_year), 1990)
    df = (
        pd.read_csv(input_eurostat)
        .query("year == @year and nrg_bal == 'PPRD'")  # Primary production
        .set_index(["country", "siec"])
        .value.unstack("siec")
    )
    return df.rename({v: k for k, v in idees_rename.items()}, axis=0), year


def biofuels_1G_from_bioliquids(
    enspreso_1G: pd.DataFrame,
    eurostat: pd.DataFrame,
    crop_energy: dict,
    efficiency: dict,
    bioliquids: dict,
    share_freed: float,
    land_multiplier: float = 1.0,
) -> pd.DataFrame:
    """
    1st-generation (1G) biofuel crop potentials (TWh of crop/y) behind the
    unsustainable bioliquids freed by their phase-out, the cropland available
    for perennialisation.

    The ENSPRESO 1G potentials (2020) give where the crops grow (NUTS2) and
    the crop and fuel mix. They are scaled with one factor, the same for all
    crops and regions, so that their fuel energy matches the Eurostat
    bioliquids of the model (biogasoline, biodiesels, bio jet kerosene and
    other liquid biofuels). Crop energy is converted to fuel energy with the
    crop-to-fuel efficiencies of technology-data (t fuel per t crop) and the
    heating values of the fuels and of the crops (ENSPRESO basis).

    Parameters
    ----------
    enspreso_1G : pd.DataFrame
        ENSPRESO 1G potentials per region (TWh of crop/y), columns
        biofuels_1G_bioethanol_cereals, _sugar and biofuels_1G_biodiesel.
    eurostat : pd.DataFrame
        Eurostat primary production (TWh) per country and fuel.
    crop_energy : dict
        MWh per tonne of crop, per biofuels_1G class.
    efficiency : dict
        Tonne of fuel per tonne of crop, per biofuels_1G class.
    bioliquids : dict
        ``biomass: unsustainable_bioliquids`` config (fuel heating values).
    share_freed : float
        Share of the unsustainable bioliquids phased out.
    land_multiplier : float, default 1.0
        Multiplier of the cropland, set by the user to project the land
        available.
    """
    fuel = {
        "biofuels_1G_bioethanol_cereals": "bioethanol",
        "biofuels_1G_bioethanol_sugar": "bioethanol",
        "biofuels_1G_biodiesel": "biodiesel",
    }

    # Switzerland has no unsustainable bioliquids (no Eurostat data)
    enspreso_1G = enspreso_1G[list(fuel)].copy()
    enspreso_1G.loc[enspreso_1G.index.str.startswith("CH")] = 0.0

    # MWh of crop per MWh of fuel
    crop_per_fuel = pd.Series(
        {
            k: crop_energy[k] / (efficiency[k] * bioliquids[f"MWh_{f}_per_t{f}"])
            for k, f in fuel.items()
        }
    )

    countries = enspreso_1G.index.str[:2].unique()
    eurostat_fuel = (
        eurostat.reindex(countries).filter(regex="R5210|R5220|R5230|R5290").sum().sum()
    )
    enspreso_fuel = (enspreso_1G / crop_per_fuel).sum().sum()
    scale = eurostat_fuel / enspreso_fuel
    logger.info(
        f"1G biofuel crops: Eurostat bioliquids {eurostat_fuel:.1f} TWh fuel, "
        f"ENSPRESO {enspreso_fuel:.1f} TWh fuel, scaled by {scale:.3f}."
    )

    return enspreso_1G * scale * share_freed * land_multiplier


def enspreso_potentials_for_year(
    enspreso_fn: str, year: int, scenario: str
) -> pd.DataFrame:
    """
    ENSPRESO potentials for a year, interpolated linearly between the
    available years and taken at 2050 after 2050.
    """
    if year > 2050:
        logger.info("No biomass potentials for years after 2050, using 2050.")
        max_year = max(AVAILABLE_BIOMASS_YEARS)
        return enspreso_biomass_potentials(enspreso_fn, max_year, scenario)

    if year not in AVAILABLE_BIOMASS_YEARS:
        before = int(np.floor(year / 10) * 10)
        after = int(np.ceil(year / 10) * 10)
        logger.info(
            f"No biomass potentials for {year}, interpolating linearly between {before} and {after}."
        )

        enspreso_before = enspreso_biomass_potentials(enspreso_fn, before, scenario)
        enspreso_after = enspreso_biomass_potentials(enspreso_fn, after, scenario)

        fraction = (year - before) / (after - before)

        return enspreso_before + fraction * (enspreso_after - enspreso_before)

    logger.info(f"Using biomass potentials for {year}.")
    return enspreso_biomass_potentials(enspreso_fn, year, scenario)


def potentials_per_region(
    enspreso: pd.DataFrame,
    pop: pd.DataFrame,
    nuts2: gpd.GeoDataFrame,
    regions: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Map ENSPRESO NUTS2 potentials onto the clustered regions.
    """
    enspreso = disaggregate_nuts0(enspreso, pop)
    df_nuts2 = gpd.GeoDataFrame(nuts2.geometry).join(enspreso)
    return convert_nuts2_to_regions(df_nuts2, regions)


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

    enspreso = enspreso_potentials_for_year(enspreso_fn, year, scenario)

    pop = build_nuts_population_data(
        snakemake.input.nuts3_population,
        snakemake.input.swiss_cantons,
        snakemake.input.swiss_population,
    )

    nuts2 = build_nuts2_shapes(snakemake.input.nuts2, snakemake.input.country_shapes)

    regions = gpd.read_file(snakemake.input.onshore_regions)

    df = potentials_per_region(enspreso, pop, nuts2, regions)

    df.to_csv(snakemake.output.biomass_potentials_all)

    perennials = snakemake.config["sector"].get("perennials", False)
    classes = resolve_biomass_classes(params["classes"], perennials)
    grouper = {v: k for k, vv in classes.items() for v in vv}
    df = df.T.groupby(grouper).sum().T

    input_eurostat = snakemake.input.eurostat
    df = add_unsustainable_potentials(
        df, input_eurostat, snakemake.params.countries, investment_year, params
    )

    df *= 1e6  # TWh/a to MWh/a
    df.index.name = "MWh/a"

    if perennials:
        # the 1G biofuel cropland freed by the phase-out of the unsustainable
        # bioliquids: ENSPRESO 2020, scaled to the same Eurostat data
        eurostat, _ = eurostat_primary_production(
            input_eurostat, snakemake.params.countries, investment_year
        )
        enspreso_1G = potentials_per_region(
            enspreso_biomass_potentials(enspreso_fn, ENSPRESO_1G_YEAR, scenario),
            pop,
            nuts2,
            regions,
        )
        enspreso_1G = enspreso_1G.T.groupby(grouper).sum().T
        costs = load_costs(snakemake.input.costs)
        efficiency = {
            "biofuels_1G_bioethanol_cereals": costs.at[
                "ethanol from wheat", "efficiency"
            ],
            "biofuels_1G_bioethanol_sugar": costs.at[
                "ethanol from sugar beet", "efficiency"
            ],
            "biofuels_1G_biodiesel": costs.at["biodiesel from rapeseed", "efficiency"],
        }
        crop_energy = crop_energy_per_tonne(snakemake.input.crop_heating_values)
        crop_energy = {
            "biofuels_1G_bioethanol_cereals": crop_energy["MINBIOCRP11"],
            "biofuels_1G_bioethanol_sugar": crop_energy["MINBIOCRP21"],
            "biofuels_1G_biodiesel": crop_energy["MINBIORPS1"],
        }
        share_freed = 1 - params["share_unsustainable_use_retained"][investment_year]
        biofuels_1G = biofuels_1G_from_bioliquids(
            enspreso_1G,
            eurostat,
            crop_energy,
            efficiency,
            params["unsustainable_bioliquids"],
            share_freed,
            snakemake.config["perennials"]["land_multiplier"],
        )
        df[biofuels_1G.columns] = biofuels_1G * 1e6  # TWh/a to MWh/a

    df.to_csv(snakemake.output.biomass_potentials)
