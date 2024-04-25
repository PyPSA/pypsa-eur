# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build total energy demands per country using JRC IDEES, eurostat, and EEA data.
"""

import logging
import multiprocessing as mp
import os
from functools import partial

import country_converter as coco
import geopandas as gpd
import numpy as np
import pandas as pd
from _helpers import configure_logging, mute_print, set_scenario_config
from tqdm import tqdm

cc = coco.CountryConverter()
logger = logging.getLogger(__name__)
idx = pd.IndexSlice


def cartesian(s1, s2):
    """
    Cartesian product of two pd.Series.
    """
    return pd.DataFrame(np.outer(s1, s2), index=s1.index, columns=s2.index)


def reverse(dictionary):
    """
    Reverses a keys and values of a dictionary.
    """
    return {v: k for k, v in dictionary.items()}


idees_rename = {"GR": "EL", "GB": "UK"}

eu28 = cc.EU28as("ISO2").ISO2.tolist()

eu28_eea = eu28.copy()
eu28_eea.remove("GB")
eu28_eea.append("UK")


to_ipcc = {
    "electricity": "1.A.1.a - Public Electricity and Heat Production",
    "residential non-elec": "1.A.4.b - Residential",
    "services non-elec": "1.A.4.a - Commercial/Institutional",
    "rail non-elec": "1.A.3.c - Railways",
    "road non-elec": "1.A.3.b - Road Transportation",
    "domestic navigation": "1.A.3.d - Domestic Navigation",
    "international navigation": "1.D.1.b - International Navigation",
    "domestic aviation": "1.A.3.a - Domestic Aviation",
    "international aviation": "1.D.1.a - International Aviation",
    "total energy": "1 - Energy",
    "industrial processes": "2 - Industrial Processes and Product Use",
    "agriculture": "3 - Agriculture",
    "agriculture, forestry and fishing": "1.A.4.c - Agriculture/Forestry/Fishing",
    "LULUCF": "4 - Land Use, Land-Use Change and Forestry",
    "waste management": "5 - Waste management",
    "other": "6 - Other Sector",
    "indirect": "ind_CO2 - Indirect CO2",
    "total wL": "Total (with LULUCF)",
    "total woL": "Total (without LULUCF)",
}


def eurostat_per_country(input_eurostat, country):
    filename = (
        f"{input_eurostat}/{country}-Energy-balance-sheets-April-2023-edition.xlsb"
    )
    sheet = pd.read_excel(
        filename,
        engine="pyxlsb",
        sheet_name=None,
        skiprows=4,
        index_col=list(range(4)),
    )
    sheet.pop("Cover")
    return pd.concat(sheet)


def build_eurostat(input_eurostat, countries, nprocesses=1, disable_progressbar=False):
    """
    Return multi-index for all countries' energy data in TWh/a.
    """
    countries = {idees_rename.get(country, country) for country in countries} - {"CH"}

    func = partial(eurostat_per_country, input_eurostat)
    tqdm_kwargs = dict(
        ascii=False,
        unit=" country",
        total=len(countries),
        desc="Build from eurostat database",
        disable=disable_progressbar,
    )
    with mute_print():
        with mp.Pool(processes=nprocesses) as pool:
            dfs = list(tqdm(pool.imap(func, countries), **tqdm_kwargs))

    index_names = ["country", "year", "lvl1", "lvl2", "lvl3", "lvl4"]
    df = pd.concat(dfs, keys=countries, names=index_names)
    df.index = df.index.set_levels(df.index.levels[1].astype(int), level=1)

    # drop columns with all NaNs
    unnamed_cols = df.columns[df.columns.astype(str).str.startswith("Unnamed")]
    df.drop(unnamed_cols, axis=1, inplace=True)
    df.drop(list(range(1990, 2022)), axis=1, inplace=True, errors="ignore")

    # make numeric values where possible
    df.replace("Z", 0, inplace=True)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.select_dtypes(include=[np.number])

    # write 'International aviation' to the lower level of the multiindex
    int_avia = df.index.get_level_values(3) == "International aviation"
    temp = df.loc[int_avia]
    temp.index = pd.MultiIndex.from_frame(
        temp.index.to_frame().fillna("International aviation")
    )
    df = pd.concat([temp, df.loc[~int_avia]])

    # Fill in missing data on "Domestic aviation" for each country.
    domestic_avia = df.index.get_level_values(4) == "Domestic aviation"
    for country in countries:
        slicer = idx[country, :, :, :, "Domestic aviation"]
        # For the Total and Fossil energy columns, fill in zeros with
        # the closest non-zero value in the year index.
        for col in ["Total", "Fossil energy"]:
            df.loc[slicer, col] = (
                df.loc[slicer, col].replace(0.0, np.nan).ffill().bfill()
            )

    # Renaming some indices
    index_rename = {
        "Households": "Residential",
        "Commercial & public services": "Services",
        "Domestic navigation": "Domestic Navigation",
        "International maritime bunkers": "Bunkers",
        "UK": "GB",
    }
    columns_rename = {"Total": "Total all products"}
    df.rename(index=index_rename, columns=columns_rename, inplace=True)
    df.sort_index(inplace=True)

    # convert to TWh/a from ktoe/a
    df *= 11.63 / 1e3

    return df


def build_swiss():
    """
    Return a pd.DataFrame of Swiss energy data in TWh/a.
    """
    fn = snakemake.input.swiss

    df = pd.read_csv(fn, index_col=[0, 1])

    df.columns = df.columns.astype(int)

    df.columns.name = "year"

    df = df.stack().unstack("item")

    df.columns.name = None

    # convert PJ/a to TWh/a
    df /= 3.6

    return df


def idees_per_country(ct, base_dir):
    ct_idees = idees_rename.get(ct, ct)
    fn_residential = f"{base_dir}/JRC-IDEES-2015_Residential_{ct_idees}.xlsx"
    fn_tertiary = f"{base_dir}/JRC-IDEES-2015_Tertiary_{ct_idees}.xlsx"
    fn_transport = f"{base_dir}/JRC-IDEES-2015_Transport_{ct_idees}.xlsx"

    ct_totals = {}

    # residential

    df = pd.read_excel(fn_residential, "RES_hh_fec", index_col=0)

    rows = ["Advanced electric heating", "Conventional electric heating"]
    ct_totals["electricity residential space"] = df.loc[rows].sum()
    ct_totals["total residential space"] = df.loc["Space heating"]
    ct_totals["total residential water"] = df.loc["Water heating"]

    assert df.index[23] == "Electricity"
    ct_totals["electricity residential water"] = df.iloc[23]

    ct_totals["total residential cooking"] = df.loc["Cooking"]

    assert df.index[30] == "Electricity"
    ct_totals["electricity residential cooking"] = df.iloc[30]

    df = pd.read_excel(fn_residential, "RES_summary", index_col=0)

    row = "Energy consumption by fuel - Eurostat structure (ktoe)"
    ct_totals["total residential"] = df.loc[row]

    assert df.index[47] == "Electricity"
    ct_totals["electricity residential"] = df.iloc[47]

    assert df.index[46] == "Derived heat"
    ct_totals["derived heat residential"] = df.iloc[46]

    assert df.index[50] == "Thermal uses"
    ct_totals["thermal uses residential"] = df.iloc[50]

    # services

    df = pd.read_excel(fn_tertiary, "SER_hh_fec", index_col=0)

    ct_totals["total services space"] = df.loc["Space heating"]

    rows = ["Advanced electric heating", "Conventional electric heating"]
    ct_totals["electricity services space"] = df.loc[rows].sum()

    ct_totals["total services water"] = df.loc["Hot water"]

    assert df.index[24] == "Electricity"
    ct_totals["electricity services water"] = df.iloc[24]

    ct_totals["total services cooking"] = df.loc["Catering"]

    assert df.index[31] == "Electricity"
    ct_totals["electricity services cooking"] = df.iloc[31]

    df = pd.read_excel(fn_tertiary, "SER_summary", index_col=0)

    row = "Energy consumption by fuel - Eurostat structure (ktoe)"
    ct_totals["total services"] = df.loc[row]

    assert df.index[50] == "Electricity"
    ct_totals["electricity services"] = df.iloc[50]

    assert df.index[49] == "Derived heat"
    ct_totals["derived heat services"] = df.iloc[49]

    assert df.index[53] == "Thermal uses"
    ct_totals["thermal uses services"] = df.iloc[53]

    # agriculture, forestry and fishing

    start = "Detailed split of energy consumption (ktoe)"
    end = "Market shares of energy uses (%)"

    df = pd.read_excel(fn_tertiary, "AGR_fec", index_col=0).loc[start:end]

    rows = [
        "Lighting",
        "Ventilation",
        "Specific electricity uses",
        "Pumping devices (electric)",
    ]
    ct_totals["total agriculture electricity"] = df.loc[rows].sum()

    rows = ["Specific heat uses", "Low enthalpy heat"]
    ct_totals["total agriculture heat"] = df.loc[rows].sum()

    rows = [
        "Motor drives",
        "Farming machine drives (diesel oil incl. biofuels)",
        "Pumping devices (diesel oil incl. biofuels)",
    ]
    ct_totals["total agriculture machinery"] = df.loc[rows].sum()

    row = "Agriculture, forestry and fishing"
    ct_totals["total agriculture"] = df.loc[row]

    # transport

    df = pd.read_excel(fn_transport, "TrRoad_ene", index_col=0)

    ct_totals["total road"] = df.loc["by fuel (EUROSTAT DATA)"]

    ct_totals["electricity road"] = df.loc["Electricity"]

    ct_totals["total two-wheel"] = df.loc["Powered 2-wheelers (Gasoline)"]

    assert df.index[19] == "Passenger cars"
    ct_totals["total passenger cars"] = df.iloc[19]

    assert df.index[30] == "Battery electric vehicles"
    ct_totals["electricity passenger cars"] = df.iloc[30]

    assert df.index[31] == "Motor coaches, buses and trolley buses"
    ct_totals["total other road passenger"] = df.iloc[31]

    assert df.index[39] == "Battery electric vehicles"
    ct_totals["electricity other road passenger"] = df.iloc[39]

    assert df.index[41] == "Light duty vehicles"
    ct_totals["total light duty road freight"] = df.iloc[41]

    assert df.index[49] == "Battery electric vehicles"
    ct_totals["electricity light duty road freight"] = df.iloc[49]

    row = "Heavy duty vehicles (Diesel oil incl. biofuels)"
    ct_totals["total heavy duty road freight"] = df.loc[row]

    assert df.index[61] == "Passenger cars"
    ct_totals["passenger car efficiency"] = df.iloc[61]

    df = pd.read_excel(fn_transport, "TrRail_ene", index_col=0)

    ct_totals["total rail"] = df.loc["by fuel (EUROSTAT DATA)"]

    ct_totals["electricity rail"] = df.loc["Electricity"]

    assert df.index[15] == "Passenger transport"
    ct_totals["total rail passenger"] = df.iloc[15]

    assert df.index[16] == "Metro and tram, urban light rail"
    assert df.index[19] == "Electric"
    assert df.index[20] == "High speed passenger trains"
    ct_totals["electricity rail passenger"] = df.iloc[[16, 19, 20]].sum()

    assert df.index[21] == "Freight transport"
    ct_totals["total rail freight"] = df.iloc[21]

    assert df.index[23] == "Electric"
    ct_totals["electricity rail freight"] = df.iloc[23]

    df = pd.read_excel(fn_transport, "TrAvia_ene", index_col=0)

    assert df.index[6] == "Passenger transport"
    ct_totals["total aviation passenger"] = df.iloc[6]

    assert df.index[10] == "Freight transport"
    ct_totals["total aviation freight"] = df.iloc[10]

    assert df.index[7] == "Domestic"
    ct_totals["total domestic aviation passenger"] = df.iloc[7]

    assert df.index[8] == "International - Intra-EU"
    assert df.index[9] == "International - Extra-EU"
    ct_totals["total international aviation passenger"] = df.iloc[[8, 9]].sum()

    assert df.index[11] == "Domestic and International - Intra-EU"
    ct_totals["total domestic aviation freight"] = df.iloc[11]

    assert df.index[12] == "International - Extra-EU"
    ct_totals["total international aviation freight"] = df.iloc[12]

    ct_totals["total domestic aviation"] = (
        ct_totals["total domestic aviation freight"]
        + ct_totals["total domestic aviation passenger"]
    )

    ct_totals["total international aviation"] = (
        ct_totals["total international aviation freight"]
        + ct_totals["total international aviation passenger"]
    )

    df = pd.read_excel(fn_transport, "TrNavi_ene", index_col=0)

    # coastal and inland
    ct_totals["total domestic navigation"] = df.loc["by fuel (EUROSTAT DATA)"]

    df = pd.read_excel(fn_transport, "TrRoad_act", index_col=0)

    assert df.index[85] == "Passenger cars"
    ct_totals["passenger cars"] = df.iloc[85]

    return pd.DataFrame(ct_totals)


def build_idees(countries):
    nprocesses = snakemake.threads
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    func = partial(idees_per_country, base_dir=snakemake.input.idees)
    tqdm_kwargs = dict(
        ascii=False,
        unit=" country",
        total=len(countries),
        desc="Build from IDEES database",
        disable=disable_progress,
    )
    with mute_print():
        with mp.Pool(processes=nprocesses) as pool:
            totals_list = list(tqdm(pool.imap(func, countries), **tqdm_kwargs))

    totals = pd.concat(
        totals_list,
        keys=countries,
        names=["country", "year"],
    )

    # convert ktoe to TWh
    exclude = totals.columns.str.fullmatch("passenger cars")
    totals.loc[:, ~exclude] *= 11.63 / 1e3

    # convert TWh/100km to kWh/km
    totals.loc[:, "passenger car efficiency"] *= 10

    return totals


def build_energy_totals(countries, eurostat, swiss, idees):
    eurostat_fuels = {"electricity": "Electricity", "total": "Total all products"}
    eurostat_countries = eurostat.index.levels[0]
    eurostat_years = eurostat.index.levels[1]

    to_drop = ["passenger cars", "passenger car efficiency"]
    new_index = pd.MultiIndex.from_product(
        [countries, eurostat_years], names=["country", "year"]
    )

    df = idees.reindex(new_index).drop(to_drop, axis=1)

    in_eurostat = df.index.levels[0].intersection(eurostat_countries)

    # add international navigation

    slicer = idx[in_eurostat, :, :, "Bunkers", :]
    fill_values = eurostat.loc[slicer, "Total all products"].groupby(level=[0, 1]).sum()
    df.loc[in_eurostat, "total international navigation"] = fill_values

    # add swiss energy data

    df = pd.concat([df.drop("CH", errors="ignore"), swiss]).sort_index()

    # get values for missing countries based on Eurostat EnergyBalances
    # divide cooking/space/water according to averages in EU28

    uses = ["space", "cooking", "water"]

    to_fill = df.index[
        df["total residential"].isna()
        & df.index.get_level_values("country").isin(eurostat_countries)
    ]
    c = to_fill.get_level_values("country")
    y = to_fill.get_level_values("year")

    for sector in ["residential", "services", "road", "rail"]:
        eurostat_sector = sector.capitalize()

        # fuel use

        for fuel in ["electricity", "total"]:
            slicer = idx[c, y, :, :, eurostat_sector]
            fill_values = (
                eurostat.loc[slicer, eurostat_fuels[fuel]].groupby(level=[0, 1]).sum()
            )
            df.loc[to_fill, f"{fuel} {sector}"] = fill_values

    for sector in ["residential", "services"]:
        # electric use

        for use in uses:
            fuel_use = df[f"electricity {sector} {use}"]
            fuel = df[f"electricity {sector}"]
            avg = fuel_use.div(fuel).mean()
            logger.debug(
                f"{sector}: average fraction of electricity for {use} is {avg:.3f}"
            )
            df.loc[to_fill, f"electricity {sector} {use}"] = (
                avg * df.loc[to_fill, f"electricity {sector}"]
            )

        # non-electric use

        for use in uses:
            nonelectric_use = (
                df[f"total {sector} {use}"] - df[f"electricity {sector} {use}"]
            )
            nonelectric = df[f"total {sector}"] - df[f"electricity {sector}"]
            avg = nonelectric_use.div(nonelectric).mean()
            logger.debug(
                f"{sector}: average fraction of non-electric for {use} is {avg:.3f}"
            )
            electric_use = df.loc[to_fill, f"electricity {sector} {use}"]
            nonelectric = (
                df.loc[to_fill, f"total {sector}"]
                - df.loc[to_fill, f"electricity {sector}"]
            )
            df.loc[to_fill, f"total {sector} {use}"] = electric_use + avg * nonelectric

    # Fix Norway space and water heating fractions
    # http://www.ssb.no/en/energi-og-industri/statistikker/husenergi/hvert-3-aar/2014-07-14
    # The main heating source for about 73 per cent of the households is based on electricity
    # => 26% is non-electric

    if "NO" in df.index:
        elec_fraction = 0.73

        no_norway = df.drop("NO")

        for sector in ["residential", "services"]:
            # assume non-electric is heating
            nonelectric = (
                df.loc["NO", f"total {sector}"] - df.loc["NO", f"electricity {sector}"]
            )
            total_heating = nonelectric / (1 - elec_fraction)

            for use in uses:
                nonelectric_use = (
                    no_norway[f"total {sector} {use}"]
                    - no_norway[f"electricity {sector} {use}"]
                )
                nonelectric = (
                    no_norway[f"total {sector}"] - no_norway[f"electricity {sector}"]
                )
                fraction = nonelectric_use.div(nonelectric).mean()
                df.loc["NO", f"total {sector} {use}"] = (
                    total_heating * fraction
                ).values
                df.loc["NO", f"electricity {sector} {use}"] = (
                    total_heating * fraction * elec_fraction
                ).values

    # Missing aviation

    slicer = idx[c, y, :, :, "Domestic aviation"]
    fill_values = eurostat.loc[slicer, "Total all products"].groupby(level=[0, 1]).sum()
    df.loc[to_fill, "total domestic aviation"] = fill_values

    slicer = idx[c, y, :, :, "International aviation"]
    fill_values = eurostat.loc[slicer, "Total all products"].groupby(level=[0, 1]).sum()
    df.loc[to_fill, "total international aviation"] = fill_values

    # missing domestic navigation

    slicer = idx[c, y, :, :, "Domestic Navigation"]
    fill_values = eurostat.loc[slicer, "Total all products"].groupby(level=[0, 1]).sum()
    df.loc[to_fill, "total domestic navigation"] = fill_values

    # split road traffic for non-IDEES
    missing = df.index[df["total passenger cars"].isna()]
    for fuel in ["total", "electricity"]:
        selection = [
            f"{fuel} passenger cars",
            f"{fuel} other road passenger",
            f"{fuel} light duty road freight",
        ]
        if fuel == "total":
            selection.extend([f"{fuel} two-wheel", f"{fuel} heavy duty road freight"])
        road = df[selection].sum()
        road_fraction = road / road.sum()
        fill_values = cartesian(df.loc[missing, f"{fuel} road"], road_fraction)
        df.loc[missing, road_fraction.index] = fill_values

    # split rail traffic for non-IDEES
    missing = df.index[df["total rail passenger"].isna()]
    for fuel in ["total", "electricity"]:
        selection = [f"{fuel} rail passenger", f"{fuel} rail freight"]
        rail = df[selection].sum()
        rail_fraction = rail / rail.sum()
        fill_values = cartesian(df.loc[missing, f"{fuel} rail"], rail_fraction)
        df.loc[missing, rail_fraction.index] = fill_values

    # split aviation traffic for non-IDEES
    missing = df.index[df["total domestic aviation passenger"].isna()]
    for destination in ["domestic", "international"]:
        selection = [
            f"total {destination} aviation passenger",
            f"total {destination} aviation freight",
        ]
        aviation = df[selection].sum()
        aviation_fraction = aviation / aviation.sum()
        fill_values = cartesian(
            df.loc[missing, f"total {destination} aviation"], aviation_fraction
        )
        df.loc[missing, aviation_fraction.index] = fill_values

    for purpose in ["passenger", "freight"]:
        attrs = [
            f"total domestic aviation {purpose}",
            f"total international aviation {purpose}",
        ]
        df.loc[missing, f"total aviation {purpose}"] = df.loc[missing, attrs].sum(
            axis=1
        )

    if "BA" in df.index:
        # fill missing data for BA (services and road energy data)
        # proportional to RS with ratio of total residential demand
        mean_BA = df.loc["BA"].loc[2014:2021, "total residential"].mean()
        mean_RS = df.loc["RS"].loc[2014:2021, "total residential"].mean()
        ratio = mean_BA / mean_RS
        df.loc["BA"] = df.loc["BA"].replace(0.0, np.nan).values
        df.loc["BA"] = df.loc["BA"].combine_first(ratio * df.loc["RS"]).values

    return df


def build_district_heat_share(countries, idees):
    # district heating share
    district_heat = idees[["derived heat residential", "derived heat services"]].sum(
        axis=1
    )
    total_heat = idees[["thermal uses residential", "thermal uses services"]].sum(
        axis=1
    )

    district_heat_share = district_heat / total_heat

    district_heat_share = district_heat_share.reindex(countries, level="country")

    # Missing district heating share
    dh_share = (
        pd.read_csv(snakemake.input.district_heat_share, index_col=0, usecols=[0, 1])
        .div(100)
        .squeeze()
    )
    # make conservative assumption and take minimum from both data sets
    district_heat_share = pd.concat(
        [district_heat_share, dh_share.reindex_like(district_heat_share)], axis=1
    ).min(axis=1)

    district_heat_share.name = "district heat share"

    # restrict to available years
    district_heat_share = (
        district_heat_share.unstack().dropna(how="all", axis=1).ffill(axis=1)
    )

    return district_heat_share


def build_eea_co2(input_co2, year=1990, emissions_scope="CO2"):
    # https://www.eea.europa.eu/data-and-maps/data/national-emissions-reported-to-the-unfccc-and-to-the-eu-greenhouse-gas-monitoring-mechanism-16
    # downloaded 201228 (modified by EEA last on 201221)
    df = pd.read_csv(input_co2, encoding="latin-1", low_memory=False)

    df.replace(dict(Year="1985-1987"), 1986, inplace=True)
    df.Year = df.Year.astype(int)
    index_col = ["Country_code", "Pollutant_name", "Year", "Sector_name"]
    df = df.set_index(index_col).sort_index()

    cts = ["CH", "EUA", "NO"] + eu28_eea

    slicer = idx[cts, emissions_scope, year, to_ipcc.values()]
    emissions = (
        df.loc[slicer, "emissions"]
        .unstack("Sector_name")
        .rename(columns=reverse(to_ipcc))
        .droplevel([1, 2])
    )

    emissions.rename(index={"EUA": "EU28", "UK": "GB"}, inplace=True)

    to_subtract = [
        "electricity",
        "services non-elec",
        "residential non-elec",
        "road non-elec",
        "rail non-elec",
        "domestic aviation",
        "international aviation",
        "domestic navigation",
        "international navigation",
        "agriculture, forestry and fishing",
    ]
    emissions["industrial non-elec"] = emissions["total energy"] - emissions[
        to_subtract
    ].sum(axis=1)

    emissions["agriculture"] += emissions["agriculture, forestry and fishing"]

    to_drop = [
        "total energy",
        "total wL",
        "total woL",
        "agriculture, forestry and fishing",
    ]
    emissions.drop(columns=to_drop, inplace=True)

    # convert from Gg to Mt
    return emissions / 1e3


def build_eurostat_co2(eurostat, year=1990):
    eurostat_year = eurostat.xs(year, level="year")

    specific_emissions = pd.Series(index=eurostat.columns, dtype=float)

    # emissions in tCO2_equiv per MWh_th
    specific_emissions["Solid fuels"] = 0.36  # Approximates coal
    specific_emissions["Oil (total)"] = 0.285  # Average of distillate and residue
    specific_emissions["Gas"] = 0.2  # For natural gas

    # oil values from https://www.eia.gov/tools/faqs/faq.cfm?id=74&t=11
    # Distillate oil (No. 2)  0.276
    # Residual oil (No. 6)  0.298
    # https://www.eia.gov/electricity/annual/html/epa_a_03.html

    return eurostat_year.multiply(specific_emissions).sum(axis=1)


def build_co2_totals(countries, eea_co2, eurostat_co2):
    co2 = eea_co2.reindex(countries)

    for ct in pd.Index(countries).intersection(["BA", "RS", "AL", "ME", "MK"]):
        mappings = {
            "electricity": (ct, "+", "Electricity & heat generation", np.nan),
            "residential non-elec": (ct, "+", "+", "Residential"),
            "services non-elec": (ct, "+", "+", "Services"),
            "road non-elec": (ct, "+", "+", "Road"),
            "rail non-elec": (ct, "+", "+", "Rail"),
            "domestic navigation": (ct, "+", "+", "Domestic Navigation"),
            "international navigation": (ct, "-", "Bunkers"),
            "domestic aviation": (ct, "+", "+", "Domestic aviation"),
            "international aviation": (ct, "-", "International aviation"),
            # does not include industrial process emissions or fuel processing/refining
            "industrial non-elec": (ct, "+", "Industry sector"),
            # does not include non-energy emissions
            "agriculture": (eurostat_co2.index.get_level_values(0) == ct)
            & eurostat_co2.index.isin(["Agriculture & forestry", "Fishing"], level=3),
        }

        for i, mi in mappings.items():
            co2.at[ct, i] = eurostat_co2.loc[mi].sum()

    return co2


def build_transport_data(countries, population, idees):
    # first collect number of cars

    transport_data = pd.DataFrame(idees["passenger cars"])

    countries_without_ch = set(countries) - {"CH"}
    new_index = pd.MultiIndex.from_product(
        [countries_without_ch, transport_data.index.levels[1]],
        names=["country", "year"],
    )

    transport_data = transport_data.reindex(index=new_index)

    # https://www.bfs.admin.ch/bfs/en/home/statistics/mobility-transport/transport-infrastructure-vehicles/vehicles/road-vehicles-stock-level-motorisation.html
    if "CH" in countries:
        fn = snakemake.input.swiss_transport
        swiss_cars = pd.read_csv(fn, index_col=0).loc[2000:2015, ["passenger cars"]]

        swiss_cars.index = pd.MultiIndex.from_product(
            [["CH"], swiss_cars.index], names=["country", "year"]
        )

        transport_data = pd.concat([transport_data, swiss_cars]).sort_index()

    transport_data.rename(columns={"passenger cars": "number cars"}, inplace=True)

    missing = transport_data.index[transport_data["number cars"].isna()]
    if not missing.empty:
        logger.info(
            f"Missing data on cars from:\n{list(missing)}\nFilling gaps with averaged data."
        )

        cars_pp = transport_data["number cars"] / population

        fill_values = {
            year: cars_pp.mean() * population for year in transport_data.index.levels[1]
        }
        fill_values = pd.DataFrame(fill_values).stack()
        fill_values = pd.DataFrame(fill_values, columns=["number cars"])
        fill_values.index.names = ["country", "year"]
        fill_values = fill_values.reindex(transport_data.index)

        transport_data = transport_data.combine_first(fill_values)

    # collect average fuel efficiency in kWh/km

    transport_data["average fuel efficiency"] = idees["passenger car efficiency"]

    missing = transport_data.index[transport_data["average fuel efficiency"].isna()]
    if not missing.empty:
        logger.info(
            f"Missing data on fuel efficiency from:\n{list(missing)}\nFilling gaps with averaged data."
        )

        fill_values = transport_data["average fuel efficiency"].mean()
        transport_data.loc[missing, "average fuel efficiency"] = fill_values

    return transport_data


def rescale_idees_from_eurostat(
    idees_countries,
    energy,
    eurostat,
):
    """
    Takes JRC IDEES data from 2015 and rescales it by the ratio of the eurostat
    data and the 2015 eurostat data.

    missing data: ['passenger car efficiency', 'passenger cars']
    """
    main_cols = ["Total all products", "Electricity"]
    # read in the eurostat data for 2015
    eurostat_2015 = eurostat.xs(2015, level="year")[main_cols]
    # calculate the ratio of the two data sets
    ratio = eurostat[main_cols] / eurostat_2015
    ratio = ratio.droplevel([2, 5])
    cols_rename = {"Total all products": "total", "Electricity": "ele"}
    index_rename = {v: k for k, v in idees_rename.items()}
    ratio.rename(columns=cols_rename, index=index_rename, inplace=True)

    mappings = {
        "Residential": {
            "total": [
                "total residential space",
                "total residential water",
                "total residential cooking",
                "total residential",
                "derived heat residential",
                "thermal uses residential",
            ],
            "elec": [
                "electricity residential space",
                "electricity residential water",
                "electricity residential cooking",
                "electricity residential",
            ],
        },
        "Services": {
            "total": [
                "total services space",
                "total services water",
                "total services cooking",
                "total services",
                "derived heat services",
                "thermal uses services",
            ],
            "elec": [
                "electricity services space",
                "electricity services water",
                "electricity services cooking",
                "electricity services",
            ],
        },
        "Agriculture & forestry": {
            "total": [
                "total agriculture heat",
                "total agriculture machinery",
                "total agriculture",
            ],
            "elec": [
                "total agriculture electricity",
            ],
        },
        "Road": {
            "total": [
                "total road",
                "total passenger cars",
                "total other road passenger",
                "total light duty road freight",
            ],
            "elec": [
                "electricity road",
                "electricity passenger cars",
                "electricity other road passenger",
                "electricity light duty road freight",
            ],
        },
        "Rail": {
            "total": [
                "total rail",
                "total rail passenger",
                "total rail freight",
            ],
            "elec": [
                "electricity rail",
                "electricity rail passenger",
                "electricity rail freight",
            ],
        },
    }

    avia_inter = [
        "total aviation passenger",
        "total aviation freight",
        "total international aviation passenger",
        "total international aviation freight",
        "total international aviation",
    ]
    avia_domestic = [
        "total domestic aviation passenger",
        "total domestic aviation freight",
        "total domestic aviation",
    ]
    navigation = [
        "total domestic navigation",
    ]

    for country in idees_countries:
        filling_years = [(2015, slice(2016, 2021)), (2000, slice(1990, 1999))]

        for source_year, target_years in filling_years:

            slicer_source = idx[country, source_year, :, :]
            slicer_target = idx[country, target_years, :, :]

            for sector, mapping in mappings.items():
                sector_ratio = ratio.loc[
                    (country, slice(None), slice(None), sector)
                ].droplevel("lvl2")

                energy.loc[slicer_target, mapping["total"]] = cartesian(
                    sector_ratio.loc[target_years, "total"],
                    energy.loc[slicer_source, mapping["total"]].squeeze(axis=0),
                ).values
                energy.loc[slicer_target, mapping["elec"]] = cartesian(
                    sector_ratio.loc[target_years, "ele"],
                    energy.loc[slicer_source, mapping["elec"]].squeeze(axis=0),
                ).values

            level_drops = ["country", "lvl2", "lvl3"]

            slicer = idx[country, :, :, "Domestic aviation"]
            avi_d = ratio.loc[slicer, "total"].droplevel(level_drops)

            slicer = idx[country, :, :, "International aviation"]
            avi_i = ratio.loc[slicer, "total"].droplevel(level_drops)

            slicer = idx[country, :, :, "Domestic Navigation"]
            nav = ratio.loc[slicer, "total"].droplevel(level_drops)

            energy.loc[slicer_target, avia_inter] = cartesian(
                avi_i.loc[target_years],
                energy.loc[slicer_source, avia_inter].squeeze(axis=0),
            ).values

            energy.loc[slicer_target, avia_domestic] = cartesian(
                avi_d.loc[target_years],
                energy.loc[slicer_source, avia_domestic].squeeze(axis=0),
            ).values

            energy.loc[slicer_target, navigation] = cartesian(
                nav.loc[target_years],
                energy.loc[slicer_source, navigation].squeeze(axis=0),
            ).values

    return energy


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_energy_totals")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params.energy

    nuts3 = gpd.read_file(snakemake.input.nuts3_shapes).set_index("index")
    population = nuts3["pop"].groupby(nuts3.country).sum()

    countries = snakemake.params.countries
    idees_countries = pd.Index(countries).intersection(eu28)

    input_eurostat = snakemake.input.eurostat
    eurostat = build_eurostat(
        input_eurostat,
        countries,
        nprocesses=snakemake.threads,
        disable_progressbar=snakemake.config["run"].get("disable_progressbar", False),
    )
    swiss = build_swiss()
    idees = build_idees(idees_countries)

    energy = build_energy_totals(countries, eurostat, swiss, idees)

    # Data from IDEES only exists from 2000-2015.
    logger.info("Extrapolate IDEES data based on eurostat for years 2015-2021.")
    energy = rescale_idees_from_eurostat(idees_countries, energy, eurostat)

    energy.to_csv(snakemake.output.energy_name)

    # use rescaled idees data to calculate district heat share
    district_heat_share = build_district_heat_share(
        countries, energy.loc[idees_countries]
    )
    district_heat_share.to_csv(snakemake.output.district_heat_share)

    base_year_emissions = params["base_emissions_year"]
    emissions_scope = snakemake.params.energy["emissions"]
    eea_co2 = build_eea_co2(snakemake.input.co2, base_year_emissions, emissions_scope)
    eurostat_co2 = build_eurostat_co2(eurostat, base_year_emissions)

    co2 = build_co2_totals(countries, eea_co2, eurostat_co2)
    co2.to_csv(snakemake.output.co2_name)

    transport = build_transport_data(countries, population, idees)
    transport.to_csv(snakemake.output.transport_name)
