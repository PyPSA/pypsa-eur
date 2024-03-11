# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build industrial energy demand per country.
"""

import multiprocessing as mp
from functools import partial

import country_converter as coco
import pandas as pd
from _helpers import set_scenario_config
from tqdm import tqdm

cc = coco.CountryConverter()

ktoe_to_twh = 0.011630

# name in JRC-IDEES Energy Balances
sector_sheets = {
    "Integrated steelworks": "cisb",
    "Electric arc": "cise",
    "Alumina production": "cnfa",
    "Aluminium - primary production": "cnfp",
    "Aluminium - secondary production": "cnfs",
    "Other non-ferrous metals": "cnfo",
    "Basic chemicals": "cbch",
    "Other chemicals": "coch",
    "Pharmaceutical products etc.": "cpha",
    "Basic chemicals feedstock": "cpch",
    "Cement": "ccem",
    "Ceramics & other NMM": "ccer",
    "Glass production": "cgla",
    "Pulp production": "cpul",
    "Paper production": "cpap",
    "Printing and media reproduction": "cprp",
    "Food, beverages and tobacco": "cfbt",
    "Transport Equipment": "ctre",
    "Machinery Equipment": "cmae",
    "Textiles and leather": "ctel",
    "Wood and wood products": "cwwp",
    "Mining and quarrying": "cmiq",
    "Construction": "ccon",
    "Non-specified": "cnsi",
}


fuels = {
    "All Products": "all",
    "Solid Fuels": "solid",
    "Total petroleum products (without biofuels)": "liquid",
    "Gases": "gas",
    "Nuclear heat": "heat",
    "Derived heat": "heat",
    "Biomass and Renewable wastes": "biomass",
    "Wastes (non-renewable)": "waste",
    "Electricity": "electricity",
}

eu28 = cc.EU28as("ISO2").ISO2.tolist()

jrc_names = {"GR": "EL", "GB": "UK"}


def industrial_energy_demand_per_country(country, year, jrc_dir):
    jrc_country = jrc_names.get(country, country)
    fn = f"{jrc_dir}/JRC-IDEES-2015_EnergyBalance_{jrc_country}.xlsx"

    sheets = list(sector_sheets.values())
    df_dict = pd.read_excel(fn, sheet_name=sheets, index_col=0)

    def get_subsector_data(sheet):
        df = df_dict[sheet][year].groupby(fuels).sum()

        df["hydrogen"] = 0.0

        df["other"] = df["all"] - df.loc[df.index != "all"].sum()

        return df

    df = pd.concat(
        {sub: get_subsector_data(sheet) for sub, sheet in sector_sheets.items()}, axis=1
    )

    sel = ["Mining and quarrying", "Construction", "Non-specified"]
    df["Other Industrial Sectors"] = df[sel].sum(axis=1)
    df["Basic chemicals"] += df["Basic chemicals feedstock"]

    df.drop(columns=sel + ["Basic chemicals feedstock"], index="all", inplace=True)

    df *= ktoe_to_twh

    return df


def separate_basic_chemicals(demand, production):

    ammonia = pd.DataFrame(
        {
            "hydrogen": production["Ammonia"] * params["MWh_H2_per_tNH3_electrolysis"],
            "electricity": production["Ammonia"]
            * params["MWh_elec_per_tNH3_electrolysis"],
        }
    ).T
    chlorine = pd.DataFrame(
        {
            "hydrogen": production["Chlorine"] * params["MWh_H2_per_tCl"],
            "electricity": production["Chlorine"] * params["MWh_elec_per_tCl"],
        }
    ).T
    methanol = pd.DataFrame(
        {
            "gas": production["Methanol"] * params["MWh_CH4_per_tMeOH"],
            "electricity": production["Methanol"] * params["MWh_elec_per_tMeOH"],
        }
    ).T

    demand["Ammonia"] = ammonia.unstack().reindex(index=demand.index, fill_value=0.0)
    demand["Chlorine"] = chlorine.unstack().reindex(index=demand.index, fill_value=0.0)
    demand["Methanol"] = methanol.unstack().reindex(index=demand.index, fill_value=0.0)

    demand["HVC"] = (
        demand["Basic chemicals"]
        - demand["Ammonia"]
        - demand["Methanol"]
        - demand["Chlorine"]
    )

    demand.drop(columns="Basic chemicals", inplace=True)

    demand["HVC"].clip(lower=0, inplace=True)

    return demand


def add_non_eu28_industrial_energy_demand(countries, demand, production):
    non_eu28 = countries.difference(eu28)
    if non_eu28.empty:
        return demand

    eu28_production = production.loc[countries.intersection(eu28)].sum()
    eu28_energy = demand.groupby(level=1).sum()
    eu28_averages = eu28_energy / eu28_production

    demand_non_eu28 = pd.concat(
        {k: v * eu28_averages for k, v in production.loc[non_eu28].iterrows()}
    )

    return pd.concat([demand, demand_non_eu28])


def industrial_energy_demand(countries, year):
    nprocesses = snakemake.threads
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    func = partial(
        industrial_energy_demand_per_country, year=year, jrc_dir=snakemake.input.jrc
    )
    tqdm_kwargs = dict(
        ascii=False,
        unit=" country",
        total=len(countries),
        desc="Build industrial energy demand",
        disable=disable_progress,
    )
    with mp.Pool(processes=nprocesses) as pool:
        demand_l = list(tqdm(pool.imap(func, countries), **tqdm_kwargs))

    return pd.concat(demand_l, keys=countries)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industrial_energy_demand_per_country_today")
    set_scenario_config(snakemake)

    params = snakemake.params.industry
    year = params.get("reference_year", 2015)
    countries = pd.Index(snakemake.params.countries)

    demand = industrial_energy_demand(countries.intersection(eu28), year)

    # output in MtMaterial/a
    production = (
        pd.read_csv(snakemake.input.industrial_production_per_country, index_col=0)
        / 1e3
    )

    demand = separate_basic_chemicals(demand, production)

    demand = add_non_eu28_industrial_energy_demand(countries, demand, production)

    # for format compatibility
    demand = demand.stack(dropna=False).unstack(level=[0, 2])

    # style and annotation
    demand.index.name = "TWh/a"
    demand.sort_index(axis=1, inplace=True)

    fn = snakemake.output.industrial_energy_demand_per_country_today
    demand.to_csv(fn)
