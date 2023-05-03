# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build industrial production per country.
"""

import logging
from functools import partial

logger = logging.getLogger(__name__)

import multiprocessing as mp

import country_converter as coco
import numpy as np
import pandas as pd
from _helpers import mute_print
from tqdm import tqdm

cc = coco.CountryConverter()

tj_to_ktoe = 0.0238845
ktoe_to_twh = 0.01163

sub_sheet_name_dict = {
    "Iron and steel": "ISI",
    "Chemicals Industry": "CHI",
    "Non-metallic mineral products": "NMM",
    "Pulp, paper and printing": "PPA",
    "Food, beverages and tobacco": "FBT",
    "Non Ferrous Metals": "NFM",
    "Transport Equipment": "TRE",
    "Machinery Equipment": "MAE",
    "Textiles and leather": "TEL",
    "Wood and wood products": "WWP",
    "Other Industrial Sectors": "OIS",
}

eu28 = cc.EU28as("ISO2").ISO2.values

jrc_names = {"GR": "EL", "GB": "UK"}

sect2sub = {
    "Iron and steel": ["Electric arc", "Integrated steelworks"],
    "Chemicals Industry": [
        "Basic chemicals",
        "Other chemicals",
        "Pharmaceutical products etc.",
    ],
    "Non-metallic mineral products": [
        "Cement",
        "Ceramics & other NMM",
        "Glass production",
    ],
    "Pulp, paper and printing": [
        "Pulp production",
        "Paper production",
        "Printing and media reproduction",
    ],
    "Food, beverages and tobacco": ["Food, beverages and tobacco"],
    "Non Ferrous Metals": [
        "Alumina production",
        "Aluminium - primary production",
        "Aluminium - secondary production",
        "Other non-ferrous metals",
    ],
    "Transport Equipment": ["Transport Equipment"],
    "Machinery Equipment": ["Machinery Equipment"],
    "Textiles and leather": ["Textiles and leather"],
    "Wood and wood products": ["Wood and wood products"],
    "Other Industrial Sectors": ["Other Industrial Sectors"],
}

sub2sect = {v: k for k, vv in sect2sub.items() for v in vv}

fields = {
    "Electric arc": "Electric arc",
    "Integrated steelworks": "Integrated steelworks",
    "Basic chemicals": "Basic chemicals (kt ethylene eq.)",
    "Other chemicals": "Other chemicals (kt ethylene eq.)",
    "Pharmaceutical products etc.": "Pharmaceutical products etc. (kt ethylene eq.)",
    "Cement": "Cement (kt)",
    "Ceramics & other NMM": "Ceramics & other NMM (kt bricks eq.)",
    "Glass production": "Glass production  (kt)",
    "Pulp production": "Pulp production (kt)",
    "Paper production": "Paper production  (kt)",
    "Printing and media reproduction": "Printing and media reproduction (kt paper eq.)",
    "Food, beverages and tobacco": "Physical output (index)",
    "Alumina production": "Alumina production (kt)",
    "Aluminium - primary production": "Aluminium - primary production",
    "Aluminium - secondary production": "Aluminium - secondary production",
    "Other non-ferrous metals": "Other non-ferrous metals (kt lead eq.)",
    "Transport Equipment": "Physical output (index)",
    "Machinery Equipment": "Physical output (index)",
    "Textiles and leather": "Physical output (index)",
    "Wood and wood products": "Physical output (index)",
    "Other Industrial Sectors": "Physical output (index)",
}

eb_names = {
    "NO": "Norway",
    "AL": "Albania",
    "BA": "Bosnia and Herzegovina",
    "MK": "FYR of Macedonia",
    "GE": "Georgia",
    "IS": "Iceland",
    "KO": "Kosovo",
    "MD": "Moldova",
    "ME": "Montenegro",
    "RS": "Serbia",
    "UA": "Ukraine",
    "TR": "Turkey",
}

eb_sectors = {
    "Iron & steel industry": "Iron and steel",
    "Chemical and Petrochemical industry": "Chemicals Industry",
    "Non-ferrous metal industry": "Non-metallic mineral products",
    "Paper, Pulp and Print": "Pulp, paper and printing",
    "Food and Tabacco": "Food, beverages and tobacco",
    "Non-metallic Minerals (Glass, pottery & building mat. Industry)": "Non Ferrous Metals",
    "Transport Equipment": "Transport Equipment",
    "Machinery": "Machinery Equipment",
    "Textile and Leather": "Textiles and leather",
    "Wood and Wood Products": "Wood and wood products",
    "Non-specified (Industry)": "Other Industrial Sectors",
}

# TODO: this should go in a csv in `data`
# Annual energy consumption in Switzerland by sector in 2015 (in TJ)
# From: Energieverbrauch in der Industrie und im Dienstleistungssektor, Der Bundesrat
# http://www.bfe.admin.ch/themen/00526/00541/00543/index.html?lang=de&dossier_id=00775
e_switzerland = pd.Series(
    {
        "Iron and steel": 7889.0,
        "Chemicals Industry": 26871.0,
        "Non-metallic mineral products": 15513.0 + 3820.0,
        "Pulp, paper and printing": 12004.0,
        "Food, beverages and tobacco": 17728.0,
        "Non Ferrous Metals": 3037.0,
        "Transport Equipment": 14993.0,
        "Machinery Equipment": 4724.0,
        "Textiles and leather": 1742.0,
        "Wood and wood products": 0.0,
        "Other Industrial Sectors": 10825.0,
        "current electricity": 53760.0,
    }
)


def find_physical_output(df):
    start = np.where(df.index.str.contains("Physical output", na=""))[0][0]
    empty_row = np.where(df.index.isnull())[0]
    end = empty_row[np.argmax(empty_row > start)]
    return slice(start, end)


def get_energy_ratio(country, eurostat_dir, jrc_dir, year):
    if country == "CH":
        e_country = e_switzerland * tj_to_ktoe
    else:
        # estimate physical output, energy consumption in the sector and country
        fn = f"{eurostat_dir}/{eb_names[country]}.XLSX"
        with mute_print():
            df = pd.read_excel(
                fn, sheet_name="2016", index_col=2, header=0, skiprows=1
            ).squeeze("columns")
        e_country = df.loc[eb_sectors.keys(), "Total all products"].rename(eb_sectors)

    fn = f"{jrc_dir}/JRC-IDEES-2015_Industry_EU28.xlsx"

    with mute_print():
        df = pd.read_excel(fn, sheet_name="Ind_Summary", index_col=0, header=0).squeeze(
            "columns"
        )

    assert df.index[48] == "by sector"
    year_i = df.columns.get_loc(year)
    e_eu28 = df.iloc[49:76, year_i]
    e_eu28.index = e_eu28.index.str.lstrip()

    e_ratio = e_country / e_eu28

    return pd.Series({k: e_ratio[v] for k, v in sub2sect.items()})


def industry_production_per_country(country, year, eurostat_dir, jrc_dir):
    def get_sector_data(sector, country):
        jrc_country = jrc_names.get(country, country)
        fn = f"{jrc_dir}/JRC-IDEES-2015_Industry_{jrc_country}.xlsx"
        sheet = sub_sheet_name_dict[sector]
        with mute_print():
            df = pd.read_excel(fn, sheet_name=sheet, index_col=0, header=0).squeeze(
                "columns"
            )

        year_i = df.columns.get_loc(year)
        df = df.iloc[find_physical_output(df), year_i]

        df = df.loc[map(fields.get, sect2sub[sector])]
        df.index = sect2sub[sector]

        return df

    ct = "EU28" if country not in eu28 else country
    demand = pd.concat([get_sector_data(s, ct) for s in sect2sub])

    if country not in eu28:
        demand *= get_energy_ratio(country, eurostat_dir, jrc_dir, year)

    demand.name = country

    return demand


def industry_production(countries, year, eurostat_dir, jrc_dir):
    nprocesses = snakemake.threads
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    func = partial(
        industry_production_per_country,
        year=year,
        eurostat_dir=eurostat_dir,
        jrc_dir=jrc_dir,
    )
    tqdm_kwargs = dict(
        ascii=False,
        unit=" country",
        total=len(countries),
        desc="Build industry production",
        disable=disable_progress,
    )
    with mp.Pool(processes=nprocesses) as pool:
        demand_l = list(tqdm(pool.imap(func, countries), **tqdm_kwargs))

    demand = pd.concat(demand_l, axis=1).T

    demand.index.name = "kton/a"

    return demand


def separate_basic_chemicals(demand, year):
    """
    Separate basic chemicals into ammonia, chlorine, methanol and HVC.
    """
    ammonia = pd.read_csv(snakemake.input.ammonia_production, index_col=0)

    there = ammonia.index.intersection(demand.index)
    missing = demand.index.symmetric_difference(there)

    logger.info(f"Following countries have no ammonia demand: {missing.tolist()}")

    demand["Ammonia"] = 0.0

    demand.loc[there, "Ammonia"] = ammonia.loc[there, str(year)]

    demand["Basic chemicals"] -= demand["Ammonia"]

    # EE, HR and LT got negative demand through subtraction - poor data
    demand["Basic chemicals"].clip(lower=0.0, inplace=True)

    # assume HVC, methanol, chlorine production proportional to non-ammonia basic chemicals
    distribution_key = demand["Basic chemicals"] / demand["Basic chemicals"].sum()
    demand["HVC"] = config["HVC_production_today"] * 1e3 * distribution_key
    demand["Chlorine"] = config["chlorine_production_today"] * 1e3 * distribution_key
    demand["Methanol"] = config["methanol_production_today"] * 1e3 * distribution_key

    demand.drop(columns=["Basic chemicals"], inplace=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industrial_production_per_country")

    logging.basicConfig(level=snakemake.config["logging"]["level"])

    countries = snakemake.config["countries"]

    year = snakemake.config["industry"]["reference_year"]

    config = snakemake.config["industry"]

    jrc_dir = snakemake.input.jrc
    eurostat_dir = snakemake.input.eurostat

    demand = industry_production(countries, year, eurostat_dir, jrc_dir)

    separate_basic_chemicals(demand, year)

    fn = snakemake.output.industrial_production_per_country
    demand.to_csv(fn, float_format="%.2f")
