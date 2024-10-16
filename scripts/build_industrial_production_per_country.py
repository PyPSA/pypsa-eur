# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This rule builds the historical industrial production per country.

Relevant Settings
-----------------

.. code:: yaml

    countries:
..

Inputs
-------
- ``resources/ammonia_production.csv``
- ``data/bundle-sector/jrc-idees-2021``
- ``data/eurostat``

Outputs
-------

- ``resources/industrial_production_per_country.csv``

Description
-------

The industrial production is taken from the `JRC-IDEES <https://joint-research-centre.ec.europa.eu/potencia-policy-oriented-tool-energy-and-climate-change-impact-assessment/jrc-idees_en)>`.
This dataset provides detailed information about the consumption of energy for various processes.
If the country is not part of the EU28, the energy consumption in the industrial sectors is taken from the `Eurostat <https://ec.europa.eu/eurostat/de/data/database>` dataset. The industrial production is calculated for the year specified in the config["industry"]["reference_year"].

The ammonia production is provided by the rule `build_ammonia_production <https://pypsa-eur.readthedocs.io/en/latest/sector.html#module-build_ammonia_production>`. Since Switzerland is not part of the EU28 nor reported by eurostat, the energy consumption in the industrial sectors is taken from the `BFE <https://pubdb.bfe.admin.ch/de/publication/download/11817> dataset.
After the industrial production is calculated, the basic chemicals are separated into ammonia, chlorine, methanol and HVC. The production of these chemicals is assumed to be proportional to the production of basic chemicals without ammonia.

The following subcategories [kton/a] are considered:
- Electric arc
- Integrated steelworks
- Other chemicals
- Pharmaceutical products etc.
- Cement
- Ceramics & other NMM
- Glass production
- Pulp production
- Paper production
- Printing and media reproduction
- Food, beverages and tobacco
- Alumina production
- Aluminium - primary production
- Aluminium - secondary production
- Other non-ferrous metals
- Transport equipment
- Machinery equipment
- Textiles and leather
- Wood and wood products
- Other industrial sectors
- Ammonia
- HVC
- Chlorine
- Methanol
"""

import logging
import multiprocessing as mp
from functools import partial

import country_converter as coco
import numpy as np
import pandas as pd
from _helpers import configure_logging, mute_print, set_scenario_config
from tqdm import tqdm

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()

tj_to_ktoe = 0.0238845
ktoe_to_twh = 0.01163

sub_sheet_name_dict = {
    "Iron and steel": "ISI",
    "Chemical industry": "CHI",
    "Non-metallic mineral products": "NMM",
    "Pulp, paper and printing": "PPA",
    "Food, beverages and tobacco": "FBT",
    "Non Ferrous Metals": "NFM",
    "Transport equipment": "TRE",
    "Machinery equipment": "MAE",
    "Textiles and leather": "TEL",
    "Wood and wood products": "WWP",
    "Other industrial sectors": "OIS",
}

eu27 = cc.EU27as("ISO2").ISO2.values

jrc_names = {"GR": "EL", "GB": "UK"}

sect2sub = {
    "Iron and steel": ["Electric arc", "Integrated steelworks"],
    "Chemical industry": [
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
    "Transport equipment": ["Transport equipment"],
    "Machinery equipment": ["Machinery equipment"],
    "Textiles and leather": ["Textiles and leather"],
    "Wood and wood products": ["Wood and wood products"],
    "Other industrial sectors": ["Other industrial sectors"],
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
    "Transport equipment": "Physical output (index)",
    "Machinery equipment": "Physical output (index)",
    "Textiles and leather": "Physical output (index)",
    "Wood and wood products": "Physical output (index)",
    "Other industrial sectors": "Physical output (index)",
}

eb_sectors = {
    "Iron & steel": "Iron and steel",
    "Chemical & petrochemical": "Chemical industry",
    "Non-ferrous metals": "Non-metallic mineral products",
    "Paper, pulp & printing": "Pulp, paper and printing",
    "Food, beverages & tobacco": "Food, beverages and tobacco",
    "Non-metallic minerals": "Non Ferrous Metals",
    "Transport equipment": "Transport equipment",
    "Machinery": "Machinery equipment",
    "Textile & leather": "Textiles and leather",
    "Wood & wood products": "Wood and wood products",
    "Not elsewhere specified (industry)": "Other industrial sectors",
}


ch_mapping = {
    "Nahrung": "Food, beverages and tobacco",
    "Textil / Leder": "Textiles and leather",
    "Papier / Druck": "Pulp, paper and printing",
    "Chemie / Pharma": "Chemical industry",
    "Zement / Beton": "Non-metallic mineral products",
    "Andere NE-Mineralien": "Other non-ferrous metals",
    "Metall / Eisen": "Iron and steel",
    "NE-Metall": "Non Ferrous Metals",
    "Metall / Geräte": "Transport equipment",
    "Maschinen": "Machinery equipment",
    "Andere Industrien": "Other industrial sectors",
}


def find_physical_output(df):
    start = np.where(df.index.str.contains("Physical output", na=""))[0][0]
    empty_row = np.where(df.index.isnull())[0]
    end = empty_row[np.argmax(empty_row > start)]
    return slice(start, end)


def get_energy_ratio(country, eurostat_dir, jrc_dir, year, snakemake):
    if country == "CH":
        # data ranges between 2014-2023
        e_country = pd.read_csv(
            snakemake.input.ch_industrial_production, index_col=0
        ).dropna()
        e_country = e_country.rename(index=ch_mapping).groupby(level=0).sum()
        e_country = e_country[str(min(2019, year))]
        e_country *= tj_to_ktoe
    else:
        ct_eurostat = country.replace("GB", "UK")
        # estimate physical output, energy consumption in the sector and country
        fn = f"{eurostat_dir}/{ct_eurostat}-Energy-balance-sheets-April-2023-edition.xlsb"
        df = pd.read_excel(
            fn,
            sheet_name=str(min(2019, year)),
            index_col=2,
            header=0,
            skiprows=4,
        )
        e_country = df.loc[eb_sectors.keys(), "Total"].rename(eb_sectors)

    fn = f"{jrc_dir}/EU27/JRC-IDEES-2021_Industry_EU27.xlsx"

    with mute_print():
        df = pd.read_excel(fn, sheet_name="Ind_Summary", index_col=0, header=0).squeeze(
            "columns"
        )

    assert df.index[49] == "by sector"
    year_i = df.columns.get_loc(year)
    e_eu27 = df.iloc[50:77, year_i]
    e_eu27.index = e_eu27.index.str.lstrip()

    e_ratio = e_country / e_eu27

    return pd.Series({k: e_ratio[v] for k, v in sub2sect.items()})


def industry_production_per_country(country, year, eurostat_dir, jrc_dir, snakemake):
    def get_sector_data(sector, country):
        jrc_country = jrc_names.get(country, country)
        fn = f"{jrc_dir}/{jrc_country}/JRC-IDEES-2021_Industry_{jrc_country}.xlsx"
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

    ct = "EU27" if country not in eu27 else country
    demand = pd.concat([get_sector_data(s, ct) for s in sect2sub])

    if country not in eu27:
        demand *= get_energy_ratio(
            country,
            eurostat_dir,
            jrc_dir,
            year,
            snakemake,
        )

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
        snakemake=snakemake,
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
    # ammonia data from 2018-2022
    ammonia = pd.read_csv(snakemake.input.ammonia_production, index_col=0)

    there = ammonia.index.intersection(demand.index)
    missing = demand.index.symmetric_difference(there)

    logger.info(f"Following countries have no ammonia demand: {missing.tolist()}")

    demand["Ammonia"] = 0.0

    year_to_use = min(max(year, 2018), 2022)
    if year_to_use != year:
        logger.info(
            f"Year {year} outside data range. Using data from {year_to_use} for ammonia production."
        )
    demand.loc[there, "Ammonia"] = ammonia.loc[there, str(year_to_use)]

    demand["Basic chemicals"] -= demand["Ammonia"]

    # EE, HR and LT got negative demand through subtraction - poor data
    col = "Basic chemicals"
    demand[col] = demand[col].clip(lower=0.0)

    # assume HVC, methanol, chlorine production proportional to non-ammonia basic chemicals
    distribution_key = (
        demand["Basic chemicals"]
        / params["basic_chemicals_without_NH3_production_today"]
        / 1e3
    )
    demand["HVC"] = params["HVC_production_today"] * 1e3 * distribution_key
    demand["Chlorine"] = params["chlorine_production_today"] * 1e3 * distribution_key
    demand["Methanol"] = params["methanol_production_today"] * 1e3 * distribution_key

    demand.drop(columns=["Basic chemicals"], inplace=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industrial_production_per_country")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    countries = snakemake.params.countries

    year = snakemake.params.industry["reference_year"]

    params = snakemake.params.industry

    jrc_dir = snakemake.input.jrc
    eurostat_dir = snakemake.input.eurostat

    demand = industry_production(countries, year, eurostat_dir, jrc_dir)

    separate_basic_chemicals(demand, year)

    demand.fillna(0.0, inplace=True)

    fn = snakemake.output.industrial_production_per_country
    demand.to_csv(fn, float_format="%.2f")
