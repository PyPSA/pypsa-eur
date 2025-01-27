# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build industrial energy demand per country.

Description
-------

This rule uses the industrial_production_per_country.csv file and the JRC-IDEES data to derive an energy demand per country and sector. If the country is not in the EU28, an average energy demand depending on the production volume is derived.
For each country and each subcategory of

- Alumina production
- Aluminium - primary production
- Aluminium - secondary production
- Ammonia
- Cement
- Ceramics & other NMM
- Chlorine
- Electric arc
- Food, beverages and tobacco
- Glass production
- HVC
- Integrated steelworks
- Machinery equipment
- Methanol
- Other industrial sectors
- Other chemicals
- Other non-ferrous metals
- Paper production
- Pharmaceutical products etc.
- Printing and media reproduction
- Pulp production
- Textiles and leather
- Transport equipment
- Wood and wood products

the output file contains the energy demand in TWh/a for the following carriers

- biomass
- electricity
- gas
- heat
- hydrogen
- liquid
- other
- solid
- waste
"""

import logging
import multiprocessing as mp
from functools import partial

import country_converter as coco
import pandas as pd
from _helpers import configure_logging, set_scenario_config
from tqdm import tqdm

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()

ktoe_to_twh = 0.011630

# name in JRC-IDEES Energy Balances
sector_sheets = {
    "Integrated steelworks": "FC_IND_IS_BF_E",
    "Electric arc": "FC_IND_IS_EA_E",
    "Alumina production": "FC_IND_NFM_AM_E",
    "Aluminium - primary production": "FC_IND_NFM_PA_E",
    "Aluminium - secondary production": "FC_IND_NFM_SA_E",
    "Other non-ferrous metals": "FC_IND_NFM_OM_E",
    "Basic chemicals": "FC_IND_CPC_BC_E",
    "Other chemicals": "FC_IND_CPC_OC_E",
    "Pharmaceutical products etc.": "FC_IND_CPC_PH_E",
    "Basic chemicals feedstock": "FC_IND_CPC_NE",
    "Cement": "FC_IND_NMM_CM_E",
    "Ceramics & other NMM": "FC_IND_NMM_CR_E",
    "Glass production": "FC_IND_NMM_GL_E",
    "Pulp production": "FC_IND_PPP_PU_E",
    "Paper production": "FC_IND_PPP_PA_E",
    "Printing and media reproduction": "FC_IND_PPP_PR_E",
    "Food, beverages and tobacco": "FC_IND_FBT_E",
    "Transport equipment": "FC_IND_TE_E",
    "Machinery equipment": "FC_IND_MAC_E",
    "Textiles and leather": "FC_IND_TL_E",
    "Wood and wood products": "FC_IND_WP_E",
    "Mining and quarrying": "FC_IND_MQ_E",
    "Construction": "FC_IND_CON_E",
    "Non-specified": "FC_IND_NSP_E",
}


fuels = {
    "Total": "all",
    "Solid fossil fuels": "solid",
    "Peat and peat products": "solid",
    "Oil shale and oil sands": "solid",
    "Oil and petroleum products": "liquid",
    "Manufactured gases": "gas",
    "Natural gas": "gas",
    "Nuclear heat": "heat",
    "Heat": "heat",
    "Renewables and biofuels": "biomass",
    "Non-renewable waste": "waste",
    "Electricity": "electricity",
}

eu27 = cc.EU27as("ISO2").ISO2.tolist()

jrc_names = {"GR": "EL", "GB": "UK"}


def industrial_energy_demand_per_country(country, year, jrc_dir, endogenous_ammonia):
    jrc_country = jrc_names.get(country, country)
    fn = f"{jrc_dir}/{jrc_country}/JRC-IDEES-2021_EnergyBalance_{jrc_country}.xlsx"

    sheets = list(sector_sheets.values())
    df_dict = pd.read_excel(fn, sheet_name=sheets, index_col=0)

    def get_subsector_data(sheet):
        df = df_dict[sheet][year].groupby(fuels).sum()

        df["hydrogen"] = 0.0

        # ammonia is handled separately
        if endogenous_ammonia:
            df["ammonia"] = 0.0

        df["other"] = df["all"] - df.loc[df.index != "all"].sum()

        return df

    df = pd.concat(
        {sub: get_subsector_data(sheet) for sub, sheet in sector_sheets.items()}, axis=1
    )

    sel = ["Mining and quarrying", "Construction", "Non-specified"]
    df["Other industrial sectors"] = df[sel].sum(axis=1)
    df["Basic chemicals"] += df["Basic chemicals feedstock"]

    df.drop(columns=sel + ["Basic chemicals feedstock"], index="all", inplace=True)

    df *= ktoe_to_twh

    return df


def separate_basic_chemicals(demand, production):
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

    demand["Chlorine"] = chlorine.unstack().reindex(index=demand.index, fill_value=0.0)
    demand["Methanol"] = methanol.unstack().reindex(index=demand.index, fill_value=0.0)

    demand["HVC"] = demand["Basic chemicals"] - demand["Methanol"] - demand["Chlorine"]

    # Deal with ammonia separately, depending on whether it is modelled endogenously.
    ammonia_exo = pd.DataFrame(
        {
            "hydrogen": production["Ammonia"] * params["MWh_H2_per_tNH3_electrolysis"],
            "electricity": production["Ammonia"]
            * params["MWh_elec_per_tNH3_electrolysis"],
        }
    ).T

    if snakemake.params.ammonia:
        ammonia = pd.DataFrame(
            {"ammonia": production["Ammonia"] * params["MWh_NH3_per_tNH3"]}
        ).T
    else:
        ammonia = ammonia_exo

    demand["Ammonia"] = ammonia.unstack().reindex(index=demand.index, fill_value=0.0)
    demand["HVC"] -= ammonia_exo.unstack().reindex(index=demand.index, fill_value=0.0)

    demand.drop(columns="Basic chemicals", inplace=True)

    demand["HVC"] = demand["HVC"].clip(lower=0)

    return demand


def add_non_eu27_industrial_energy_demand(countries, demand, production):
    non_eu27 = countries.difference(eu27)
    if non_eu27.empty:
        return demand

    eu27_production = production.loc[countries.intersection(eu27)].sum()
    eu27_energy = demand.groupby(level=1).sum()
    eu27_averages = eu27_energy / eu27_production

    demand_non_eu27 = pd.concat(
        {k: v * eu27_averages for k, v in production.loc[non_eu27].iterrows()}
    )

    return pd.concat([demand, demand_non_eu27])


def industrial_energy_demand(countries, year):
    nprocesses = snakemake.threads
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    func = partial(
        industrial_energy_demand_per_country,
        year=year,
        jrc_dir=snakemake.input.jrc,
        endogenous_ammonia=snakemake.params.ammonia,
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


def add_coke_ovens(demand, fn, year, factor=0.75):
    """
    Adds the energy consumption of coke ovens to the energy demand for
    integrated steelworks.

    This function reads the energy consumption data for coke ovens from a
    CSV file, processes it to match the structure of the `demand` DataFrame,
    and then adds a specified share of  this energy consumption to the energy
    demand for integrated steelworks.

    The `factor`  parameter controls what proportion of the coke ovens' energy
    consumption should be attributed to the iron and steel production.
    The default value of 75% is based on https://doi.org/10.1016/j.erss.2022.102565

    Parameters
    ----------
    demand (pd.DataFrame): A pandas DataFrame containing energy demand data
                           with a multi-level column index where one of the
                           levels corresponds to "Integrated steelworks".
    fn (str): The file path to the CSV file containing the coke ovens energy
              consumption data.
    year (int): The year for which the coke ovens data should be selected.
    factor (float, optional): The proportion of coke ovens energy consumption to add to the
                              integrated steelworks demand. Defaults to 0.75.

    Returns
    -------
    pd.DataFrame: The updated `demand` DataFrame with the coke ovens energy
    consumption added to the integrated steelworks energy demand.
    """

    df = pd.read_csv(fn, index_col=[0, 1]).xs(year, level=1)
    df = df.rename(columns={"Total all products": "Total"})[fuels.keys()]
    df = df.rename(columns=fuels).T.groupby(level=0).sum().T
    df["other"] = df["all"] - df.loc[:, df.columns != "all"].sum(axis=1)
    df = df.T.reindex_like(demand.xs("Integrated steelworks", axis=1, level=1)).fillna(
        0
    )
    sel = demand.columns.get_level_values(1) == "Integrated steelworks"
    demand.loc[:, sel] += factor * df.values

    return demand


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industrial_energy_demand_per_country_today")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params.industry
    year = params.get("reference_year", 2019)
    countries = pd.Index(snakemake.params.countries)

    demand = industrial_energy_demand(countries.intersection(eu27), year)

    # output in MtMaterial/a
    production = (
        pd.read_csv(snakemake.input.industrial_production_per_country, index_col=0)
        / 1e3
    )

    demand = separate_basic_chemicals(demand, production)

    demand = add_non_eu27_industrial_energy_demand(countries, demand, production)

    # for format compatibility
    demand = demand.stack(future_stack=True).unstack(level=[0, 2])

    # add energy consumption of coke ovens
    demand = add_coke_ovens(demand, snakemake.input.transformation_output_coke, year)

    # style and annotation
    demand.index.name = "TWh/a"
    demand.sort_index(axis=1, inplace=True)

    fn = snakemake.output.industrial_energy_demand_per_country_today
    demand.to_csv(fn)
