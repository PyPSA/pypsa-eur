# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build specific energy consumption per industrial subsector and country with process heat split into temperature bands in MWh/t.

Today's ratios are derived per country by dividing the energy demand from [build_industrial_energy_demand_per_country_today][] by the production from [build_industrial_production_per_country][]. For each subsector listed in `data/industry_heat_temperature_shares.csv`, all fuels and purchased heat are summed and split into the bands below 100 °C (`heat`), 100-200 °C, 200-500 °C and above 500 °C, while electricity keeps today's ratio. Feedstocks are not heat demand: liquid fuels of HVC keep their exogenous ratio, and the natural gas for ammonia, which the statistics report with basic chemicals, is removed. The heat demand is given as fuel input, which overestimates useful heat by the losses of today's boilers and furnaces. The steel sector is not covered.

References
----------
- [Fleiter et al. (2025)](https://doi.org/10.1002/ente.202300981). Hydrogen Infrastructure in the Future CO2-Neutral European Energy System.
"""

import logging

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

FUELS = ["solid", "liquid", "gas", "heat", "biomass", "waste", "other"]
REPLACED = ["coal", "coke", "biomass", "methane", "hydrogen", "naphtha"]


def build_industry_sector_ratios_endogenous(
    demand: pd.DataFrame,
    production: pd.DataFrame,
    shares: pd.DataFrame,
    gas_per_ammonia: float,
) -> pd.DataFrame:
    """
    Split today's fuel use per tonne of product into temperature bands.

    Parameters
    ----------
    demand : pd.DataFrame
        Today's energy demand (TWh/a) with fuels as index and country and
        subsector as columns.
    production : pd.DataFrame
        Today's production (kt/a) with countries as index and subsectors as
        columns.
    shares : pd.DataFrame
        Share of heat demand per temperature band with subsectors as index.
    gas_per_ammonia : float
        Natural gas demand of ammonia production (MWh/t).

    Returns
    -------
    pd.DataFrame
        Ratios (MWh/t) for the subsectors in `shares`. Carriers replaced by
        heat are zero. Missing values keep the exogenous ratio.
    """
    demand = demand.copy()
    hvc_gas = demand.loc["gas"].xs("HVC", level=1)
    ammonia_gas = production["Ammonia"] * gas_per_ammonia / 1e3
    hvc_gas = hvc_gas.sub(ammonia_gas.reindex(hvc_gas.index), fill_value=0).clip(0)
    demand.loc["gas", (hvc_gas.index, "HVC")] = hvc_gas.values

    production = (production / 1e3).stack()
    production.index.names = [None, None]

    ratios = demand.div(production, axis=1).replace([np.inf, -np.inf], np.nan)
    ratios = ratios.loc[:, ratios.columns.get_level_values(1).isin(shares.index)]
    # subsectors without production today keep the exogenous ratios
    ratios = ratios.dropna(how="all", axis=1).fillna(0.0)

    fuels = ratios.loc[FUELS]
    feedstock = fuels.columns.get_level_values(1) == "HVC"
    fuels.loc["liquid", feedstock] = 0.0

    shares = shares.div(shares.sum(axis=1), axis=0).rename(columns={"heat<100": "heat"})
    sector_shares = shares.T.reindex(columns=ratios.columns, level=1)

    endogenous = pd.concat(
        [
            ratios.loc[["electricity"]].rename({"electricity": "elec"}),
            pd.DataFrame(0.0, index=REPLACED, columns=ratios.columns),
            sector_shares * fuels.sum(),
        ]
    )
    endogenous.loc["naphtha", feedstock] = np.nan
    endogenous.index.name = "MWh/t"

    return endogenous


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_industry_sector_ratios_endogenous")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    demand = pd.read_csv(
        snakemake.input.industrial_energy_demand_per_country_today,
        header=[0, 1],
        index_col=0,
    )
    production = pd.read_csv(
        snakemake.input.industrial_production_per_country, index_col=0
    )
    shares = pd.read_csv(snakemake.input.temperature_shares, index_col=0)
    shares = shares.drop(columns="source").drop(snakemake.params.exclude_sectors)

    endogenous = build_industry_sector_ratios_endogenous(
        demand, production, shares, snakemake.params.gas_per_ammonia
    )
    endogenous.to_csv(snakemake.output.industry_sector_ratios_endogenous)
