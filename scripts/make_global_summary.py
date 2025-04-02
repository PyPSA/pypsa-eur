# SPDX-FileCopyrightText: : 2025 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Create summary CSV files for all scenario runs including costs, capacities,
capacity factors, curtailment, energy balances, prices and other metrics.
"""

import logging

import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

INDEX_COLS = {
    "nodal_costs": 4,
    "nodal_capacities": 3,
    "nodal_capacity_factors": 3,
    "capacity_factors": 2,
    "costs": 3,
    "capacities": 2,
    "curtailment": 1,
    "energy": 2,
    "energy_balance": 3,
    "nodal_energy_balance": 4,
    "prices": 1,
    "weighted_prices": 1,
    "market_values": 1,
    "metrics": 1,
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("make_global_summary")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    for kind in snakemake.output.keys():
        logger.info(f"Creating global summary for {kind}")

        summaries_dict = {
            (cluster, opt + sector_opt, planning_horizon): "results/"
            + snakemake.params.RDIR
            + f"csvs/individual/{kind}_s_{cluster}_{opt}_{sector_opt}_{planning_horizon}.csv"
            for cluster in snakemake.params.scenario["clusters"]
            for opt in snakemake.params.scenario["opts"]
            for sector_opt in snakemake.params.scenario["sector_opts"]
            for planning_horizon in snakemake.params.scenario["planning_horizons"]
        }

        columns = pd.MultiIndex.from_tuples(
            summaries_dict.keys(),
            names=["cluster", "opt", "planning_horizon"],
        )

        summaries = []
        for _, filename in summaries_dict.items():
            s = pd.read_csv(filename, index_col=list(range(INDEX_COLS[kind])))
            summaries.append(s)

        summaries = pd.concat(summaries, axis=1)

        summaries.columns = columns

        summaries.sort_index().to_csv(snakemake.output[kind])

        del summaries
