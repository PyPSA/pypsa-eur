# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
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
    "costs": [0, 1, 2],
    "supply_energy": [0, 1, 2],
    "metrics": [0],
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("make_global_summary", kind="costs")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    for kind in ["costs", "supply_energy", "metrics"]:

        logger.info(f"Creating global summary for {kind}")

        summaries_dict = {
            (cluster, ll, opt + sector_opt, planning_horizon): "results/"
            + snakemake.params.RDIR
            + f"csvs/individual/{kind}_s_{cluster}_l{ll}_{opt}_{sector_opt}_{planning_horizon}.csv"
            for cluster in snakemake.params.scenario["clusters"]
            for opt in snakemake.params.scenario["opts"]
            for sector_opt in snakemake.params.scenario["sector_opts"]
            for ll in snakemake.params.scenario["ll"]
            for planning_horizon in snakemake.params.scenario["planning_horizons"]
        }

        columns = pd.MultiIndex.from_tuples(
            summaries_dict.keys(),
            names=["cluster", "ll", "opt", "planning_horizon"],
        )

        summaries = []
        for _, filename in summaries_dict.items():
            s = pd.read_csv(filename, index_col=INDEX_COLS[kind])
            s.index.names = [None] * s.index.nlevels
            summaries.append(s)

        summaries = pd.concat(summaries, axis=1)

        summaries.columns = columns

        summaries.sort_index().to_csv(snakemake.output[kind])

        del summaries
