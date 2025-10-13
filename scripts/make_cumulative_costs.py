# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Calculate cumulative costs for myopic foresight scenarios.
"""

import logging

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

idx = pd.IndexSlice
logger = logging.getLogger(__name__)


def calculate_cumulative_cost(costs, planning_horizons):
    cumulative_cost = pd.DataFrame(
        index=costs.sum().index,
        columns=pd.Series(data=np.arange(0, 0.1, 0.01), name="social discount rate"),
    )

    # discount cost and express them in money value of planning_horizons[0]
    for r in cumulative_cost.columns:
        cumulative_cost[r] = [
            costs.sum()[index] / ((1 + r) ** (index[-1] - planning_horizons[0]))
            for index in cumulative_cost.index
        ]

    # integrate cost throughout the transition path
    for r in cumulative_cost.columns:
        for cluster in cumulative_cost.index.get_level_values(level=0).unique():
            for sector_opts in cumulative_cost.index.get_level_values(level=1).unique():
                cumulative_cost.loc[(cluster, sector_opts, "cumulative cost"), r] = (
                    np.trapz(
                        cumulative_cost.loc[
                            idx[cluster, sector_opts, planning_horizons], r
                        ].values,
                        x=planning_horizons,
                    )
                )

    return cumulative_cost


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "make_cumulative_costs",
            clusters="5",
            opts="",
            sector_opts="",
            configfiles="config/test/config.myopic.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    costs = pd.read_csv(snakemake.input.costs, index_col=[0, 1, 2], header=[0, 1, 2])

    # clean multiindex: handle empty opts/sector_opts & make planning horizons numeric
    costs.columns = pd.MultiIndex.from_tuples(
        [
            (lvl0, "" if lvl1.startswith("Unnamed:") else lvl1, pd.to_numeric(lvl2))
            for lvl0, lvl1, lvl2 in costs.columns
        ],
        names=["cluster", "opt", "planning_horizon"],
    )

    planning_horizons = snakemake.params.scenario["planning_horizons"]

    cumulative_cost = calculate_cumulative_cost(costs, planning_horizons)
    cumulative_cost.to_csv(snakemake.output[0])
