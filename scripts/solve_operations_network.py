# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Solves linear optimal dispatch in hourly resolution using the capacities of
previous capacity expansion in rule :mod:`solve_network`.
"""

import logging

import numpy as np
import pypsa
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from solve_network import prepare_network, solve_network

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_operations_network",
            configfiles="test/config.electricity.yaml",
            opts="",
            clusters="5",
            ll="v1.5",
            sector_opts="",
            planning_horizons="",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solve_opts = snakemake.params.options

    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)

    n.optimize.fix_optimal_capacities()
    n = prepare_network(n, solve_opts, config=snakemake.config)
    n = solve_network(
        n,
        config=snakemake.config,
        params=snakemake.params,
        solving=snakemake.params.solving,
        log_fn=snakemake.log.solver,
    )

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
