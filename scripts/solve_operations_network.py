# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Solves linear optimal dispatch in hourly resolution using the capacities of a
previous capacity expansion in rule [solve_network][].

Loads `networks/solved_{horizon}.nc`, fixes the optimal capacities, and
re-solves dispatch only (optionally with rolling horizon). Note that
`load_shedding` is incompatible with operational solving, as its components
already exist in the solved network.
"""

import logging
from functools import partial
from pathlib import Path

import numpy as np
import pypsa

import scripts.solve_network as solve_network
from scripts._benchmark import memory_logger
from scripts._helpers import configure_logging, set_scenario_config
from scripts.solve_network import (
    collect_kwargs,
    create_optimization_model,
    extra_functionality,
    prepare_network,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_operations_network",
            configfiles="config/test/config.electricity.yaml",
            horizon=2030,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    solve_network.snakemake = snakemake

    cf_solving = snakemake.params.solving["options"]
    cf_operations = snakemake.params.solving["operations"]
    np.random.seed(cf_solving.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)
    planning_horizons = snakemake.wildcards["horizon"]
    resource_dir = Path(snakemake.input.network).resolve().parents[1]

    rolling_horizon = cf_operations["rolling_horizon"]

    n.optimize.fix_optimal_capacities()

    prepare_network(
        n,
        solve_opts=cf_solving,
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params["sector"]["limit_max_growth"],
        resource_dir=resource_dir,
        rolling_horizon=rolling_horizon,
    )

    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=logging_frequency
    ) as mem:
        if rolling_horizon:
            logger.info("Solving operations network with rolling horizon...")
            all_kwargs, _ = collect_kwargs(
                snakemake.config,
                snakemake.params.solving,
                planning_horizons,
                log_fn=snakemake.log.solver,
                mode="rolling_horizon",
                horizon=cf_operations["horizon"],
                overlap=cf_operations["overlap"],
            )
            n.config = snakemake.config
            n.params = snakemake.params
            all_kwargs["extra_functionality"] = partial(
                extra_functionality, planning_horizons=planning_horizons
            )
            n.optimize.optimize_with_rolling_horizon(**all_kwargs)
        else:
            logger.info("Solving operations network...")
            model_kwargs, solve_kwargs = collect_kwargs(
                snakemake.config,
                snakemake.params.solving,
                planning_horizons,
                log_fn=snakemake.log.solver,
                mode="single",
            )
            create_optimization_model(
                n,
                config=snakemake.config,
                params=snakemake.params,
                model_kwargs=model_kwargs,
                solve_kwargs=solve_kwargs,
                planning_horizons=planning_horizons,
            )
            n.optimize.solve_model(**solve_kwargs)

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)
