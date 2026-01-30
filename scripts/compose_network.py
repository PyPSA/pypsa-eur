# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Compose network by combining all electricity and sector components.
"""

import logging

import pypsa

from scripts._helpers import (
    configure_logging,
    load_costs,
    sanitize_custom_columns,
    set_scenario_config,
)
from scripts.add_brownfield import main as apply_brownfield
from scripts.add_electricity import main as add_electricity_components
from scripts.add_electricity import (
    sanitize_carriers,
    sanitize_locations,
)
from scripts.add_existing_baseyear import main as add_existing_capacities
from scripts.prepare_network import (
    apply_co2_budget_constraints,
    apply_temporal_aggregation,
    maybe_adjust_costs_and_potentials,
)
from scripts.prepare_network import main as prepare_network_for_solving
from scripts.prepare_perfect_foresight import main as prepare_perfect_foresight
from scripts.prepare_sector_network import (
    main as add_sector_components,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "compose_network",
            run="test-run",
            horizon="2050",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Extract configuration and parameters
    config = snakemake.config
    params = snakemake.params
    inputs = snakemake.input

    current_horizon = int(snakemake.wildcards.horizon)
    horizons = config["planning_horizons"]
    foresight = config["foresight"]
    sector_mode = params.sector["enabled"]

    is_first_horizon = current_horizon == horizons[0]

    n = pypsa.Network(snakemake.input.network)
    n_previous = None if is_first_horizon else pypsa.Network(inputs.network_previous)
    nyears = n.snapshot_weightings.objective.sum() / 8760.0
    costs = load_costs(snakemake.input.tech_costs)

    add_electricity_components(n, inputs, params, costs)

    if sector_mode:
        add_sector_components(n, inputs, params, costs, nyears, current_horizon)

    apply_temporal_aggregation(n, inputs, params)

    if foresight != "overnight" and is_first_horizon:
        add_existing_capacities(n, inputs, params, costs)

    if foresight == "myopic" and not is_first_horizon:
        apply_brownfield(n, n_previous, inputs, params, current_horizon)

    prepare_network_for_solving(n, inputs, params, costs, nyears)

    if foresight == "perfect":
        n = prepare_perfect_foresight(n, n_previous, params, current_horizon)

    apply_co2_budget_constraints(
        n, inputs=inputs, params=params, nyears=nyears, current_horizon=current_horizon
    )

    adjustments = params.adjustments
    maybe_adjust_costs_and_potentials(n, adjustments["electricity"], current_horizon)
    maybe_adjust_costs_and_potentials(n, adjustments["sector"], current_horizon)

    sanitize_custom_columns(n)
    sanitize_carriers(n, config)
    sanitize_locations(n)
    n.meta = dict(config, **dict(wildcards=dict(snakemake.wildcards)))
    n.consistency_check()

    logger.info(f"Exporting composed network for horizon {current_horizon}")
    n.export_to_netcdf(snakemake.output[0])
