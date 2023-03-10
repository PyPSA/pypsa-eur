# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Solves linear optimal dispatch in hourly resolution using the capacities of
previous capacity expansion in rule :mod:`solve_network`.

Relevant Settings
-----------------

.. code:: yaml

    solving:
        tmpdir:
        options:
            formulation:
            clip_p_max_pu:
            load_shedding:
            noisy_costs:
            nhours:
            min_iterations:
            max_iterations:
        solver:
            name:
            (solveroptions):

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`solving_cf`

Inputs
------

- ``networks/elec_s{simpl}_{clusters}.nc``: confer :ref:`cluster`
- ``results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: confer :ref:`solve`

Outputs
-------

- ``results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op.nc``: Solved PyPSA network for optimal dispatch including optimisation results

Description
-----------
"""

import logging

import numpy as np
import pypsa
from _helpers import (
    configure_logging,
    override_component_attrs,
    update_config_with_sector_opts,
)
from solve_network import prepare_network, solve_network
from vresutils.benchmark import memory_logger

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_operations_network",
            configfiles="test/config.test1.yaml",
            simpl="",
            opts="",
            clusters="5",
            ll="v1.5",
            sector_opts="",
            planning_horizons="",
        )

    configure_logging(snakemake)
    update_config_with_sector_opts(snakemake.config, snakemake.wildcards.sector_opts)

    opts = (snakemake.wildcards.opts + "-" + snakemake.wildcards.sector_opts).split("-")
    opts = [o for o in opts if o != ""]
    solve_opts = snakemake.config["solving"]["options"]

    np.random.seed(solve_opts.get("seed", 123))

    fn = getattr(snakemake.log, "memory", None)
    with memory_logger(filename=fn, interval=30.0) as mem:
        if "overrides" in snakemake.input:
            overrides = override_component_attrs(snakemake.input.overrides)
            n = pypsa.Network(
                snakemake.input.network, override_component_attrs=overrides
            )
        else:
            n = pypsa.Network(snakemake.input.network)

        n.optimize.fix_optimal_capacities()
        n = prepare_network(n, solve_opts, config=snakemake.config)
        n = solve_network(
            n, config=snakemake.config, opts=opts, log_fn=snakemake.log.solver
        )

        n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
