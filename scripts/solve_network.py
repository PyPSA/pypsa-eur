"""
Solves linear optimal power flow for a network iteratively while updating reactances.

Relevant Settings
-----------------

.. code:: yaml

    (electricity:)
        (BAU_mincapacities:)
        (SAFE_reservemargin:)

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

    (plotting:)
        (conv_techs:)

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`electricity_cf`, :ref:`solving_cf`, :ref:`plotting_cf`

Inputs
------

- ``networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: confer :ref:`prepare`

Outputs
-------

- ``results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: Solved PyPSA network including optimisation results

    .. image:: ../img/results.png
        :scale: 40 %

Description
-----------

Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.
Additionaly some extra constraints from :mod:`prepare_network` are added.

Solving the network in multiple iterations is motivated through the dependence of transmission line capacities and impedances.
As lines are expanded their electrical parameters change, which renders the optimisation bilinear even if the power flow
equations are linearized.
To retain the computational advantage of continuous linear programming, a sequential linear programming technique
is used, where in between iterations the line impedances are updated.
Details (and errors made through this heuristic) are discussed in the paper

- Fabian Neumann and Tom Brown. `Heuristics for Transmission Expansion Planning in Low-Carbon Energy System Models <https://arxiv.org/abs/1907.10548>`_), *16th International Conference on the European Energy Market*, 2019. `arXiv:1907.10548 <https://arxiv.org/abs/1907.10548>`_.

.. warning::
    Capital costs of existing network components are not included in the objective function,
    since for the optimisation problem they are just a constant term (no influence on optimal result).

    Therefore, these capital costs are not included in ``network.objective``!

    If you want to calculate the full total annual system costs add these to the objective value.

.. tip::
    The rule :mod:`solve_all_networks` runs
    for all ``scenario`` s in the configuration file
    the rule :mod:`solve_network`.

"""

import logging
logger = logging.getLogger(__name__)
from _helpers import configure_logging

import numpy as np
import pandas as pd
import gc

import pypsa
from pypsa.linopf import (get_var, define_constraints, linexpr, join_exprs,
                          lopf, ilopf)

from vresutils.benchmark import memory_logger


def prepare_network(n, solve_opts=None):
    if solve_opts is None:
        solve_opts = snakemake.config['solving']['options']

    if 'clip_p_max_pu' in solve_opts:
        for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df>solve_opts['clip_p_max_pu'], other=0., inplace=True)

    if solve_opts.get('noisy_costs'):
        for t in n.iterate_components(n.one_port_components):
            #if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if 'marginal_cost' in t.df:
                t.df['marginal_cost'] += 1e-2 + 2e-3*(np.random.random(len(t.df)) - 0.5)

        for t in n.iterate_components(['Line', 'Link']):
            t.df['capital_cost'] += (1e-1 +
                2e-2*(np.random.random(len(t.df)) - 0.5)) * t.df['length']

    if solve_opts.get('nhours'):
        nhours = solve_opts['nhours']
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760./nhours

    return n


def add_opts_constraints(n, opts=None):
    if opts is None:
        opts = snakemake.wildcards.opts.split('-')

    if 'BAU' in opts:
        mincaps = snakemake.config['electricity']['BAU_mincapacities']
        lhs = (linexpr((1, get_var(n, 'Generator', 'p_nom')))
               .groupby(n.generators.carrier).apply(join_exprs))
        define_constraints(n, lhs, '<=', mincaps, 'Carrier', 'bau_mincaps')


    if 'SAFE' in opts:
        peakdemand = (1. + snakemake.config['electricity']['SAFE_reservemargin']) *\
                      n.loads_t.p_set.sum(axis=1).max()
        conv_techs = snakemake.config['plotting']['conv_techs']
        exist_conv_caps = n.generators.query('~p_nom_extendable & carrier in @conv_techs')\
                           .p_nom.sum()
        ext_gens_i = n.generators.query('carrier in @conv_techs & p_nom_extendable').index
        lhs = linexpr((1, get_var('n', 'Generator', 'p_nom')[ext_gens_i])).sum()
        rhs = peakdemand - exist_conv_caps
        define_constraints(n, lhs, '>=', rhs, 'Safe', 'mintotalcap')

    # Add constraints on the per-carrier capacity in each country
    if 'CCL' in opts:
        agg_p_nom_limits = snakemake.config['electricity'].get('agg_p_nom_limits')

        try:
            agg_p_nom_minmax = pd.read_csv(agg_p_nom_limits, index_col=list(range(2)))
        except IOError:
            logger.exception("Need to specify the path to a .csv file containing "
                             "aggregate capacity limits per country in "
                             "config['electricity']['agg_p_nom_limit'].")
        logger.info("Adding per carrier generation capacity constraints for "
                    "individual countries")

        gen_country = n.generators.bus.map(n.buses.country)
        # min, cc means country and carrier
        p_nom_per_cc = (pd.DataFrame(
                       {'p_nom': linexpr((1, get_var(n, 'Generator', 'p_nom'))),
                        'country': gen_country, 'carrier': n.generators.carrier})
                       .groupby(['country', 'carrier']).p_nom
                       .apply(join_exprs))


        def agg_p_nom_min_rule(model, country, carrier):
            min = agg_p_nom_minmax.at[(country, carrier), 'min']
            return ((sum(model.generator_p_nom[gen]
                        for gen in n.generators.index[(gen_country == country) & (n.generators.carrier == carrier)])
                    >= min)
                    if np.isfinite(min) else pypsa.opt.Constraint.Skip)

        def agg_p_nom_max_rule(model, country, carrier):
            max = agg_p_nom_minmax.at[(country, carrier), 'max']
            return ((sum(model.generator_p_nom[gen]
                        for gen in n.generators.index[(gen_country == country) & (n.generators.carrier == carrier)])
                    <= max)
                    if np.isfinite(max) else pypsa.opt.Constraint.Skip)

        n.model.agg_p_nom_min = pypsa.opt.Constraint(list(agg_p_nom_minmax.index), rule=agg_p_nom_min_rule)
        n.model.agg_p_nom_max = pypsa.opt.Constraint(list(agg_p_nom_minmax.index), rule=agg_p_nom_max_rule)



#def add_eps_storage_constraint(n):
#    if not hasattr(n, 'epsilon'):
#        n.epsilon = 1e-5
#    fix_sus_i = n.storage_units.index[~ n.storage_units.p_nom_extendable]
#    n.model.objective.expr += sum(n.epsilon * n.model.state_of_charge[su, n.snapshots[0]] for su in fix_sus_i)


def solve_network(n, config=None, solver_log=None, skip_iterating=False,
                  **kwargs):
    if config is None:
        config = snakemake.config['solving']
#    solve_opts = config['options']

    solver_options = config['solver'].copy()
    if solver_log is None:
        solver_log = snakemake.log.solver
    solver_name = solver_options.pop('name')
    if skip_iterating:
        lopf(n, solver_name=solver_name, solver_options=solver_options, **kwargs)
    else:
        ilopf(n, solver_name=solver_name, solver_options=solver_options, **kwargs)
    return n

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('solve_network', network='elec', simpl='',
                                  clusters='5', ll='copt', opts='Co2L-24H')
    configure_logging(snakemake)

    with memory_logger(filename=getattr(snakemake.log, 'memory', None),
                       interval=30.) as mem:
        n = pypsa.Network(snakemake.input[0])
        n = prepare_network(n)
        n = solve_network(n)
        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
