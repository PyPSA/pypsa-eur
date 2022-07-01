"""Solve network."""

import pypsa

import numpy as np
import pandas as pd

from pypsa.linopt import get_var, linexpr, define_constraints

from pypsa.linopf import network_lopf, ilopf

from vresutils.benchmark import memory_logger

from helper import override_component_attrs

import logging
logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def add_land_use_constraint(n):

    if 'm' in snakemake.wildcards.clusters:
        _add_land_use_constraint_m(n)
    else:
        _add_land_use_constraint(n)


def _add_land_use_constraint(n):
    #warning: this will miss existing offwind which is not classed AC-DC and has carrier 'offwind'

    for carrier in ['solar', 'onwind', 'offwind-ac', 'offwind-dc']:
        existing = n.generators.loc[n.generators.carrier==carrier,"p_nom"].groupby(n.generators.bus.map(n.buses.location)).sum()
        existing.index += " " + carrier + "-" + snakemake.wildcards.planning_horizons
        n.generators.loc[existing.index,"p_nom_max"] -= existing

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def _add_land_use_constraint_m(n):
    # if generators clustering is lower than network clustering, land_use accounting is at generators clusters

    planning_horizons = snakemake.config["scenario"]["planning_horizons"]
    grouping_years = snakemake.config["existing_capacities"]["grouping_years"]
    current_horizon = snakemake.wildcards.planning_horizons

    for carrier in ['solar', 'onwind', 'offwind-ac', 'offwind-dc']:

        existing = n.generators.loc[n.generators.carrier==carrier,"p_nom"]
        ind = list(set([i.split(sep=" ")[0] + ' ' + i.split(sep=" ")[1] for i in existing.index]))

        previous_years = [
            str(y) for y in
            planning_horizons + grouping_years
            if y < int(snakemake.wildcards.planning_horizons)
        ]

        for p_year in previous_years:
            ind2 = [i for i in ind if i + " " + carrier + "-" + p_year in existing.index]
            sel_current = [i + " " + carrier + "-" + current_horizon for i in ind2]
            sel_p_year = [i + " " + carrier + "-" + p_year for i in ind2]
            n.generators.loc[sel_current, "p_nom_max"] -= existing.loc[sel_p_year].rename(lambda x: x[:-4] + current_horizon)

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def prepare_network(n, solve_opts=None):

    if 'clip_p_max_pu' in solve_opts:
        for df in (n.generators_t.p_max_pu, n.generators_t.p_min_pu, n.storage_units_t.inflow):
            df.where(df>solve_opts['clip_p_max_pu'], other=0., inplace=True)

    if solve_opts.get('load_shedding'):
        n.add("Carrier", "Load")
        n.madd("Generator", n.buses.index, " load",
               bus=n.buses.index,
               carrier='load',
               sign=1e-3, # Adjust sign to measure p and p_nom in kW instead of MW
               marginal_cost=1e2, # Eur/kWh
               # intersect between macroeconomic and surveybased
               # willingness to pay
               # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
               p_nom=1e9 # kW
        )

    if solve_opts.get('noisy_costs'):
        for t in n.iterate_components():
            #if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if 'marginal_cost' in t.df:
                np.random.seed(174)
                t.df['marginal_cost'] += 1e-2 + 2e-3 * (np.random.random(len(t.df)) - 0.5)

        for t in n.iterate_components(['Line', 'Link']):
            np.random.seed(123)
            t.df['capital_cost'] += (1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)) * t.df['length']

    if solve_opts.get('nhours'):
        nhours = solve_opts['nhours']
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760./nhours

    if snakemake.config['foresight'] == 'myopic':
        add_land_use_constraint(n)

    return n


def add_battery_constraints(n):

    chargers_b = n.links.carrier.str.contains("battery charger")
    chargers = n.links.index[chargers_b & n.links.p_nom_extendable]
    dischargers = chargers.str.replace("charger", "discharger")

    if chargers.empty or ('Link', 'p_nom') not in n.variables.index:
        return

    link_p_nom = get_var(n, "Link", "p_nom")

    lhs = linexpr((1,link_p_nom[chargers]),
                  (-n.links.loc[dischargers, "efficiency"].values,
                   link_p_nom[dischargers].values))

    define_constraints(n, lhs, "=", 0, 'Link', 'charger_ratio')


def add_chp_constraints(n):

    electric_bool = (n.links.index.str.contains("urban central")
                     & n.links.index.str.contains("CHP")
                     & n.links.index.str.contains("electric"))
    heat_bool = (n.links.index.str.contains("urban central")
                 & n.links.index.str.contains("CHP")
                 & n.links.index.str.contains("heat"))

    electric = n.links.index[electric_bool]
    heat = n.links.index[heat_bool]

    electric_ext = n.links.index[electric_bool & n.links.p_nom_extendable]
    heat_ext = n.links.index[heat_bool & n.links.p_nom_extendable]

    electric_fix = n.links.index[electric_bool & ~n.links.p_nom_extendable]
    heat_fix = n.links.index[heat_bool & ~n.links.p_nom_extendable]

    link_p = get_var(n, "Link", "p")

    if not electric_ext.empty:

        link_p_nom = get_var(n, "Link", "p_nom")

        #ratio of output heat to electricity set by p_nom_ratio
        lhs = linexpr((n.links.loc[electric_ext, "efficiency"]
                       *n.links.loc[electric_ext, "p_nom_ratio"],
                       link_p_nom[electric_ext]),
                      (-n.links.loc[heat_ext, "efficiency"].values,
                       link_p_nom[heat_ext].values))

        define_constraints(n, lhs, "=", 0, 'chplink', 'fix_p_nom_ratio')

        #top_iso_fuel_line for extendable
        lhs = linexpr((1,link_p[heat_ext]),
                      (1,link_p[electric_ext].values),
                      (-1,link_p_nom[electric_ext].values))

        define_constraints(n, lhs, "<=", 0, 'chplink', 'top_iso_fuel_line_ext')

    if not electric_fix.empty:

        #top_iso_fuel_line for fixed
        lhs = linexpr((1,link_p[heat_fix]),
                      (1,link_p[electric_fix].values))

        rhs = n.links.loc[electric_fix, "p_nom"].values

        define_constraints(n, lhs, "<=", rhs, 'chplink', 'top_iso_fuel_line_fix')

    if not electric.empty:

        #backpressure
        lhs = linexpr((n.links.loc[electric, "c_b"].values
                       *n.links.loc[heat, "efficiency"],
                       link_p[heat]),
                      (-n.links.loc[electric, "efficiency"].values,
                       link_p[electric].values))

        define_constraints(n, lhs, "<=", 0, 'chplink', 'backpressure')

def basename(x):
     return x.split("-2")[0]

def add_pipe_retrofit_constraint(n):
    """Add constraint for retrofitting existing CH4 pipelines to H2 pipelines."""

    gas_pipes_i = n.links.query("carrier == 'gas pipeline' and p_nom_extendable").index
    h2_retrofitted_i = n.links.query("carrier == 'H2 pipeline retrofitted' and p_nom_extendable").index

    if h2_retrofitted_i.empty or gas_pipes_i.empty: return

    link_p_nom = get_var(n, "Link", "p_nom")

    CH4_per_H2 = 1 / n.config["sector"]["H2_retrofit_capacity_per_CH4"]
    fr = "H2 pipeline retrofitted"
    to = "gas pipeline"

    pipe_capacity = n.links.loc[gas_pipes_i, 'p_nom'].rename(basename)

    lhs = linexpr(
        (CH4_per_H2, link_p_nom.loc[h2_retrofitted_i].rename(index=lambda x: x.replace(fr, to))),
        (1, link_p_nom.loc[gas_pipes_i])
    )

    lhs.rename(basename, inplace=True)
    define_constraints(n, lhs, "=", pipe_capacity, 'Link', 'pipe_retrofit')


def add_co2_sequestration_limit(n, sns):

    co2_stores = n.stores.loc[n.stores.carrier=='co2 stored'].index

    if co2_stores.empty or ('Store', 'e') not in n.variables.index:
        return

    vars_final_co2_stored = get_var(n, 'Store', 'e').loc[sns[-1], co2_stores]

    lhs = linexpr((1, vars_final_co2_stored)).sum()

    limit = n.config["sector"].get("co2_sequestration_potential", 200) * 1e6
    for o in opts:
        if not "seq" in o: continue
        limit = float(o[o.find("seq")+3:])
        break

    name = 'co2_sequestration_limit'
    sense = "<="

    n.add("GlobalConstraint", name, sense=sense, constant=limit,
          type=np.nan, carrier_attribute=np.nan)

    define_constraints(n, lhs, sense, limit, 'GlobalConstraint',
                       'mu', axes=pd.Index([name]), spec=name)


def extra_functionality(n, snapshots):
    add_battery_constraints(n)
    add_pipe_retrofit_constraint(n)
    add_co2_sequestration_limit(n, snapshots)


def solve_network(n, config, opts='', **kwargs):
    solver_options = config['solving']['solver'].copy()
    solver_name = solver_options.pop('name')
    cf_solving = config['solving']['options']
    track_iterations = cf_solving.get('track_iterations', False)
    min_iterations = cf_solving.get('min_iterations', 4)
    max_iterations = cf_solving.get('max_iterations', 6)
    keep_shadowprices = cf_solving.get('keep_shadowprices', True)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    if cf_solving.get('skip_iterations', False):
        network_lopf(n, solver_name=solver_name, solver_options=solver_options,
                     extra_functionality=extra_functionality,
                     keep_shadowprices=keep_shadowprices, **kwargs)
    else:
        ilopf(n, solver_name=solver_name, solver_options=solver_options,
              track_iterations=track_iterations,
              min_iterations=min_iterations,
              max_iterations=max_iterations,
              extra_functionality=extra_functionality,
              keep_shadowprices=keep_shadowprices,
              **kwargs)
    return n


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'solve_network',
            simpl='',
            opts="",
            clusters="37",
            lv=1.0,
            sector_opts='168H-T-H-B-I-A-solar+p3-dist1',
            planning_horizons="2030",
        )

    logging.basicConfig(filename=snakemake.log.python,
                        level=snakemake.config['logging_level'])

    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        from pathlib import Path
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
    opts = snakemake.wildcards.opts.split('-')
    solve_opts = snakemake.config['solving']['options']

    fn = getattr(snakemake.log, 'memory', None)
    with memory_logger(filename=fn, interval=30.) as mem:

        overrides = override_component_attrs(snakemake.input.overrides)
        n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

        n = prepare_network(n, solve_opts)

        n = solve_network(n, config=snakemake.config, opts=opts,
                          solver_dir=tmpdir,
                          solver_logfile=snakemake.log.solver)

        if "lv_limit" in n.global_constraints.index:
            n.line_volume_limit = n.global_constraints.at["lv_limit", "constant"]
            n.line_volume_limit_dual = n.global_constraints.at["lv_limit", "mu"]

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
