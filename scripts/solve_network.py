# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Solves optimal operation and capacity for a network with the option to
iteratively optimize while updating line reactances.

This script is used for optimizing the electrical network as well as the
sector coupled network.

Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning)
is provided in the
[documentation of PyPSA](https://docs.pypsa.org/latest/user-guide/network-optimization/).

The optimization is based on the `network.optimize` function.
Additionally, the custom math in `scripts/math/` is stored on the network as
PyPSA overrides.
"""

import importlib
import logging
import os
import re
import sys
from collections.abc import Callable
from functools import partial
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pypsa
import yaml
from linopy.constants import SolverStatus, TerminationCondition

from scripts._benchmark import memory_logger
from scripts._helpers import (
    configure_logging,
    get,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

pypsa.network.power_flow.logger.setLevel(logging.WARNING)

MATH = Path(__file__).parent / "math"


class ObjectiveValueError(Exception):
    pass


# TODO: move to compose_network
def add_land_use_constraint(n: pypsa.Network) -> None:
    """
    Add global constraints for tech capacity limit.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance

    Returns
    -------
    pypsa.Network
        Network with added land use constraints
    """
    logger.info("Add land-use constraint for perfect foresight")

    def check_p_min_p_max(p_nom_max: pd.Series) -> None:
        p_nom_min = n.generators[ext_i].groupby(grouper).sum().p_nom_min
        p_nom_min = p_nom_min.reindex(p_nom_max.index)
        check = (
            p_nom_min.groupby(level=[0, 1]).sum()
            > p_nom_max.groupby(level=[0, 1]).min()
        )
        if check.sum():
            logger.warning(
                f"summed p_min_pu values at node larger than technical potential {check[check].index}"
            )

    grouper = [n.generators.carrier, n.generators.bus, n.generators.build_year]
    ext_i = n.generators.p_nom_extendable
    # get technical limit per node and investment period
    p_nom_max = n.generators[ext_i].groupby(grouper).min().p_nom_max
    # drop carriers without tech limit
    p_nom_max = p_nom_max[~p_nom_max.isin([np.inf, np.nan])]
    carriers = p_nom_max.index.get_level_values(0).unique()
    gen_i = n.generators[(n.generators.carrier.isin(carriers)) & (ext_i)].index
    n.generators.loc[gen_i, "p_nom_min"] = 0
    # check minimum capacities
    check_p_min_p_max(p_nom_max)
    # adjust name to fit syntax of nominal constraint per bus
    df = p_nom_max.reset_index()
    df["name"] = df.apply(
        lambda row: (
            f"nom_max_{row['carrier']}"
            + (f"_{row['build_year']}" if row["build_year"] is not None else "")
        ),
        axis=1,
    )

    for name in df.name.unique():
        df_carrier = df[df.name == name]
        bus = df_carrier.bus
        n.buses.loc[bus, name] = df_carrier.p_nom_max.values


def solar_potential(
    n: pypsa.Network, config: dict, planning_horizons: str | None
) -> bool:
    """
    Share the solar potential of a bus among fixed-tilt and tracking solar.

    Extendable fixed-tilt solar holds the potential in its `p_nom_max`.
    Tracking solar takes more land per unit of capacity, by the ratio of the
    `capacity_per_sqkm` of both technologies.
    """
    solar = {"solar", "solar-hsat"}
    electricity = config["electricity"]
    if not (
        solar.issubset(electricity["renewable_carriers"])
        and solar.issubset(electricity["extendable_carriers"]["Generator"])
    ):
        return False
    renewable = config["renewable"]
    land_use = (
        renewable["solar"]["capacity_per_sqkm"]
        / renewable["solar-hsat"]["capacity_per_sqkm"]
    )
    gens = n.generators
    potential = gens.carrier.eq("solar") & gens.p_nom_extendable
    hsat = gens.carrier.eq("solar-hsat")
    n.generators["solar_potential"] = potential
    n.generators["solar_land_use"] = np.where(
        potential, 1.0, np.where(hsat, land_use, np.nan)
    )
    return True


def add_co2_sequestration_limit(
    n: pypsa.Network,
    limit_dict: dict[str, float],
    planning_horizons: str | None,
) -> None:
    """
    Add a global constraint on the amount of Mt CO2 that can be sequestered.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    limit_dict : dict[str, float]
        CO2 sequestration potential limit constraints by year.
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight
    """

    if not n.investment_periods.empty:
        nyears = n.snapshot_weightings.groupby(level="period").generators.sum() / 8760
        periods = n.investment_periods
        limit = pd.Series(
            {period: nyears[period] * get(limit_dict, period) for period in periods}
        )
        limit.index = limit.index.map(lambda s: f"co2_sequestration_limit-{s}")
        names = limit.index
    else:
        nyears = n.snapshot_weightings.generators.sum() / 8760
        limit = get(limit_dict, int(planning_horizons)) * nyears
        periods = np.nan
        names = "co2_sequestration_limit"

    n.add(
        "GlobalConstraint",
        names,
        sense=">=",
        constant=-limit * 1e6,
        type="operational_limit",
        carrier_attribute="co2 sequestered",
        investment_period=periods,
    )


def emission_stores(n: pypsa.Network, emissions: pd.Index) -> pd.DataFrame:
    """Return the non-cyclic stores that accumulate emissions on their bus."""
    bus_carrier = n.stores.bus.map(n.buses.carrier)
    stores = n.stores[bus_carrier.isin(emissions) & ~n.stores.e_cyclic]
    seeded = stores.index[stores.e_initial.ne(0) | stores.e_initial_per_period]
    if not seeded.empty:
        raise ValueError(
            f"Emission stores {list(seeded)} start from a non-zero level "
            "(e_initial or e_initial_per_period). The CO2 constraints read the "
            "store level as cumulative emissions, which only holds from zero."
        )
    return stores


def carbon(n: pypsa.Network, config: dict, planning_horizons: str | None) -> bool:
    """
    Bound the emissions accumulated in the emission stores.

    A `co2_atmosphere` row bounds the emissions of its investment period, or
    of the horizon without one. A `Co2Budget` row bounds the emissions
    accumulated until the end of its investment period, times the years of
    the period.
    """
    glcs = n.global_constraints.query("type in ['co2_atmosphere', 'Co2Budget']")
    if glcs.empty:
        return False
    attributes = glcs.carrier_attribute.unique()
    if len(attributes) > 1:
        raise ValueError(
            f"The CO2 limits {list(glcs.index)} name several carrier attributes "
            f"{list(attributes)}, but the emission stores they bound must be the same."
        )
    emissions = n.carriers.index[n.carriers[attributes[0]] != 0]
    stores = emission_stores(n, emissions)
    if stores.empty:
        return False
    n.stores["emission_store"] = n.stores.index.isin(stores.index)
    return True


def add_max_growth(n: pypsa.Network, opts: dict) -> None:
    """
    Add maximum growth rates for different carriers.
    """

    # take maximum yearly difference between investment periods since historic growth is per year
    factor = n.investment_period_weightings.years.max() * opts["factor"]
    for carrier in opts["max_growth"]:
        max_per_period = opts["max_growth"][carrier] * factor
        logger.info(
            f"set maximum growth rate per investment period of {carrier} to {max_per_period} GW."
        )
        n.carriers.loc[carrier, "max_growth"] = max_per_period * 1e3

    for carrier in opts["max_relative_growth"]:
        max_r_per_period = opts["max_relative_growth"][carrier]
        logger.info(
            f"set maximum relative growth per investment period of {carrier} to {max_r_per_period}."
        )
        n.carriers.loc[carrier, "max_relative_growth"] = max_r_per_period


def gas_boiler_retrofit(
    n: pypsa.Network, config: dict, planning_horizons: str | None
) -> bool:
    """
    Let existing gas boilers be retrofitted to H2 boilers that follow the heat demand.

    Each existing gas boiler and its retrofitted H2 boiler together follow
    the heat demand at the boiler's heat bus, scaled to the gas boiler's
    capacity at peak demand. Only with investment periods.
    """
    if n.investment_periods.empty:
        return False
    links = n.links
    gas_i = links.index[
        links.carrier.str.contains("gas boiler") & ~links.p_nom_extendable
    ]
    h2_i = gas_i.str.replace("gas boiler", "retrofitted H2 boiler")
    active = n.c.links.active_assets
    paired = h2_i.isin(active) & gas_i.isin(active)
    if not paired.any():
        return False

    logger.info("Add constraint for retrofitting gas boilers to H2 boilers.")
    gas_i, h2_i = gas_i[paired], h2_i[paired]
    p_set = n.loads_t.p_set
    heat = p_set.columns[
        p_set.columns.str.contains("heat")
        & ~p_set.columns.str.contains("industry")
        & ~p_set.columns.str.contains("agriculture")
    ]
    profile = p_set[heat].div(p_set[heat].groupby(level=0).max(), level=0).fillna(0)
    profile = profile.rename(columns=n.loads.bus).reindex(columns=links.bus1[gas_i])
    n.links_t["heat_demand_pu"] = profile.set_axis(gas_i, axis=1)
    n.links["H2_boiler_retrofit"] = links.index.map(dict(zip(gas_i, h2_i)))
    return True


def enforce_autarky(n: pypsa.Network, only_crossborder: bool = False) -> None:
    """
    Remove transmission lines to enforce autarky.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    only_crossborder : bool
        If True, only remove cross-border connections. If False, remove all lines.
    """
    if only_crossborder:
        lines_rm = n.lines.loc[
            n.lines.bus0.map(n.buses.country) != n.lines.bus1.map(n.buses.country)
        ].index
        links_rm = n.links.loc[
            n.links.bus0.map(n.buses.country) != n.links.bus1.map(n.buses.country)
        ].index
    else:
        lines_rm = n.lines.index
        links_rm = n.links.loc[n.links.carrier == "DC"].index

    logger.info(
        f"Enforcing autarky: removing {len(lines_rm)} lines and {len(links_rm)} links"
    )
    n.remove("Line", lines_rm)
    n.remove("Link", links_rm)


def add_load_balance_components(
    n: pypsa.Network, config: dict, sign: float = 1
) -> None:
    """
    Add load shedding or load sinks to the network with carrier 'load'.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to be modified.
    config : dict
        The load shedding or load sinks settings.
    sign : float
        Direction of the added generators. Positive for load shedding, negative for load sinks.

    Returns
    -------
    None
        Modifies PyPSA network in place.
    """
    if "load" not in n.carriers.index:
        n.add("Carrier", "load")

    carriers = config.get("carriers", {})
    default_cost = config.get("default_cost")
    balance_comp = "shedding" if sign > 0 else "sink"

    logger.info(
        f"Add load {balance_comp} for {'all carriers' if config.get('all_carriers') else ', '.join(carriers)}."
    )

    for bus_carrier, price in carriers.items():
        buses_i = n.buses[n.buses.carrier == bus_carrier].index
        n.add(
            "Generator",
            buses_i,
            f" load {balance_comp}",
            bus=buses_i,
            carrier="load",
            marginal_cost=price,
            p_nom=np.inf,
            sign=sign,
        )

    if config.get("all_carriers", False):
        buses_rest_i = n.buses[~n.buses.carrier.isin(carriers)].index
        n.add(
            "Generator",
            buses_rest_i,
            f" load {balance_comp}",
            bus=buses_rest_i,
            carrier="load",
            marginal_cost=default_cost,
            p_nom=np.inf,
            sign=sign,
        )


def prepare_network(
    n: pypsa.Network,
    config: dict,
    solve_opts: dict,
    foresight: str,
    planning_horizons: str | None,
    co2_sequestration_potential: dict[str, float],
    limit_max_growth: dict[str, Any] | None = None,
    rolling_horizon: bool = False,
) -> None:
    """
    Prepare network with various constraints and modifications.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    config : dict
        Configuration dictionary
    solve_opts : Dict
        Dictionary of solving options containing clip_p_max_pu, load_shedding etc.
    foresight : str
        Planning foresight type ('myopic' or 'perfect')
    planning_horizons : str or None
        The current planning horizon year or None for perfect foresight
    co2_sequestration_potential : Dict[str, float]
        CO2 sequestration potential constraints by year

    Returns
    -------
    pypsa.Network
        Modified PyPSA network with added constraints
    """
    if "clip_p_max_pu" in solve_opts:
        for df in (
            n.generators_t.p_max_pu,
            n.generators_t.p_min_pu,
            n.links_t.p_max_pu,
            n.links_t.p_min_pu,
            n.storage_units_t.inflow,
        ):
            df.where(df.abs() > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    if (load_shedding := solve_opts.get("load_shedding", {})).get("enable", False):
        # intersect between macroeconomic and surveybased willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
        add_load_balance_components(n, load_shedding)

    if (load_sinks := solve_opts.get("load_sinks", {})).get("enable", False):
        add_load_balance_components(n, load_sinks, sign=-1)

    if solve_opts.get("curtailment_mode"):
        n.add("Carrier", "curtailment", color="#fedfed", nice_name="Curtailment")
        n.generators_t.p_min_pu = n.generators_t.p_max_pu
        buses_i = n.buses.query("carrier == 'AC'").index
        n.add(
            "Generator",
            buses_i,
            suffix=" curtailment",
            bus=buses_i,
            p_min_pu=-1,
            p_max_pu=0,
            marginal_cost=-0.1,
            carrier="curtailment",
            p_nom=1e6,
        )

    if solve_opts.get("noisy_costs"):
        for t in n.components:
            # if 'capital_cost' in t.static:
            #    t.static['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.static)) - 0.5)
            if "marginal_cost" in t.static:
                t.static["marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.static)) - 0.5
                )

        for t in n.components[["Line", "Link"]]:
            if t.static.empty:
                continue
            t.static["capital_cost"] += (
                1e-1 + 2e-2 * (np.random.random(len(t.static)) - 0.5)
            ) * t.static["length"]

    if solve_opts.get("nhours"):
        nhours = solve_opts["nhours"]
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760.0 / nhours

    if foresight == "perfect":
        add_land_use_constraint(n)
        if limit_max_growth is not None and limit_max_growth["enable"]:
            add_max_growth(n, limit_max_growth)

    if n.stores.carrier.eq("co2 sequestered").any():
        limit_dict = co2_sequestration_potential
        add_co2_sequestration_limit(
            n, limit_dict=limit_dict, planning_horizons=planning_horizons
        )

    if config["solving"]["constraints"]["BAU"] and n.generators.p_nom_extendable.any():
        add_BAU_constraints(n, config)

    # Apply autarky constraint if configured
    autarky_cfg = config["electricity"]["autarky"]
    if autarky_cfg["enable"]:
        only_crossborder = autarky_cfg["by_country"]
        enforce_autarky(n, only_crossborder=only_crossborder)

    # rolling horizon disables cyclic storage
    if rolling_horizon:
        n.storage_units.cyclic_state_of_charge = False
        n.storage_units.state_of_charge_initial = 0
        n.stores.e_cyclic = False
        n.stores.e_initial = 0

    register_math(n, config, planning_horizons)


def ccl(n: pypsa.Network, config: dict, planning_horizons: str | None) -> bool:
    """
    Bound the generator capacity per country and carrier (CCL).

    Add minimum and maximum levels of generator nominal capacity per carrier
    for individual countries. Opts and path for agg_p_nom_minmax.csv must be defined
    in config.yaml. Default file is available at data/agg_p_nom_minmax.csv.

    Examples
    --------
    scenario:
        opts: [Co2L-CCL-24h]
    electricity:
        agg_p_nom_limits: data/agg_p_nom_minmax.csv
    """
    gens = n.generators
    if not (config["solving"]["constraints"]["CCL"] and gens.p_nom_extendable.any()):
        return False
    assert planning_horizons is not None, (
        "The CCL constraints are not implemented for perfect foresight, yet"
    )

    options = config["solving"]["agg_p_nom_limits"]
    limits = pd.read_csv(options["file"], index_col=[0, 1], header=[0, 1])
    limits = limits[planning_horizons].set_axis(limits.index.map(" ".join))

    rename = {}
    if options["agg_offwind"]:
        offwind = ["offwind-ac", "offwind-dc", "offwind-float", "offwind"]
        rename |= dict.fromkeys(offwind, "offwind-all")
    if options["agg_solar"]:
        rename |= dict.fromkeys(["solar", "solar-hsat", "solar rooftop"], "solar-all")
    group = gens.bus.map(n.buses.country) + " " + gens.carrier.replace(rename)

    if options["include_existing"]:
        existing = gens[~gens.p_nom_extendable]
        existing = existing[
            existing.build_year + existing.lifetime >= int(planning_horizons)
        ]
        existing = existing.p_nom.groupby(group).sum().reindex(limits.index)
        limits = limits.sub(existing.fillna(0), axis=0).clip(lower=0)

    extendable = gens.index[gens.p_nom_extendable].intersection(
        n.c.generators.active_assets
    )
    limits = limits[limits.index.isin(group[extendable])]
    if limits.empty:
        return False

    logger.info("Adding generation capacity constraints per carrier and country")
    for bound, sense in [("min", ">="), ("max", "<=")]:
        rhs = limits[bound].dropna()
        names = f"CCL {bound} " + rhs.index
        n.add(
            "GlobalConstraint",
            names,
            type="CCL",
            sense=sense,
            constant=rhs.values,
            overwrite=True,
        )
        row = pd.Series(names, index=rhs.index)
        n.generators[f"CCL_{bound}"] = group.map(row).where(gens.index.isin(extendable))
    return True


def equity(n: pypsa.Network, config: dict, planning_horizons: str | None) -> bool:
    """
    Require each country or node to produce a share of its demand itself.

    Currently this is only implemented for the electricity sector.

    Examples
    --------
    scenario:
        opts: [Co2L-EQ0.7-24h]

    EQ0.7c demands each country to produce on average at least 70% of its
    consumption; EQ0.7 demands each node to produce on average at least 70%
    of its consumption.
    """
    if not (o := config["solving"]["constraints"]["EQ"]):
        return False
    level = float(re.findall(r"[0-9]*\.?[0-9]+", o)[0])

    def region(bus: pd.Series) -> pd.Series:
        return bus.map(n.buses.country) if o[-1] == "c" else bus

    inflow = n.storage_units.index.isin(n.storage_units_t.inflow.columns)
    n.generators["equity_region"] = region(n.generators.bus)
    n.storage_units["equity_region"] = region(n.storage_units.bus).where(inflow)
    n.loads["equity_region"] = region(n.loads.bus)
    n.loads["equity_level"] = level
    return True


def add_BAU_constraints(n: pypsa.Network, config: dict) -> None:
    """
    Add business-as-usual (BAU) minimum capacities per carrier.

    Each minimum is a global constraint of type `tech_capacity_expansion_limit`,
    which bounds the extendable capacity of all components of the carrier.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network instance
    config : dict
        Configuration dictionary containing BAU minimum capacities
    """
    mincaps = pd.Series(config["electricity"]["BAU_mincapacities"])
    n.add(
        "GlobalConstraint",
        "BAU minimum " + mincaps.index,
        type="tech_capacity_expansion_limit",
        carrier_attribute=mincaps.index,
        sense=">=",
        constant=mincaps.values,
    )


# TODO: think about removing or make per country
def safe(n: pypsa.Network, config: dict, planning_horizons: str | None) -> bool:
    """
    Add a capacity reserve margin of a certain fraction above the peak demand.

    Only conventional generators count. Renewable generators and storage do
    not contribute. Ignores the network.

    Examples
    --------
    scenario:
        opts: [Co2L-SAFE-24h]
    electricity:
        SAFE_reservemargin: 0.1
    Which sets a reserve margin of 10% above the peak demand.
    """
    gens = n.generators
    if not (config["solving"]["constraints"]["SAFE"] and gens.p_nom_extendable.any()):
        return False
    peak = n.loads_t.p_set.sum(axis=1).max()
    margin = 1.0 + config["electricity"]["SAFE_reservemargin"]
    conventional = gens.carrier.isin(config["electricity"]["conventional_carriers"])
    n.add(
        "GlobalConstraint",
        "SAFE",
        type="SAFE",
        sense=">=",
        constant=peak * margin,
        overwrite=True,
    )
    n.generators["SAFE"] = pd.Series("SAFE", gens.index).where(conventional)
    return True


def operational_reserve(
    n: pypsa.Network, config: dict, planning_horizons: str | None
) -> bool:
    """
    Build reserve margin constraints based on the formulation given in
    https://genxproject.github.io/GenX.jl/stable/Model_Reference/core/#Operational-Reserves.

    Examples
    --------
    electricity:
        operational_reserve:
            activate: true
            epsilon_load: 0.02 # percentage of load at each snapshot
            epsilon_vres: 0.02 # percentage of VRES at each snapshot
            contingency: 4000 # MW
    """
    reserve = config["electricity"].get("operational_reserve", {})
    if not reserve.get("activate"):
        return False
    if config["sector"]["enabled"]:
        logger.warning(
            "Operational reserve margin only considers electricity generators "
            "and ignores sector-coupling components."
        )
    vres = [
        "ror" if c == "hydro" else c
        for c in config["electricity"]["renewable_carriers"]
    ]
    carrier = n.generators.carrier
    n.add(
        "GlobalConstraint",
        "operational_reserve",
        type="operational_reserve",
        sense=">=",
        constant=reserve["contingency"],
        overwrite=True,
    )
    n.generators["operational_reserve"] = carrier != "load"
    n.generators["operational_reserve_share"] = np.where(
        carrier.isin(vres), reserve["epsilon_vres"], np.nan
    )
    n.loads["operational_reserve_share"] = reserve["epsilon_load"]
    return True


def reversed_links(n: pypsa.Network) -> pd.Series:
    """Return whether each link is the reversed link of a lossy bidirectional connection."""
    return n.links.get("reversed", pd.Series(False, n.links.index)).astype(bool)


def lossy_bidirectional_links(
    n: pypsa.Network, config: dict, planning_horizons: str | None
) -> bool:
    """
    Tie the capacity of each extendable reversed link to its forward link.
    """
    links = n.links
    extendable = links.index[links.p_nom_extendable].intersection(
        n.c.links.active_assets
    )
    backwards = extendable[reversed_links(n)[extendable]]
    forwards = backwards.str.replace("-reversed", "")
    paired = forwards.isin(extendable)
    if not paired.any():
        return False
    n.links["reversal_of"] = links.index.map(
        dict(zip(backwards[paired], forwards[paired]))
    )
    return True


def pipe_retrofit(
    n: pypsa.Network, config: dict, planning_horizons: str | None
) -> bool:
    """
    Let extendable gas pipelines be retrofitted to hydrogen pipelines.
    """
    links = n.links
    candidates = links.index[links.p_nom_extendable & ~reversed_links(n)]
    candidates = candidates.intersection(n.c.links.active_assets)
    carrier = links.carrier[candidates]
    gas_i = candidates[carrier.eq("gas pipeline")]
    h2_i = gas_i.str.replace("gas pipeline", "H2 pipeline retrofitted")
    paired = h2_i.isin(candidates[carrier.eq("H2 pipeline retrofitted")])
    if not paired.any():
        return False
    retrofit = links.index.map(dict(zip(gas_i[paired], h2_i[paired])))
    n.links["H2_pipe_retrofit"] = retrofit
    capacity_per_CH4 = config["sector"]["H2_retrofit_capacity_per_CH4"]
    n.links["H2_retrofit_capacity_per_CH4"] = np.where(
        retrofit.notna(), capacity_per_CH4, np.nan
    )
    return True


def flexible_egs(n: pypsa.Network, config: dict, planning_horizons: str | None) -> bool:
    """
    Limit each geothermal reservoir to the capacity of the wells at its surface bus.
    """
    egs = config["sector"]["enhanced_geothermal"]
    if not (egs["enable"] and egs["flexible"]):
        return False
    wells = n.links.carrier.eq("geothermal heat")
    reservoirs = n.storage_units.carrier.eq("geothermal heat")
    reservoirs &= n.storage_units.p_nom_extendable
    n.links["geothermal_well"] = n.links.bus1.where(wells)
    n.storage_units["geothermal_reservoir"] = n.storage_units.bus.where(reservoirs)
    return True


def import_limit(n: pypsa.Network, config: dict, planning_horizons: str | None) -> bool:
    """
    Limit green energy imports (synthetic and biomass).

    Does not include fossil fuel imports. Links count with the fuel they deliver.
    """
    imports = config["sector"]["imports"]
    if not imports["enable"]:
        return False
    gens = n.generators.carrier.str.contains("import")
    links = n.links.carrier.str.contains("import")
    if not (gens.any() or links.any()) or not np.isfinite(imports["limit"]):
        return False

    nyears = n.snapshot_weightings.generators.sum() / 8760
    n.add(
        "GlobalConstraint",
        "import_limit",
        type="import_limit",
        sense=imports["limit_sense"],
        constant=imports["limit"] * 1e6 * nyears,
        overwrite=True,
    )
    n.generators["import_limit"] = pd.Series("import_limit", gens.index).where(gens)
    n.links["import_limit"] = pd.Series("import_limit", links.index).where(links)
    n.links["import_efficiency"] = n.links.efficiency.where(links)
    return True


MathWriter = Callable[[pypsa.Network, dict, str | None], bool]

MATH_WRITERS: tuple[MathWriter, ...] = (
    lossy_bidirectional_links,
    pipe_retrofit,
    carbon,
    safe,
    ccl,
    equity,
    solar_potential,
    operational_reserve,
    gas_boiler_retrofit,
    flexible_egs,
    import_limit,
)


def register_math(
    n: pypsa.Network, config: dict, planning_horizons: str | None
) -> None:
    """
    Store the custom math of the active constraints on the network.

    Each writer in `MATH_WRITERS` checks its switch, writes the data its math
    in `scripts/math/<writer>.yaml` reads onto the network and returns whether
    the math applies. The network stores the math as overrides, so every build
    applies it, also in iterative and rolling-horizon solves and after a
    reload.
    """
    overrides = n.optimize.overrides
    for writer in MATH_WRITERS:
        name = writer.__name__
        if name in overrides:
            overrides.unregister(name)
        if writer(n, config, planning_horizons):
            logger.info(f"Adding custom math {name}.")
            overrides.register(MATH / f"{name}.yaml")


def extra_functionality(n: pypsa.Network, snapshots: pd.Index, snakemake=None) -> None:
    """
    Add the constraints of `custom_extra_functionality` to the model.

    The constraints of PyPSA-Eur are custom math in `scripts/math`, which
    `prepare_network` stores on the network. If you want to enforce
    additional custom constraints, add math there, or point
    `solving: custom_extra_functionality` to a function with linopy code.
    """
    if source_path := n.params.custom_extra_functionality:
        assert os.path.exists(source_path), f"{source_path} does not exist"
        sys.path.append(os.path.dirname(source_path))
        module_name = os.path.splitext(os.path.basename(source_path))[0]
        module = importlib.import_module(module_name)
        custom_extra_functionality = getattr(module, module_name)
        custom_extra_functionality(n, snapshots, snakemake)


def check_objective_value(n: pypsa.Network, solving: dict, horizon: str) -> None:
    """
    Check if total system cost matches expected value within tolerance.

    The total system cost is the sum of capital and operational expenditure
    from `n.statistics`.

    Parameters
    ----------
    n : pypsa.Network
        Solved network
    solving : Dict
        Dictionary containing objective checking parameters
    horizon : str
        Planning horizon of the current solve, used to select the expected
        value for myopic foresight (per-horizon mapping)

    Raises
    ------
    ObjectiveValueError
        If total system cost differs from expected value beyond tolerance
    """
    check_objective = solving["check_objective"]
    if not check_objective["enable"]:
        return
    atol = check_objective["atol"]
    rtol = check_objective["rtol"]
    expected_value = check_objective["expected_value"]
    if isinstance(expected_value, dict):
        expected_value = expected_value[int(horizon)]
    system_cost = np.nansum(n.statistics.capex()) + np.nansum(n.statistics.opex())
    if not np.isclose(system_cost, expected_value, atol=atol, rtol=rtol):
        raise ObjectiveValueError(
            f"Total system cost {system_cost} differs from expected value "
            f"{expected_value} by more than {atol}."
        )


def collect_kwargs(
    config: dict,
    solving: dict,
    planning_horizons: str | None = None,
    log_fn: str | None = None,
    mode: str = "single",
    horizon: int | None = None,
    overlap: int | None = None,
) -> tuple[dict, dict]:
    """
    Prepare keyword arguments separated for model creation and model solving.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing solver settings
    solving : dict
        Dictionary of solving options and configuration
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight
    log_fn : str, optional
        Path to solver log file
    mode : str, optional
        Optimization mode: 'single', 'rolling_horizon', or 'iterative'
        Default is 'single'

    Returns
    -------
    tuple[dict, dict]
        Two dictionaries: (model_kwargs, solve_kwargs)
        - model_kwargs: Arguments for n.optimize.create_model()
        - solve_kwargs: Arguments for n.optimize.solve_model()
        For 'rolling_horizon' and 'iterative' modes, returns merged kwargs
        with additional mode-specific parameters
    """
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    # Model creation kwargs
    model_kwargs = {}
    model_kwargs["multi_investment_periods"] = config["foresight"] == "perfect"
    model_kwargs["transmission_losses"] = cf_solving.get("transmission_losses", False)
    model_kwargs["linearized_unit_commitment"] = cf_solving.get(
        "linearized_unit_commitment", False
    )
    model_kwargs["include_objective_constant"] = False

    # Solve kwargs
    solver_name = solving["solver"]["name"]
    solver_options = solving["solver_options"][set_of_options] if set_of_options else {}

    solve_kwargs = {}
    solve_kwargs["solver_name"] = solver_name
    solve_kwargs["solver_options"] = solver_options
    solve_kwargs["assign_all_duals"] = cf_solving.get("assign_all_duals", False)
    solve_kwargs["io_api"] = cf_solving.get("io_api", None)
    solve_kwargs["keep_files"] = cf_solving.get("keep_files", False)

    if log_fn:
        solve_kwargs["log_fn"] = log_fn

    if solver_name == "gurobi":
        logging.getLogger("gurobipy").setLevel(logging.CRITICAL)

    # Handle special modes
    if mode == "rolling_horizon":
        all_kwargs = {**model_kwargs, **solve_kwargs}
        all_kwargs["horizon"] = (
            horizon if horizon is not None else cf_solving["horizon"]
        )
        all_kwargs["overlap"] = (
            overlap if overlap is not None else cf_solving["overlap"]
        )
        return all_kwargs, {}

    elif mode == "iterative":
        all_kwargs = {**model_kwargs, **solve_kwargs}
        all_kwargs["track_iterations"] = cf_solving["track_iterations"]
        all_kwargs["min_iterations"] = cf_solving["min_iterations"]
        all_kwargs["max_iterations"] = cf_solving["max_iterations"]

        if cf_solving["post_discretization"].get("enable", False):
            logger.info("Add post-discretization parameters.")
            cf_solving["post_discretization"].pop("enable", None)
            all_kwargs.update(cf_solving["post_discretization"])

        return all_kwargs, {}

    return model_kwargs, solve_kwargs


def create_optimization_model(
    n: pypsa.Network,
    config: dict,
    params: dict,
    model_kwargs: dict,
    solve_kwargs: dict,
    snakemake=None,
) -> None:
    """
    Prepare optimization problem by creating model and adding extra functionality.

    This function:
    1. Attaches config and params to network for extra_functionality
    2. Creates the optimization model
    3. Adds extra functionality (custom constraints)

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    config : dict
        Configuration dictionary containing solver settings
    params : dict
        Dictionary of solving parameters
    model_kwargs : dict
        Arguments for n.optimize.create_model()
    solve_kwargs : dict
        Arguments for n.optimize.solve_model()
    """
    # Add config and params to network for extra_functionality
    n.config = config
    n.params = params

    # Create optimization model
    logger.info("Creating optimization model...")
    n.optimize.create_model(**model_kwargs)

    # Add extra functionality (custom constraints)
    logger.info("Adding extra functionality (custom constraints)...")
    extra_functionality(n, n.snapshots, snakemake=snakemake)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_network",
            configfiles="config/test/config.overnight.yaml",
            horizon=2030,
        )
    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    solve_opts = snakemake.params.solving["options"]
    cf_solving = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    # Load network
    n = pypsa.Network(snakemake.input.network)
    planning_horizons = snakemake.wildcards.horizon

    # Prepare network (settings before solving)
    prepare_network(
        n,
        config=snakemake.config,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params["sector"]["limit_max_growth"],
        rolling_horizon=cf_solving["rolling_horizon"],
    )

    # Determine solve mode
    rolling_horizon = cf_solving["rolling_horizon"]
    skip_iterations = cf_solving.get("skip_iterations", False)

    if not n.lines.s_nom_extendable.any():
        skip_iterations = True
        logger.info("No expandable lines found. Skipping iterative solving.")

    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )

    # Solve network based on mode
    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=logging_frequency
    ) as mem:
        if rolling_horizon:
            logger.info("Using rolling horizon optimization...")
            all_kwargs, _ = collect_kwargs(
                snakemake.config,
                snakemake.params.solving,
                planning_horizons,
                log_fn=snakemake.log.solver,
                mode="rolling_horizon",
            )

            n.config = snakemake.config
            n.params = snakemake.params
            all_kwargs["extra_functionality"] = partial(
                extra_functionality,
                snakemake=snakemake,
            )
            n.optimize.optimize_with_rolling_horizon(**all_kwargs)
            status, condition = "", ""

        elif skip_iterations:
            logger.info("Using single-pass optimization...")
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
                snakemake=snakemake,
            )

            logger.info("Solving model...")
            status, condition = n.optimize.solve_model(**solve_kwargs)

        else:
            logger.info("Using iterative transmission expansion optimization...")

            all_kwargs, _ = collect_kwargs(
                snakemake.config,
                snakemake.params.solving,
                planning_horizons,
                log_fn=snakemake.log.solver,
                mode="iterative",
            )

            n.config = snakemake.config
            n.params = snakemake.params
            all_kwargs["extra_functionality"] = partial(
                extra_functionality,
                snakemake=snakemake,
            )
            status, condition = n.optimize.optimize_transmission_expansion_iteratively(
                **all_kwargs
            )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    # Check results
    if not rolling_horizon:
        if status != SolverStatus.ok:
            logger.warning(
                f"Solving status '{status}' with termination condition '{condition}'"
            )
        check_objective_value(n, snakemake.params.solving, snakemake.wildcards.horizon)

    if condition in [
        TerminationCondition.infeasible,
        TerminationCondition.infeasible_or_unbounded,
    ]:
        labels = n.model.compute_infeasibilities()
        logger.info(f"Labels:\n{labels}")
        n.model.print_infeasibilities()
        raise RuntimeError("Solving status 'infeasible'. Infeasibilities computed.")

    if status == SolverStatus.warning:
        raise RuntimeError("Solving status 'warning'. Discarding solution.")

    n.meta["snakemake_config"] = dict(
        snakemake.config, wildcards=dict(snakemake.wildcards)
    )
    n.export_to_netcdf(snakemake.output.network)

    if snakemake.output.get("model"):
        n.model.to_netcdf(snakemake.output.model)

    config_output = getattr(snakemake.output, "config", None)
    if config_output:
        with open(config_output, "w") as file:
            yaml.dump(
                n.meta["snakemake_config"],
                file,
                default_flow_style=False,
                allow_unicode=True,
                sort_keys=False,
            )
