# -*- coding: utf-8 -*-
<<<<<<< HEAD
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Creates summaries of aggregated energy and costs as ``.csv`` files.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        year:
        version:
        fill_values:
        marginal_cost:
        capital_cost:

    electricity:
        max_hours:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`costs_cf`, :ref:`electricity_cf`

Inputs
------

Outputs
-------

Description
-----------

The following rule can be used to summarize the results in separate .csv files:

.. code:: bash

    snakemake results/summaries/elec_s_all_lall_Co2L-3H_all
                                         clusters
                                             line volume or cost cap
                                                - options
                                                        - all countries

the line volume/cost cap field can be set to one of the following:
* ``lv1.25`` for a particular line volume extension by 25%
* ``lc1.25`` for a line cost extension by 25 %
* ``lall`` for all evaluated caps
* ``lvall`` for all line volume caps
* ``lcall`` for all line cost caps

Replacing '/summaries/' with '/plots/' creates nice colored maps of the results.
"""

import logging
import os

import pandas as pd
import pypsa
from _helpers import configure_logging
from add_electricity import load_costs, update_transmission_costs

idx = pd.IndexSlice

logger = logging.getLogger(__name__)

opt_name = {"Store": "e", "Line": "s", "Transformer": "s"}


def _add_indexed_rows(df, raw_index):
    new_index = df.index.union(pd.MultiIndex.from_product(raw_index))
    if isinstance(new_index, pd.Index):
        new_index = pd.MultiIndex.from_tuples(new_index)

    return df.reindex(new_index)


def assign_carriers(n):
    if "carrier" not in n.loads:
        n.loads["carrier"] = "electricity"
        for carrier in ["transport", "heat", "urban heat"]:
            n.loads.loc[n.loads.index.str.contains(carrier), "carrier"] = carrier

    n.storage_units["carrier"].replace(
        {"hydro": "hydro+PHS", "PHS": "hydro+PHS"}, inplace=True
    )

    if "carrier" not in n.lines:
        n.lines["carrier"] = "AC"

    n.lines["carrier"].replace({"AC": "lines"}, inplace=True)

    if n.links.empty:
        n.links["carrier"] = pd.Series(dtype=str)
    n.links["carrier"].replace({"DC": "lines"}, inplace=True)

    if (
        "EU gas store" in n.stores.index
        and n.stores.loc["EU gas Store", "carrier"] == ""
    ):
        n.stores.loc["EU gas Store", "carrier"] = "gas Store"


def calculate_costs(n, label, costs):
    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        capital_costs = c.df.capital_cost * c.df[opt_name.get(c.name, "p") + "_nom_opt"]
        capital_costs_grouped = capital_costs.groupby(c.df.carrier).sum()

        # Index tuple(s) indicating the newly to-be-added row(s)
        raw_index = tuple(
            [[c.list_name], ["capital"], list(capital_costs_grouped.index)]
        )
        costs = _add_indexed_rows(costs, raw_index)

        costs.loc[idx[raw_index], label] = capital_costs_grouped.values
=======
import logging

logger = logging.getLogger(__name__)

import sys

import numpy as np
import pandas as pd
import pypsa
import yaml
from helper import override_component_attrs
from prepare_sector_network import prepare_costs

idx = pd.IndexSlice

opt_name = {"Store": "e", "Line": "s", "Transformer": "s"}


def assign_carriers(n):
    if "carrier" not in n.lines:
        n.lines["carrier"] = "AC"


def assign_locations(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.unique():
            names = ifind.index[ifind == i]
            if i == -1:
                c.df.loc[names, "location"] = ""
            else:
                c.df.loc[names, "location"] = names.str[:i]


def calculate_nodal_cfs(n, label, nodal_cfs):
    # Beware this also has extraneous locations for country (e.g. biomass) or continent-wide (e.g. fossil gas/oil) stuff
    for c in n.iterate_components(
        (n.branch_components ^ {"Line", "Transformer"})
        | n.controllable_one_port_components ^ {"Load", "StorageUnit"}
    ):
        capacities_c = c.df.groupby(["location", "carrier"])[
            opt_name.get(c.name, "p") + "_nom_opt"
        ].sum()

        if c.name == "Link":
            p = c.pnl.p0.abs().mean()
        elif c.name == "Generator":
            p = c.pnl.p.abs().mean()
        elif c.name == "Store":
            p = c.pnl.e.abs().mean()
        else:
            sys.exit()

        c.df["p"] = p
        p_c = c.df.groupby(["location", "carrier"])["p"].sum()

        cf_c = p_c / capacities_c

        index = pd.MultiIndex.from_tuples(
            [(c.list_name,) + t for t in cf_c.index.to_list()]
        )
        nodal_cfs = nodal_cfs.reindex(index.union(nodal_cfs.index))
        nodal_cfs.loc[index, label] = cf_c.values

    return nodal_cfs


def calculate_cfs(n, label, cfs):
    for c in n.iterate_components(
        n.branch_components
        | n.controllable_one_port_components ^ {"Load", "StorageUnit"}
    ):
        capacities_c = (
            c.df[opt_name.get(c.name, "p") + "_nom_opt"].groupby(c.df.carrier).sum()
        )

        if c.name in ["Link", "Line", "Transformer"]:
            p = c.pnl.p0.abs().mean()
        elif c.name == "Store":
            p = c.pnl.e.abs().mean()
        else:
            p = c.pnl.p.abs().mean()

        p_c = p.groupby(c.df.carrier).sum()

        cf_c = p_c / capacities_c

        cf_c = pd.concat([cf_c], keys=[c.list_name])

        cfs = cfs.reindex(cf_c.index.union(cfs.index))

        cfs.loc[cf_c.index, label] = cf_c

    return cfs


def calculate_nodal_costs(n, label, nodal_costs):
    # Beware this also has extraneous locations for country (e.g. biomass) or continent-wide (e.g. fossil gas/oil) stuff
    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        c.df["capital_costs"] = (
            c.df.capital_cost * c.df[opt_name.get(c.name, "p") + "_nom_opt"]
        )
        capital_costs = c.df.groupby(["location", "carrier"])["capital_costs"].sum()
        index = pd.MultiIndex.from_tuples(
            [(c.list_name, "capital") + t for t in capital_costs.index.to_list()]
        )
        nodal_costs = nodal_costs.reindex(index.union(nodal_costs.index))
        nodal_costs.loc[index, label] = capital_costs.values
>>>>>>> pypsa-eur-sec/master

        if c.name == "Link":
            p = c.pnl.p0.multiply(n.snapshot_weightings.generators, axis=0).sum()
        elif c.name == "Line":
            continue
        elif c.name == "StorageUnit":
            p_all = c.pnl.p.multiply(n.snapshot_weightings.generators, axis=0)
            p_all[p_all < 0.0] = 0.0
            p = p_all.sum()
        else:
            p = c.pnl.p.multiply(n.snapshot_weightings.generators, axis=0).sum()

<<<<<<< HEAD
=======
        # correct sequestration cost
        if c.name == "Store":
            items = c.df.index[
                (c.df.carrier == "co2 stored") & (c.df.marginal_cost <= -100.0)
            ]
            c.df.loc[items, "marginal_cost"] = -20.0

        c.df["marginal_costs"] = p * c.df.marginal_cost
        marginal_costs = c.df.groupby(["location", "carrier"])["marginal_costs"].sum()
        index = pd.MultiIndex.from_tuples(
            [(c.list_name, "marginal") + t for t in marginal_costs.index.to_list()]
        )
        nodal_costs = nodal_costs.reindex(index.union(nodal_costs.index))
        nodal_costs.loc[index, label] = marginal_costs.values

    return nodal_costs


def calculate_costs(n, label, costs):
    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        capital_costs = c.df.capital_cost * c.df[opt_name.get(c.name, "p") + "_nom_opt"]
        capital_costs_grouped = capital_costs.groupby(c.df.carrier).sum()

        capital_costs_grouped = pd.concat([capital_costs_grouped], keys=["capital"])
        capital_costs_grouped = pd.concat([capital_costs_grouped], keys=[c.list_name])

        costs = costs.reindex(capital_costs_grouped.index.union(costs.index))

        costs.loc[capital_costs_grouped.index, label] = capital_costs_grouped

        if c.name == "Link":
            p = c.pnl.p0.multiply(n.snapshot_weightings.generators, axis=0).sum()
        elif c.name == "Line":
            continue
        elif c.name == "StorageUnit":
            p_all = c.pnl.p.multiply(n.snapshot_weightings.generators, axis=0)
            p_all[p_all < 0.0] = 0.0
            p = p_all.sum()
        else:
            p = c.pnl.p.multiply(n.snapshot_weightings.generators, axis=0).sum()

        # correct sequestration cost
        if c.name == "Store":
            items = c.df.index[
                (c.df.carrier == "co2 stored") & (c.df.marginal_cost <= -100.0)
            ]
            c.df.loc[items, "marginal_cost"] = -20.0

>>>>>>> pypsa-eur-sec/master
        marginal_costs = p * c.df.marginal_cost

        marginal_costs_grouped = marginal_costs.groupby(c.df.carrier).sum()

<<<<<<< HEAD
        costs = costs.reindex(
            costs.index.union(
                pd.MultiIndex.from_product(
                    [[c.list_name], ["marginal"], marginal_costs_grouped.index]
                )
            )
        )

        costs.loc[
            idx[c.list_name, "marginal", list(marginal_costs_grouped.index)], label
        ] = marginal_costs_grouped.values
=======
        marginal_costs_grouped = pd.concat([marginal_costs_grouped], keys=["marginal"])
        marginal_costs_grouped = pd.concat([marginal_costs_grouped], keys=[c.list_name])

        costs = costs.reindex(marginal_costs_grouped.index.union(costs.index))

        costs.loc[marginal_costs_grouped.index, label] = marginal_costs_grouped

    # add back in all hydro
    # costs.loc[("storage_units", "capital", "hydro"),label] = (0.01)*2e6*n.storage_units.loc[n.storage_units.group=="hydro", "p_nom"].sum()
    # costs.loc[("storage_units", "capital", "PHS"),label] = (0.01)*2e6*n.storage_units.loc[n.storage_units.group=="PHS", "p_nom"].sum()
    # costs.loc[("generators", "capital", "ror"),label] = (0.02)*3e6*n.generators.loc[n.generators.group=="ror", "p_nom"].sum()
>>>>>>> pypsa-eur-sec/master

    return costs


<<<<<<< HEAD
=======
def calculate_cumulative_cost():
    planning_horizons = snakemake.config["scenario"]["planning_horizons"]

    cumulative_cost = pd.DataFrame(
        index=df["costs"].sum().index,
        columns=pd.Series(data=np.arange(0, 0.1, 0.01), name="social discount rate"),
    )

    # discount cost and express them in money value of planning_horizons[0]
    for r in cumulative_cost.columns:
        cumulative_cost[r] = [
            df["costs"].sum()[index] / ((1 + r) ** (index[-1] - planning_horizons[0]))
            for index in cumulative_cost.index
        ]

    # integrate cost throughout the transition path
    for r in cumulative_cost.columns:
        for cluster in cumulative_cost.index.get_level_values(level=0).unique():
            for lv in cumulative_cost.index.get_level_values(level=1).unique():
                for sector_opts in cumulative_cost.index.get_level_values(
                    level=2
                ).unique():
                    cumulative_cost.loc[
                        (cluster, lv, sector_opts, "cumulative cost"), r
                    ] = np.trapz(
                        cumulative_cost.loc[
                            idx[cluster, lv, sector_opts, planning_horizons], r
                        ].values,
                        x=planning_horizons,
                    )

    return cumulative_cost


def calculate_nodal_capacities(n, label, nodal_capacities):
    # Beware this also has extraneous locations for country (e.g. biomass) or continent-wide (e.g. fossil gas/oil) stuff
    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        nodal_capacities_c = c.df.groupby(["location", "carrier"])[
            opt_name.get(c.name, "p") + "_nom_opt"
        ].sum()
        index = pd.MultiIndex.from_tuples(
            [(c.list_name,) + t for t in nodal_capacities_c.index.to_list()]
        )
        nodal_capacities = nodal_capacities.reindex(index.union(nodal_capacities.index))
        nodal_capacities.loc[index, label] = nodal_capacities_c.values

    return nodal_capacities


def calculate_capacities(n, label, capacities):
    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        capacities_grouped = (
            c.df[opt_name.get(c.name, "p") + "_nom_opt"].groupby(c.df.carrier).sum()
        )
        capacities_grouped = pd.concat([capacities_grouped], keys=[c.list_name])

        capacities = capacities.reindex(
            capacities_grouped.index.union(capacities.index)
        )

        capacities.loc[capacities_grouped.index, label] = capacities_grouped

    return capacities


>>>>>>> pypsa-eur-sec/master
def calculate_curtailment(n, label, curtailment):
    avail = (
        n.generators_t.p_max_pu.multiply(n.generators.p_nom_opt)
        .sum()
        .groupby(n.generators.carrier)
        .sum()
    )
    used = n.generators_t.p.sum().groupby(n.generators.carrier).sum()

    curtailment[label] = (((avail - used) / avail) * 100).round(3)

    return curtailment


def calculate_energy(n, label, energy):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
<<<<<<< HEAD
        if c.name in {"Generator", "Load", "ShuntImpedance"}:
=======
        if c.name in n.one_port_components:
>>>>>>> pypsa-eur-sec/master
            c_energies = (
                c.pnl.p.multiply(n.snapshot_weightings.generators, axis=0)
                .sum()
                .multiply(c.df.sign)
                .groupby(c.df.carrier)
                .sum()
            )
<<<<<<< HEAD
        elif c.name in {"StorageUnit", "Store"}:
            c_energies = (
                c.pnl.p.multiply(n.snapshot_weightings.stores, axis=0)
                .sum()
                .multiply(c.df.sign)
                .groupby(c.df.carrier)
                .sum()
            )
        else:
            c_energies = (
                (
                    -c.pnl.p1.multiply(n.snapshot_weightings.generators, axis=0).sum()
                    - c.pnl.p0.multiply(n.snapshot_weightings.generators, axis=0).sum()
                )
                .groupby(c.df.carrier)
                .sum()
            )

        energy = include_in_summary(energy, [c.list_name], label, c_energies)
=======
        else:
            c_energies = pd.Series(0.0, c.df.carrier.unique())
            for port in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
                totals = (
                    c.pnl["p" + port]
                    .multiply(n.snapshot_weightings.generators, axis=0)
                    .sum()
                )
                # remove values where bus is missing (bug in nomopyomo)
                no_bus = c.df.index[c.df["bus" + port] == ""]
                totals.loc[no_bus] = n.component_attrs[c.name].loc[
                    "p" + port, "default"
                ]
                c_energies -= totals.groupby(c.df.carrier).sum()

        c_energies = pd.concat([c_energies], keys=[c.list_name])

        energy = energy.reindex(c_energies.index.union(energy.index))

        energy.loc[c_energies.index, label] = c_energies
>>>>>>> pypsa-eur-sec/master

    return energy


<<<<<<< HEAD
def include_in_summary(summary, multiindexprefix, label, item):
    # Index tuple(s) indicating the newly to-be-added row(s)
    raw_index = tuple([multiindexprefix, list(item.index)])
    summary = _add_indexed_rows(summary, raw_index)

    summary.loc[idx[raw_index], label] = item.values

    return summary


def calculate_capacity(n, label, capacity):
    for c in n.iterate_components(n.one_port_components):
        if "p_nom_opt" in c.df.columns:
            c_capacities = (
                abs(c.df.p_nom_opt.multiply(c.df.sign)).groupby(c.df.carrier).sum()
            )
            capacity = include_in_summary(capacity, [c.list_name], label, c_capacities)
        elif "e_nom_opt" in c.df.columns:
            c_capacities = (
                abs(c.df.e_nom_opt.multiply(c.df.sign)).groupby(c.df.carrier).sum()
            )
            capacity = include_in_summary(capacity, [c.list_name], label, c_capacities)

    for c in n.iterate_components(n.passive_branch_components):
        c_capacities = c.df["s_nom_opt"].groupby(c.df.carrier).sum()
        capacity = include_in_summary(capacity, [c.list_name], label, c_capacities)

    for c in n.iterate_components(n.controllable_branch_components):
        c_capacities = c.df.p_nom_opt.groupby(c.df.carrier).sum()
        capacity = include_in_summary(capacity, [c.list_name], label, c_capacities)

    return capacity


def calculate_supply(n, label, supply):
    """
    calculate the max dispatch of each component at the buses where the loads
    are attached.
    """
    load_types = n.buses.carrier.unique()

    for i in load_types:
        buses = n.buses.query("carrier == @i").index

        bus_map = pd.Series(False, index=n.buses.index)

        bus_map.loc[buses] = True

        for c in n.iterate_components(n.one_port_components):
            items = c.df.index[c.df.bus.map(bus_map)]

            if len(items) == 0 or c.pnl.p.empty:
=======
def calculate_supply(n, label, supply):
    """
    Calculate the max dispatch of each component at the buses aggregated by
    carrier.
    """

    bus_carriers = n.buses.carrier.unique()

    for i in bus_carriers:
        bus_map = n.buses.carrier == i
        bus_map.at[""] = False

        for c in n.iterate_components(n.one_port_components):
            items = c.df.index[c.df.bus.map(bus_map).fillna(False)]

            if len(items) == 0:
>>>>>>> pypsa-eur-sec/master
                continue

            s = (
                c.pnl.p[items]
                .max()
                .multiply(c.df.loc[items, "sign"])
                .groupby(c.df.loc[items, "carrier"])
                .sum()
            )
<<<<<<< HEAD

            # Index tuple(s) indicating the newly to-be-added row(s)
            raw_index = tuple([[i], [c.list_name], list(s.index)])
            supply = _add_indexed_rows(supply, raw_index)

            supply.loc[idx[raw_index], label] = s.values

        for c in n.iterate_components(n.branch_components):
            for end in ["0", "1"]:
                items = c.df.index[c.df["bus" + end].map(bus_map)]

                if len(items) == 0 or c.pnl["p" + end].empty:
=======
            s = pd.concat([s], keys=[c.list_name])
            s = pd.concat([s], keys=[i])

            supply = supply.reindex(s.index.union(supply.index))
            supply.loc[s.index, label] = s

        for c in n.iterate_components(n.branch_components):
            for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
                items = c.df.index[c.df["bus" + end].map(bus_map).fillna(False)]

                if len(items) == 0:
>>>>>>> pypsa-eur-sec/master
                    continue

                # lots of sign compensation for direction and to do maximums
                s = (-1) ** (1 - int(end)) * (
                    (-1) ** int(end) * c.pnl["p" + end][items]
                ).max().groupby(c.df.loc[items, "carrier"]).sum()
<<<<<<< HEAD

                supply = supply.reindex(
                    supply.index.union(
                        pd.MultiIndex.from_product([[i], [c.list_name], s.index])
                    )
                )
                supply.loc[idx[i, c.list_name, list(s.index)], label] = s.values
=======
                s.index = s.index + end
                s = pd.concat([s], keys=[c.list_name])
                s = pd.concat([s], keys=[i])

                supply = supply.reindex(s.index.union(supply.index))
                supply.loc[s.index, label] = s
>>>>>>> pypsa-eur-sec/master

    return supply


def calculate_supply_energy(n, label, supply_energy):
    """
<<<<<<< HEAD
    calculate the total dispatch of each component at the buses where the loads
    are attached.
    """
    load_types = n.buses.carrier.unique()

    for i in load_types:
        buses = n.buses.query("carrier == @i").index

        bus_map = pd.Series(False, index=n.buses.index)

        bus_map.loc[buses] = True

        for c in n.iterate_components(n.one_port_components):
            items = c.df.index[c.df.bus.map(bus_map)]

            if len(items) == 0 or c.pnl.p.empty:
=======
    Calculate the total energy supply/consuption of each component at the buses
    aggregated by carrier.
    """

    bus_carriers = n.buses.carrier.unique()

    for i in bus_carriers:
        bus_map = n.buses.carrier == i
        bus_map.at[""] = False

        for c in n.iterate_components(n.one_port_components):
            items = c.df.index[c.df.bus.map(bus_map).fillna(False)]

            if len(items) == 0:
>>>>>>> pypsa-eur-sec/master
                continue

            s = (
                c.pnl.p[items]
<<<<<<< HEAD
=======
                .multiply(n.snapshot_weightings.generators, axis=0)
>>>>>>> pypsa-eur-sec/master
                .sum()
                .multiply(c.df.loc[items, "sign"])
                .groupby(c.df.loc[items, "carrier"])
                .sum()
            )
<<<<<<< HEAD

            # Index tuple(s) indicating the newly to-be-added row(s)
            raw_index = tuple([[i], [c.list_name], list(s.index)])
            supply_energy = _add_indexed_rows(supply_energy, raw_index)

            supply_energy.loc[idx[raw_index], label] = s.values

        for c in n.iterate_components(n.branch_components):
            for end in ["0", "1"]:
                items = c.df.index[c.df["bus" + end].map(bus_map)]

                if len(items) == 0 or c.pnl["p" + end].empty:
                    continue

                s = (-1) * c.pnl["p" + end][items].sum().groupby(
                    c.df.loc[items, "carrier"]
                ).sum()

                supply_energy = supply_energy.reindex(
                    supply_energy.index.union(
                        pd.MultiIndex.from_product([[i], [c.list_name], s.index])
                    )
                )
                supply_energy.loc[idx[i, c.list_name, list(s.index)], label] = s.values
=======
            s = pd.concat([s], keys=[c.list_name])
            s = pd.concat([s], keys=[i])

            supply_energy = supply_energy.reindex(s.index.union(supply_energy.index))
            supply_energy.loc[s.index, label] = s

        for c in n.iterate_components(n.branch_components):
            for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
                items = c.df.index[c.df["bus" + str(end)].map(bus_map).fillna(False)]

                if len(items) == 0:
                    continue

                s = (-1) * c.pnl["p" + end][items].multiply(
                    n.snapshot_weightings.generators, axis=0
                ).sum().groupby(c.df.loc[items, "carrier"]).sum()
                s.index = s.index + end
                s = pd.concat([s], keys=[c.list_name])
                s = pd.concat([s], keys=[i])

                supply_energy = supply_energy.reindex(
                    s.index.union(supply_energy.index)
                )

                supply_energy.loc[s.index, label] = s
>>>>>>> pypsa-eur-sec/master

    return supply_energy


def calculate_metrics(n, label, metrics):
<<<<<<< HEAD
    metrics = metrics.reindex(
        metrics.index.union(
            pd.Index(
                [
                    "line_volume",
                    "line_volume_limit",
                    "line_volume_AC",
                    "line_volume_DC",
                    "line_volume_shadow",
                    "co2_shadow",
                ]
            )
        )
    )
=======
    metrics_list = [
        "line_volume",
        "line_volume_limit",
        "line_volume_AC",
        "line_volume_DC",
        "line_volume_shadow",
        "co2_shadow",
    ]

    metrics = metrics.reindex(pd.Index(metrics_list).union(metrics.index))
>>>>>>> pypsa-eur-sec/master

    metrics.at["line_volume_DC", label] = (n.links.length * n.links.p_nom_opt)[
        n.links.carrier == "DC"
    ].sum()
    metrics.at["line_volume_AC", label] = (n.lines.length * n.lines.s_nom_opt).sum()
    metrics.at["line_volume", label] = metrics.loc[
        ["line_volume_AC", "line_volume_DC"], label
    ].sum()

    if hasattr(n, "line_volume_limit"):
        metrics.at["line_volume_limit", label] = n.line_volume_limit
<<<<<<< HEAD

    if hasattr(n, "line_volume_limit_dual"):
=======
>>>>>>> pypsa-eur-sec/master
        metrics.at["line_volume_shadow", label] = n.line_volume_limit_dual

    if "CO2Limit" in n.global_constraints.index:
        metrics.at["co2_shadow", label] = n.global_constraints.at["CO2Limit", "mu"]

    return metrics


def calculate_prices(n, label, prices):
<<<<<<< HEAD
    bus_type = pd.Series(n.buses.index.str[3:], n.buses.index).replace(
        "", "electricity"
    )

    prices = prices.reindex(prices.index.union(bus_type.value_counts().index))

    logger.warning("Prices are time-averaged, not load-weighted")
    prices[label] = n.buses_t.marginal_price.mean().groupby(bus_type).mean()
=======
    prices = prices.reindex(prices.index.union(n.buses.carrier.unique()))

    # WARNING: this is time-averaged, see weighted_prices for load-weighted average
    prices[label] = n.buses_t.marginal_price.mean().groupby(n.buses.carrier).mean()
>>>>>>> pypsa-eur-sec/master

    return prices


def calculate_weighted_prices(n, label, weighted_prices):
<<<<<<< HEAD
    logger.warning("Weighted prices don't include storage units as loads")
=======
    # Warning: doesn't include storage units as loads
>>>>>>> pypsa-eur-sec/master

    weighted_prices = weighted_prices.reindex(
        pd.Index(
            [
                "electricity",
                "heat",
                "space heat",
                "urban heat",
                "space urban heat",
                "gas",
                "H2",
            ]
        )
    )

    link_loads = {
        "electricity": [
            "heat pump",
            "resistive heater",
            "battery charger",
            "H2 Electrolysis",
        ],
        "heat": ["water tanks charger"],
        "urban heat": ["water tanks charger"],
        "space heat": [],
        "space urban heat": [],
        "gas": ["OCGT", "gas boiler", "CHP electric", "CHP heat"],
        "H2": ["Sabatier", "H2 Fuel Cell"],
    }

    for carrier in link_loads:
        if carrier == "electricity":
            suffix = ""
        elif carrier[:5] == "space":
            suffix = carrier[5:]
        else:
            suffix = " " + carrier

        buses = n.buses.index[n.buses.index.str[2:] == suffix]

        if buses.empty:
            continue

        if carrier in ["H2", "gas"]:
            load = pd.DataFrame(index=n.snapshots, columns=buses, data=0.0)
        elif carrier[:5] == "space":
            load = heat_demand_df[buses.str[:2]].rename(
                columns=lambda i: str(i) + suffix
            )
        else:
            load = n.loads_t.p_set[buses]

        for tech in link_loads[carrier]:
            names = n.links.index[n.links.index.to_series().str[-len(tech) :] == tech]

            if names.empty:
                continue

            load += (
<<<<<<< HEAD
                n.links_t.p0[names]
                .groupby(n.links.loc[names, "bus0"], axis=1)
                .sum(axis=1)
            )

        # Add H2 Store when charging
        if carrier == "H2":
            stores = (
                n.stores_t.p[buses + " Store"]
                .groupby(n.stores.loc[buses + " Store", "bus"], axis=1)
                .sum(axis=1)
            )
            stores[stores > 0.0] = 0.0
            load += -stores
=======
                n.links_t.p0[names].groupby(n.links.loc[names, "bus0"], axis=1).sum()
            )

        # Add H2 Store when charging
        # if carrier == "H2":
        #    stores = n.stores_t.p[buses+ " Store"].groupby(n.stores.loc[buses+ " Store", "bus"],axis=1).sum(axis=1)
        #    stores[stores > 0.] = 0.
        #    load += -stores
>>>>>>> pypsa-eur-sec/master

        weighted_prices.loc[carrier, label] = (
            load * n.buses_t.marginal_price[buses]
        ).sum().sum() / load.sum().sum()

<<<<<<< HEAD
        if carrier[:5] == "space":
            print(load * n.buses_t.marginal_price[buses])
=======
        # still have no idea what this is for, only for debug reasons.
        if carrier[:5] == "space":
            logger.debug(load * n.buses_t.marginal_price[buses])
>>>>>>> pypsa-eur-sec/master

    return weighted_prices


<<<<<<< HEAD
outputs = [
    "costs",
    "curtailment",
    "energy",
    "capacity",
    "supply",
    "supply_energy",
    "prices",
    "weighted_prices",
    "metrics",
]


def make_summaries(networks_dict, paths, config, country="all"):
    columns = pd.MultiIndex.from_tuples(
        networks_dict.keys(), names=["simpl", "clusters", "ll", "opts"]
    )

    dfs = {}

    for output in outputs:
        dfs[output] = pd.DataFrame(columns=columns, dtype=float)

    for label, filename in networks_dict.items():
        print(label, filename)
        if not os.path.exists(filename):
            print("does not exist!!")
            continue

        try:
            n = pypsa.Network(filename)
        except OSError:
            logger.warning("Skipping {filename}".format(filename=filename))
            continue

        if country != "all":
            n = n[n.buses.country == country]

        Nyears = n.snapshot_weightings.objective.sum() / 8760.0
        costs = load_costs(paths[0], config["costs"], config["electricity"], Nyears)
        update_transmission_costs(n, costs)

        assign_carriers(n)

        for output in outputs:
            dfs[output] = globals()["calculate_" + output](n, label, dfs[output])

    return dfs


def to_csv(dfs, dir):
    os.makedirs(dir, exist_ok=True)
    for key, df in dfs.items():
        df.to_csv(os.path.join(dir, f"{key}.csv"))
=======
def calculate_market_values(n, label, market_values):
    # Warning: doesn't include storage units

    carrier = "AC"

    buses = n.buses.index[n.buses.carrier == carrier]

    ## First do market value of generators ##

    generators = n.generators.index[n.buses.loc[n.generators.bus, "carrier"] == carrier]

    techs = n.generators.loc[generators, "carrier"].value_counts().index

    market_values = market_values.reindex(market_values.index.union(techs))

    for tech in techs:
        gens = generators[n.generators.loc[generators, "carrier"] == tech]

        dispatch = (
            n.generators_t.p[gens]
            .groupby(n.generators.loc[gens, "bus"], axis=1)
            .sum()
            .reindex(columns=buses, fill_value=0.0)
        )

        revenue = dispatch * n.buses_t.marginal_price[buses]

        market_values.at[tech, label] = revenue.sum().sum() / dispatch.sum().sum()

    ## Now do market value of links ##

    for i in ["0", "1"]:
        all_links = n.links.index[n.buses.loc[n.links["bus" + i], "carrier"] == carrier]

        techs = n.links.loc[all_links, "carrier"].value_counts().index

        market_values = market_values.reindex(market_values.index.union(techs))

        for tech in techs:
            links = all_links[n.links.loc[all_links, "carrier"] == tech]

            dispatch = (
                n.links_t["p" + i][links]
                .groupby(n.links.loc[links, "bus" + i], axis=1)
                .sum()
                .reindex(columns=buses, fill_value=0.0)
            )

            revenue = dispatch * n.buses_t.marginal_price[buses]

            market_values.at[tech, label] = revenue.sum().sum() / dispatch.sum().sum()

    return market_values


def calculate_price_statistics(n, label, price_statistics):
    price_statistics = price_statistics.reindex(
        price_statistics.index.union(
            pd.Index(["zero_hours", "mean", "standard_deviation"])
        )
    )

    buses = n.buses.index[n.buses.carrier == "AC"]

    threshold = 0.1  # higher than phoney marginal_cost of wind/solar

    df = pd.DataFrame(data=0.0, columns=buses, index=n.snapshots)

    df[n.buses_t.marginal_price[buses] < threshold] = 1.0

    price_statistics.at["zero_hours", label] = df.sum().sum() / (
        df.shape[0] * df.shape[1]
    )

    price_statistics.at["mean", label] = (
        n.buses_t.marginal_price[buses].unstack().mean()
    )

    price_statistics.at["standard_deviation", label] = (
        n.buses_t.marginal_price[buses].unstack().std()
    )

    return price_statistics


def make_summaries(networks_dict):
    outputs = [
        "nodal_costs",
        "nodal_capacities",
        "nodal_cfs",
        "cfs",
        "costs",
        "capacities",
        "curtailment",
        "energy",
        "supply",
        "supply_energy",
        "prices",
        "weighted_prices",
        "price_statistics",
        "market_values",
        "metrics",
    ]

    columns = pd.MultiIndex.from_tuples(
        networks_dict.keys(), names=["cluster", "lv", "opt", "planning_horizon"]
    )

    df = {}

    for output in outputs:
        df[output] = pd.DataFrame(columns=columns, dtype=float)

    for label, filename in networks_dict.items():
        logger.info(f"Make summary for scenario {label}, using {filename}")

        overrides = override_component_attrs(snakemake.input.overrides)
        n = pypsa.Network(filename, override_component_attrs=overrides)

        assign_carriers(n)
        assign_locations(n)

        for output in outputs:
            df[output] = globals()["calculate_" + output](n, label, df[output])

    return df


def to_csv(df):
    for key in df:
        df[key].to_csv(snakemake.output[key])
>>>>>>> pypsa-eur-sec/master


if __name__ == "__main__":
    if "snakemake" not in globals():
<<<<<<< HEAD
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "make_summary",
            simpl="",
            clusters="5",
            ll="copt",
            opts="Co2L-24H",
            country="all",
        )
        network_dir = os.path.join(
            "..", "results", "networks", snakemake.config["run"]["name"]
        )
    else:
        network_dir = os.path.join(
            "results", "networks", snakemake.config["run"]["name"]
        )
    configure_logging(snakemake)

    config = snakemake.config
    wildcards = snakemake.wildcards

    def expand_from_wildcard(key, config):
        w = getattr(wildcards, key)
        return config["scenario"][key] if w == "all" else [w]

    if wildcards.ll.endswith("all"):
        ll = config["scenario"]["ll"]
        if len(wildcards.ll) == 4:
            ll = [l for l in ll if l[0] == wildcards.ll[0]]
    else:
        ll = [wildcards.ll]

    networks_dict = {
        (simpl, clusters, l, opts): os.path.join(
            network_dir, f"elec_s{simpl}_" f"{clusters}_ec_l{l}_{opts}.nc"
        )
        for simpl in expand_from_wildcard("simpl", config)
        for clusters in expand_from_wildcard("clusters", config)
        for l in ll
        for opts in expand_from_wildcard("opts", config)
    }

    dfs = make_summaries(
        networks_dict, snakemake.input, config, country=wildcards.country
    )

    to_csv(dfs, snakemake.output[0])
=======
        from helper import mock_snakemake

        snakemake = mock_snakemake("make_summary")

    logging.basicConfig(level=snakemake.config["logging_level"])

    networks_dict = {
        (cluster, lv, opt + sector_opt, planning_horizon): snakemake.config[
            "results_dir"
        ]
        + snakemake.config["run"]
        + f"/postnetworks/elec_s{simpl}_{cluster}_lv{lv}_{opt}_{sector_opt}_{planning_horizon}.nc"
        for simpl in snakemake.config["scenario"]["simpl"]
        for cluster in snakemake.config["scenario"]["clusters"]
        for opt in snakemake.config["scenario"]["opts"]
        for sector_opt in snakemake.config["scenario"]["sector_opts"]
        for lv in snakemake.config["scenario"]["lv"]
        for planning_horizon in snakemake.config["scenario"]["planning_horizons"]
    }

    Nyears = 1

    costs_db = prepare_costs(
        snakemake.input.costs,
        snakemake.config["costs"]["USD2013_to_EUR2013"],
        snakemake.config["costs"]["discountrate"],
        Nyears,
        snakemake.config["costs"]["lifetime"],
    )

    df = make_summaries(networks_dict)

    df["metrics"].loc["total costs"] = df["costs"].sum()

    to_csv(df)

    if snakemake.config["foresight"] == "myopic":
        cumulative_cost = calculate_cumulative_cost()
        cumulative_cost.to_csv(
            snakemake.config["summary_dir"]
            + "/"
            + snakemake.config["run"]
            + "/csvs/cumulative_cost.csv"
        )
>>>>>>> pypsa-eur-sec/master
