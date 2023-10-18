# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Create summary CSV files for all scenario runs with perfect foresight including
costs, capacities, capacity factors, curtailment, energy balances, prices and
other metrics.
"""


import numpy as np
import pandas as pd
import pypsa
from make_summary import (
    assign_carriers,
    assign_locations,
    calculate_cfs,
    calculate_nodal_cfs,
    calculate_nodal_costs,
)
from prepare_sector_network import prepare_costs
from pypsa.descriptors import get_active_assets, nominal_attrs
from six import iteritems

idx = pd.IndexSlice

opt_name = {"Store": "e", "Line": "s", "Transformer": "s"}


def calculate_costs(n, label, costs):
    investments = n.investment_periods
    cols = pd.MultiIndex.from_product(
        [
            costs.columns.levels[0],
            costs.columns.levels[1],
            costs.columns.levels[2],
            investments,
        ],
        names=costs.columns.names[:3] + ["year"],
    )
    costs = costs.reindex(cols, axis=1)

    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        capital_costs = c.df.capital_cost * c.df[opt_name.get(c.name, "p") + "_nom_opt"]
        active = pd.concat(
            [
                get_active_assets(n, c.name, inv_p).rename(inv_p)
                for inv_p in investments
            ],
            axis=1,
        ).astype(int)
        capital_costs = active.mul(capital_costs, axis=0)
        discount = (
            n.investment_period_weightings["objective"]
            / n.investment_period_weightings["years"]
        )
        capital_costs_grouped = capital_costs.groupby(c.df.carrier).sum().mul(discount)

        capital_costs_grouped = pd.concat([capital_costs_grouped], keys=["capital"])
        capital_costs_grouped = pd.concat([capital_costs_grouped], keys=[c.list_name])

        costs = costs.reindex(capital_costs_grouped.index.union(costs.index))

        costs.loc[capital_costs_grouped.index, label] = capital_costs_grouped.values

        if c.name == "Link":
            p = (
                c.pnl.p0.multiply(n.snapshot_weightings.generators, axis=0)
                .groupby(level=0)
                .sum()
            )
        elif c.name == "Line":
            continue
        elif c.name == "StorageUnit":
            p_all = c.pnl.p.multiply(n.snapshot_weightings.stores, axis=0)
            p_all[p_all < 0.0] = 0.0
            p = p_all.groupby(level=0).sum()
        else:
            p = (
                round(c.pnl.p, ndigits=2)
                .multiply(n.snapshot_weightings.generators, axis=0)
                .groupby(level=0)
                .sum()
            )

        # correct sequestration cost
        if c.name == "Store":
            items = c.df.index[
                (c.df.carrier == "co2 stored") & (c.df.marginal_cost <= -100.0)
            ]
            c.df.loc[items, "marginal_cost"] = -20.0

        marginal_costs = p.mul(c.df.marginal_cost).T
        # marginal_costs = active.mul(marginal_costs, axis=0)
        marginal_costs_grouped = (
            marginal_costs.groupby(c.df.carrier).sum().mul(discount)
        )

        marginal_costs_grouped = pd.concat([marginal_costs_grouped], keys=["marginal"])
        marginal_costs_grouped = pd.concat([marginal_costs_grouped], keys=[c.list_name])

        costs = costs.reindex(marginal_costs_grouped.index.union(costs.index))

        costs.loc[marginal_costs_grouped.index, label] = marginal_costs_grouped.values

    # add back in all hydro
    # costs.loc[("storage_units","capital","hydro"),label] = (0.01)*2e6*n.storage_units.loc[n.storage_units.group=="hydro","p_nom"].sum()
    # costs.loc[("storage_units","capital","PHS"),label] = (0.01)*2e6*n.storage_units.loc[n.storage_units.group=="PHS","p_nom"].sum()
    # costs.loc[("generators","capital","ror"),label] = (0.02)*3e6*n.generators.loc[n.generators.group=="ror","p_nom"].sum()

    return costs


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
    investments = n.investment_periods
    cols = pd.MultiIndex.from_product(
        [
            capacities.columns.levels[0],
            capacities.columns.levels[1],
            capacities.columns.levels[2],
            investments,
        ],
        names=capacities.columns.names[:3] + ["year"],
    )
    capacities = capacities.reindex(cols, axis=1)

    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        active = pd.concat(
            [
                get_active_assets(n, c.name, inv_p).rename(inv_p)
                for inv_p in investments
            ],
            axis=1,
        ).astype(int)
        caps = c.df[opt_name.get(c.name, "p") + "_nom_opt"]
        caps = active.mul(caps, axis=0)
        capacities_grouped = (
            caps.groupby(c.df.carrier).sum().drop("load", errors="ignore")
        )
        capacities_grouped = pd.concat([capacities_grouped], keys=[c.list_name])

        capacities = capacities.reindex(
            capacities_grouped.index.union(capacities.index)
        )

        capacities.loc[capacities_grouped.index, label] = capacities_grouped.values

    return capacities


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
    investments = n.investment_periods
    cols = pd.MultiIndex.from_product(
        [
            energy.columns.levels[0],
            energy.columns.levels[1],
            energy.columns.levels[2],
            investments,
        ],
        names=energy.columns.names[:3] + ["year"],
    )
    energy = energy.reindex(cols, axis=1)

    for c in n.iterate_components(n.one_port_components | n.branch_components):
        if c.name in n.one_port_components:
            c_energies = (
                c.pnl.p.multiply(n.snapshot_weightings.generators, axis=0)
                .groupby(level=0)
                .sum()
                .multiply(c.df.sign)
                .groupby(c.df.carrier, axis=1)
                .sum()
            )
        else:
            c_energies = pd.DataFrame(
                0.0, columns=c.df.carrier.unique(), index=n.investment_periods
            )
            for port in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
                totals = (
                    c.pnl["p" + port]
                    .multiply(n.snapshot_weightings.generators, axis=0)
                    .groupby(level=0)
                    .sum()
                )
                # remove values where bus is missing (bug in nomopyomo)
                no_bus = c.df.index[c.df["bus" + port] == ""]
                totals[no_bus] = float(
                    n.component_attrs[c.name].loc["p" + port, "default"]
                )
                c_energies -= totals.groupby(c.df.carrier, axis=1).sum()

        c_energies = pd.concat([c_energies.T], keys=[c.list_name])

        energy = energy.reindex(c_energies.index.union(energy.index))

        energy.loc[c_energies.index, label] = c_energies.values

    return energy


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
                continue

            s = (
                c.pnl.p[items]
                .max()
                .multiply(c.df.loc[items, "sign"])
                .groupby(c.df.loc[items, "carrier"])
                .sum()
            )
            s = pd.concat([s], keys=[c.list_name])
            s = pd.concat([s], keys=[i])

            supply = supply.reindex(s.index.union(supply.index))
            supply.loc[s.index, label] = s

        for c in n.iterate_components(n.branch_components):
            for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
                items = c.df.index[c.df["bus" + end].map(bus_map).fillna(False)]

                if len(items) == 0:
                    continue

                # lots of sign compensation for direction and to do maximums
                s = (-1) ** (1 - int(end)) * (
                    (-1) ** int(end) * c.pnl["p" + end][items]
                ).max().groupby(c.df.loc[items, "carrier"]).sum()
                s.index = s.index + end
                s = pd.concat([s], keys=[c.list_name])
                s = pd.concat([s], keys=[i])

                supply = supply.reindex(s.index.union(supply.index))
                supply.loc[s.index, label] = s

    return supply


def calculate_supply_energy(n, label, supply_energy):
    """
    Calculate the total energy supply/consuption of each component at the buses
    aggregated by carrier.
    """

    investments = n.investment_periods
    cols = pd.MultiIndex.from_product(
        [
            supply_energy.columns.levels[0],
            supply_energy.columns.levels[1],
            supply_energy.columns.levels[2],
            investments,
        ],
        names=supply_energy.columns.names[:3] + ["year"],
    )
    supply_energy = supply_energy.reindex(cols, axis=1)

    bus_carriers = n.buses.carrier.unique()

    for i in bus_carriers:
        bus_map = n.buses.carrier == i
        bus_map.at[""] = False

        for c in n.iterate_components(n.one_port_components):
            items = c.df.index[c.df.bus.map(bus_map).fillna(False)]

            if len(items) == 0:
                continue

            if c.name == "Generator":
                weightings = n.snapshot_weightings.generators
            else:
                weightings = n.snapshot_weightings.stores

            if i in ["oil", "co2", "H2"]:
                if c.name == "Load":
                    c.df.loc[items, "carrier"] = [
                        load.split("-202")[0] for load in items
                    ]
                if i == "oil" and c.name == "Generator":
                    c.df.loc[items, "carrier"] = "imported oil"
            s = (
                c.pnl.p[items]
                .multiply(weightings, axis=0)
                .groupby(level=0)
                .sum()
                .multiply(c.df.loc[items, "sign"])
                .groupby(c.df.loc[items, "carrier"], axis=1)
                .sum()
                .T
            )
            s = pd.concat([s], keys=[c.list_name])
            s = pd.concat([s], keys=[i])

            supply_energy = supply_energy.reindex(
                s.index.union(supply_energy.index, sort=False)
            )
            supply_energy.loc[s.index, label] = s.values

        for c in n.iterate_components(n.branch_components):
            for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
                items = c.df.index[c.df[f"bus{str(end)}"].map(bus_map).fillna(False)]

                if len(items) == 0:
                    continue

                s = (
                    (-1)
                    * c.pnl["p" + end]
                    .reindex(items, axis=1)
                    .multiply(n.snapshot_weightings.objective, axis=0)
                    .groupby(level=0)
                    .sum()
                    .groupby(c.df.loc[items, "carrier"], axis=1)
                    .sum()
                ).T
                s.index = s.index + end
                s = pd.concat([s], keys=[c.list_name])
                s = pd.concat([s], keys=[i])

                supply_energy = supply_energy.reindex(
                    s.index.union(supply_energy.index, sort=False)
                )

                supply_energy.loc[s.index, label] = s.values

    return supply_energy


def calculate_metrics(n, label, metrics):
    metrics = metrics.reindex(
        pd.Index(
            [
                "line_volume",
                "line_volume_limit",
                "line_volume_AC",
                "line_volume_DC",
                "line_volume_shadow",
                "co2_shadow",
            ]
        ).union(metrics.index)
    )

    metrics.at["line_volume_DC", label] = (n.links.length * n.links.p_nom_opt)[
        n.links.carrier == "DC"
    ].sum()
    metrics.at["line_volume_AC", label] = (n.lines.length * n.lines.s_nom_opt).sum()
    metrics.at["line_volume", label] = metrics.loc[
        ["line_volume_AC", "line_volume_DC"], label
    ].sum()

    if hasattr(n, "line_volume_limit"):
        metrics.at["line_volume_limit", label] = n.line_volume_limit
        metrics.at["line_volume_shadow", label] = n.line_volume_limit_dual

    if "CO2Limit" in n.global_constraints.index:
        metrics.at["co2_shadow", label] = n.global_constraints.at["CO2Limit", "mu"]

    return metrics


def calculate_prices(n, label, prices):
    prices = prices.reindex(prices.index.union(n.buses.carrier.unique()))

    # WARNING: this is time-averaged, see weighted_prices for load-weighted average
    prices[label] = n.buses_t.marginal_price.mean().groupby(n.buses.carrier).mean()

    return prices


def calculate_weighted_prices(n, label, weighted_prices):
    # Warning: doesn't include storage units as loads

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

    for carrier, value in link_loads.items():
        if carrier == "electricity":
            suffix = ""
        elif carrier[:5] == "space":
            suffix = carrier[5:]
        else:
            suffix = " " + carrier

        buses = n.buses.index[n.buses.index.str[2:] == suffix]

        if buses.empty:
            continue

        load = (
            pd.DataFrame(index=n.snapshots, columns=buses, data=0.0)
            if carrier in ["H2", "gas"]
            else n.loads_t.p_set.reindex(buses, axis=1)
        )
        for tech in value:
            names = n.links.index[n.links.index.to_series().str[-len(tech) :] == tech]

            if names.empty:
                continue

            load += (
                n.links_t.p0[names].groupby(n.links.loc[names, "bus0"], axis=1).sum()
            )

        # Add H2 Store when charging
        # if carrier == "H2":
        #    stores = n.stores_t.p[buses+ " Store"].groupby(n.stores.loc[buses+ " Store","bus"],axis=1).sum(axis=1)
        #    stores[stores > 0.] = 0.
        #    load += -stores

        weighted_prices.loc[carrier, label] = (
            load * n.buses_t.marginal_price[buses]
        ).sum().sum() / load.sum().sum()

        if carrier[:5] == "space":
            print(load * n.buses_t.marginal_price[buses])

    return weighted_prices


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

    price_statistics.at["mean", label] = n.buses_t.marginal_price[buses].mean().mean()

    price_statistics.at["standard_deviation", label] = (
        n.buses_t.marginal_price[buses].droplevel(0).unstack().std()
    )

    return price_statistics


def calculate_co2_emissions(n, label, df):
    carattr = "co2_emissions"
    emissions = n.carriers.query(f"{carattr} != 0")[carattr]

    if emissions.empty:
        return

    weightings = n.snapshot_weightings.generators.mul(
        n.investment_period_weightings["years"]
        .reindex(n.snapshots)
        .fillna(method="bfill")
        .fillna(1.0),
        axis=0,
    )

    # generators
    gens = n.generators.query("carrier in @emissions.index")
    if not gens.empty:
        em_pu = gens.carrier.map(emissions) / gens.efficiency
        em_pu = (
            weightings["generators"].to_frame("weightings")
            @ em_pu.to_frame("weightings").T
        )
        emitted = n.generators_t.p[gens.index].mul(em_pu)

        emitted_grouped = (
            emitted.groupby(level=0).sum().groupby(n.generators.carrier, axis=1).sum().T
        )

        df = df.reindex(emitted_grouped.index.union(df.index))

        df.loc[emitted_grouped.index, label] = emitted_grouped.values

    if any(n.stores.carrier == "co2"):
        co2_i = n.stores[n.stores.carrier == "co2"].index
        df[label] = n.stores_t.e.groupby(level=0).last()[co2_i].iloc[:, 0]

    return df


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
    "co2_emissions",
]


def make_summaries(networks_dict):
    columns = pd.MultiIndex.from_tuples(
        networks_dict.keys(), names=["cluster", "lv", "opt"]
    )
    df = {}

    for output in outputs:
        df[output] = pd.DataFrame(columns=columns, dtype=float)

    for label, filename in iteritems(networks_dict):
        print(label, filename)
        try:
            n = pypsa.Network(filename)
        except OSError:
            print(label, " not solved yet.")
            continue
            # del networks_dict[label]

        if not hasattr(n, "objective"):
            print(label, " not solved correctly. Check log if infeasible or unbounded.")
            continue
        assign_carriers(n)
        assign_locations(n)

        for output in outputs:
            df[output] = globals()["calculate_" + output](n, label, df[output])

    return df


def to_csv(df):
    for key in df:
        df[key] = df[key].apply(lambda x: pd.to_numeric(x))
        df[key].to_csv(snakemake.output[key])


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("make_summary_perfect")

    run = snakemake.config["run"]["name"]
    if run != "":
        run += "/"

    networks_dict = {
        (clusters, lv, opts + sector_opts): "results/"
        + run
        + f"postnetworks/elec_s{simpl}_{clusters}_l{lv}_{opts}_{sector_opts}_brownfield_all_years.nc"
        for simpl in snakemake.config["scenario"]["simpl"]
        for clusters in snakemake.config["scenario"]["clusters"]
        for opts in snakemake.config["scenario"]["opts"]
        for sector_opts in snakemake.config["scenario"]["sector_opts"]
        for lv in snakemake.config["scenario"]["ll"]
    }

    print(networks_dict)

    nyears = 1
    costs_db = prepare_costs(
        snakemake.input.costs,
        snakemake.config["costs"],
        nyears,
    )

    df = make_summaries(networks_dict)

    df["metrics"].loc["total costs"] = df["costs"].sum().groupby(level=[0, 1, 2]).sum()

    to_csv(df)
