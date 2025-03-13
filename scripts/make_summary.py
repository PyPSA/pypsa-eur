# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create summary CSV files for all scenario runs including costs, capacities,
capacity factors, curtailment, energy balances, prices and other metrics.
"""

import logging
import sys

import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging, set_scenario_config

idx = pd.IndexSlice
logger = logging.getLogger(__name__)
OPT_NAME = {"Store": "e", "Line": "s", "Transformer": "s"}

OUTPUTS = [
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
    "nodal_supply_energy",
    "prices",
    "weighted_prices",
    "price_statistics",
    "market_values",
    "metrics",
]


def assign_carriers(n):
    if "carrier" not in n.lines:
        n.lines["carrier"] = "AC"


def assign_locations(n):
    for c in n.iterate_components(n.one_port_components):
        c.df["location"] = c.df.bus.map(n.buses.location)

    for c in n.iterate_components(n.branch_components):
        c_bus_cols = c.df.filter(regex="^bus")
        locs = c_bus_cols.apply(lambda c: c.map(n.buses.location))
        # take the longest location string for each row;
        # links and lines inherit location of highest resolved bus;
        # regional locations "BE0 0" are longer than "EU" by design
        c.df["location"] = locs.apply(lambda row: max(row.dropna(), key=len), axis=1)

def calculate_nodal_cfs(n):
    nodal_cfs = []

    # Beware this also has extraneous locations for country (e.g. biomass) or continent-wide (e.g. fossil gas/oil) stuff
    for c in n.iterate_components(
        (n.branch_components ^ {"Line", "Transformer"})
        | n.controllable_one_port_components ^ {"Load", "StorageUnit"}
    ):
        capacities_c = c.df.groupby(["location", "carrier"])[
            OPT_NAME.get(c.name, "p") + "_nom_opt"
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
        cf_c = pd.Series(cf_c.values, index=index)
        nodal_cfs.append(cf_c)

    return pd.concat(nodal_cfs).sort_index()


def calculate_cfs(n):
    cfs = []

    for c in n.iterate_components(
        n.branch_components
        | n.controllable_one_port_components ^ {"Load", "StorageUnit"}
    ):
        capacities_c = (
            c.df[OPT_NAME.get(c.name, "p") + "_nom_opt"].groupby(c.df.carrier).sum()
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

        cfs.append(cf_c)

    return pd.concat(cfs).sort_index()


def calculate_nodal_costs(n):
    nodal_costs = []

    # Beware this also has extraneous locations for country (e.g. biomass) or continent-wide (e.g. fossil gas/oil) stuff
    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        c.df["capital_costs"] = (
            c.df.capital_cost * c.df[OPT_NAME.get(c.name, "p") + "_nom_opt"]
        )
        capital_costs = c.df.groupby(["location", "carrier"])["capital_costs"].sum()
        index = pd.MultiIndex.from_tuples(
            [(c.list_name, "capital") + t for t in capital_costs.index.to_list()]
        )
        capital_costs = pd.Series(capital_costs.values, index=index)
        nodal_costs.append(capital_costs)

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

        c.df["marginal_costs"] = p * c.df.marginal_cost
        marginal_costs = c.df.groupby(["location", "carrier"])["marginal_costs"].sum()
        index = pd.MultiIndex.from_tuples(
            [(c.list_name, "marginal") + t for t in marginal_costs.index.to_list()]
        )
        marginal_costs = pd.Series(marginal_costs.values, index=index)
        nodal_costs.append(marginal_costs)

    return pd.concat(nodal_costs).sort_index()


def calculate_costs(n):
    costs = []
    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        capital_costs = c.df.capital_cost * c.df[OPT_NAME.get(c.name, "p") + "_nom_opt"]
        capital_costs_grouped = capital_costs.groupby(c.df.carrier).sum()

        capital_costs_grouped = pd.concat([capital_costs_grouped], keys=["capital"])
        capital_costs_grouped = pd.concat([capital_costs_grouped], keys=[c.list_name])

        costs.append(capital_costs_grouped)

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

        marginal_costs = p * c.df.marginal_cost

        marginal_costs_grouped = marginal_costs.groupby(c.df.carrier).sum()

        marginal_costs_grouped = pd.concat([marginal_costs_grouped], keys=["marginal"])
        marginal_costs_grouped = pd.concat([marginal_costs_grouped], keys=[c.list_name])

        costs.append(marginal_costs_grouped)

    return pd.concat(costs).sort_index()


def calculate_nodal_capacities(n):
    # Beware this also has extraneous locations for country (e.g. biomass) or continent-wide (e.g. fossil gas/oil) stuff
    nodal_capacities = []
    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        nodal_capacities_c = c.df.groupby(["location", "carrier"])[
            OPT_NAME.get(c.name, "p") + "_nom_opt"
        ].sum()
        index = pd.MultiIndex.from_tuples(
            [(c.list_name,) + t for t in nodal_capacities_c.index.to_list()]
        )
        nodal_capacities_c = pd.Series(nodal_capacities_c.values, index=index)
        nodal_capacities.append(nodal_capacities_c)

    return pd.concat(nodal_capacities).sort_index()


def calculate_capacities(n):
    capacities = []

    for c in n.iterate_components(
        n.branch_components | n.controllable_one_port_components ^ {"Load"}
    ):
        capacities_grouped = (
            c.df[OPT_NAME.get(c.name, "p") + "_nom_opt"].groupby(c.df.carrier).sum()
        )
        capacities_grouped = pd.concat([capacities_grouped], keys=[c.list_name])

        capacities.append(capacities_grouped)

    return pd.concat(capacities).sort_index()


def calculate_curtailment(n):
    avail = (
        n.generators_t.p_max_pu.multiply(n.generators.p_nom_opt)
        .sum()
        .groupby(n.generators.carrier)
        .sum()
    )
    used = n.generators_t.p.sum().groupby(n.generators.carrier).sum()

    curtailment = (((avail - used) / avail) * 100).round(3)

    return curtailment.sort_index()


def calculate_energy(n):
    energy = []
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        if c.name in n.one_port_components:
            c_energies = (
                c.pnl.p.multiply(n.snapshot_weightings.generators, axis=0)
                .sum()
                .multiply(c.df.sign)
                .groupby(c.df.carrier)
                .sum()
            )
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
                totals.loc[no_bus] = float(
                    n.component_attrs[c.name].loc["p" + port, "default"]
                )
                c_energies -= totals.groupby(c.df.carrier).sum()

        c_energies = pd.concat([c_energies], keys=[c.list_name])

        energy.append(c_energies)

    return pd.concat(energy)


def calculate_supply(n):
    """
    Calculate the max dispatch of each component at the buses aggregated by
    carrier.
    """
    bus_carriers = n.buses.carrier.unique()

    supply = []

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

            supply.append(s)

        for c in n.iterate_components(n.branch_components):
            for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
                items = c.df.index[c.df["bus" + end].map(bus_map).fillna(False)]

                if len(items) == 0 or c.pnl["p" + end].empty:
                    continue

                # lots of sign compensation for direction and to do maximums
                s = (-1) ** (1 - int(end)) * (
                    (-1) ** int(end) * c.pnl["p" + end][items]
                ).max().groupby(c.df.loc[items, "carrier"]).sum()
                s.index = s.index + end
                s = pd.concat([s], keys=[c.list_name])
                s = pd.concat([s], keys=[i])

                supply.append(s)

    return pd.concat(supply).sort_index()


def calculate_supply_energy(n):
    """
    Calculate the total energy supply/consumption of each component at the buses
    aggregated by carrier.
    """
    bus_carriers = n.buses.carrier.unique()

    supply_energy = []

    for i in bus_carriers:
        bus_map = n.buses.carrier == i
        bus_map.at[""] = False

        for c in n.iterate_components(n.one_port_components):
            items = c.df.index[c.df.bus.map(bus_map).fillna(False)]

            if len(items) == 0:
                continue

            s = (
                c.pnl.p[items]
                .multiply(n.snapshot_weightings.generators, axis=0)
                .sum()
                .multiply(c.df.loc[items, "sign"])
                .groupby(c.df.loc[items, "carrier"])
                .sum()
            )
            s = pd.concat([s], keys=[c.list_name])
            s = pd.concat([s], keys=[i])

            supply_energy.append(s)

        for c in n.iterate_components(n.branch_components):
            for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
                items = c.df.index[c.df[f"bus{str(end)}"].map(bus_map).fillna(False)]

                if len(items) == 0 or c.pnl["p" + end].empty:
                    continue

                s = (-1) * c.pnl["p" + end][items].multiply(
                    n.snapshot_weightings.generators, axis=0
                ).sum().groupby(c.df.loc[items, "carrier"]).sum()
                s.index = s.index + end
                s = pd.concat([s], keys=[c.list_name])
                s = pd.concat([s], keys=[i])

                supply_energy.append(s)

    return pd.concat(supply_energy)


def calculate_nodal_supply_energy(n):
    """
    Calculate the total energy supply/consumption of each component at the
    buses aggregated by carrier and node.
    """

    bus_carriers = n.buses.carrier.unique()

    nodal_supply_energy = []

    for i in bus_carriers:
        bus_map = n.buses.carrier == i
        bus_map.at[""] = False

        for c in n.iterate_components(n.one_port_components):
            items = c.df.index[c.df.bus.map(bus_map).fillna(False)]

            if len(items) == 0:
                continue

            s = (
                pd.concat(
                    [
                        (
                            c.pnl.p[items]
                            .multiply(n.snapshot_weightings.generators, axis=0)
                            .sum()
                            .multiply(c.df.loc[items, "sign"])
                        ),
                        c.df.loc[items][["bus", "carrier"]],
                    ],
                    axis=1,
                )
                .groupby(by=["bus", "carrier"])
                .sum()[0]
            )
            s = pd.concat([s], keys=[c.list_name])
            s = pd.concat([s], keys=[i])

            nodal_supply_energy.append(s)

        for c in n.iterate_components(n.branch_components):
            for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
                items = c.df.index[c.df["bus" + str(end)].map(bus_map).fillna(False)]

                if (len(items) == 0) or c.pnl["p" + end].empty:
                    continue

                s = (
                    pd.concat(
                        [
                            (
                                (-1)
                                * c.pnl["p" + end][items]
                                .multiply(n.snapshot_weightings.generators, axis=0)
                                .sum()
                            ),
                            c.df.loc[items][["bus0", "carrier"]],
                        ],
                        axis=1,
                    )
                    .groupby(by=["bus0", "carrier"])
                    .sum()[0]
                )

                s.index = s.index.map(lambda x: (x[0], x[1] + end))
                s = pd.concat([s], keys=[c.list_name])
                s = pd.concat([s], keys=[i])

                nodal_supply_energy.append(s)

    return pd.concat(nodal_supply_energy).sort_index()


def calculate_metrics(n):
    metrics = {}

    metrics["line_volume_DC"] = (n.links.length * n.links.p_nom_opt)[
        n.links.carrier == "DC"
    ].sum()
    metrics["line_volume_AC"] = (n.lines.length * n.lines.s_nom_opt).sum()
    metrics["line_volume"] = metrics["line_volume_AC"] + metrics["line_volume_DC"]
    metrics["total costs"] = n.statistics.capex().sum() + n.statistics.opex().sum()

    if "lv_limit" in n.global_constraints.index:
        metrics["line_volume_limit"] = n.global_constraints.at["lv_limit", "constant"]
        metrics["line_volume_shadow"] = n.global_constraints.at["lv_limit", "mu"]

    if "CO2Limit" in n.global_constraints.index:
        metrics["co2_shadow"] = n.global_constraints.at["CO2Limit", "mu"]

    if "co2_sequestration_limit" in n.global_constraints.index:
        metrics["co2_storage_shadow"] = n.global_constraints.at[
            "co2_sequestration_limit", "mu"
        ]

    return pd.Series(metrics).sort_index()


def calculate_prices(n):
    # WARNING: this is time-averaged, see weighted_prices for load-weighted average
    prices = n.buses_t.marginal_price.mean().groupby(n.buses.carrier).mean()

    return prices


def calculate_weighted_prices(n):
    carriers = n.buses.carrier.unique()

    weighted_prices = {}

    for carrier in carriers:
        load = n.statistics.withdrawal(
            groupby=pypsa.statistics.groupers["bus", "carrier"],
            aggregate_time=False,
            nice_names=False,
            bus_carrier=carrier,
        )

        if not load.empty and load.sum().sum() > 0:
            load = load.groupby(level="bus").sum().T.fillna(0)

            price = n.buses_t.marginal_price.loc[:, n.buses.carrier == carrier]
            price = price.reindex(columns=load.columns, fill_value=1)

            weighted_prices[carrier] = (load * price).sum().sum() / load.sum().sum()

    return pd.Series(weighted_prices).sort_index()


def calculate_market_values(n):
    # Warning: doesn't include storage units

    carrier = "AC"

    market_values = {}

    buses = n.buses.index[n.buses.carrier == carrier]

    ## First do market value of generators ##

    generators = n.generators.index[n.buses.loc[n.generators.bus, "carrier"] == carrier]

    techs = n.generators.loc[generators, "carrier"].value_counts().index

    for tech in techs:
        gens = generators[n.generators.loc[generators, "carrier"] == tech]

        dispatch = (
            n.generators_t.p[gens]
            .T.groupby(n.generators.loc[gens, "bus"])
            .sum()
            .T.reindex(columns=buses, fill_value=0.0)
        )
        revenue = dispatch * n.buses_t.marginal_price[buses]

        if total_dispatch := dispatch.sum().sum():
            market_values[tech] = revenue.sum().sum() / total_dispatch
        else:
            market_values[tech] = np.nan

    ## Now do market value of links ##

    for i in ["0", "1"]:
        all_links = n.links.index[n.buses.loc[n.links["bus" + i], "carrier"] == carrier]

        techs = n.links.loc[all_links, "carrier"].value_counts().index

        for tech in techs:
            links = all_links[n.links.loc[all_links, "carrier"] == tech]

            dispatch = (
                n.links_t["p" + i][links]
                .T.groupby(n.links.loc[links, "bus" + i])
                .sum()
                .T.reindex(columns=buses, fill_value=0.0)
            )

            revenue = dispatch * n.buses_t.marginal_price[buses]

            if total_dispatch := dispatch.sum().sum():
                market_values[tech] = revenue.sum().sum() / total_dispatch
            else:
                market_values[tech] = np.nan

    return pd.Series(market_values).sort_index()


def calculate_price_statistics(n):
    price_statistics = {}

    buses = n.buses.index[n.buses.carrier == "AC"]

    threshold = 0.1  # higher than phoney marginal_cost of wind/solar

    df = pd.DataFrame(data=0.0, columns=buses, index=n.snapshots)

    df[n.buses_t.marginal_price[buses] < threshold] = 1.0

    price_statistics["zero_hours"] = df.sum().sum() / (df.shape[0] * df.shape[1])

    price_statistics["mean"] = n.buses_t.marginal_price[buses].unstack().mean()

    price_statistics["standard_deviation"] = (
        n.buses_t.marginal_price[buses].unstack().std()
    )

    return pd.Series(price_statistics).sort_index()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "make_summary",
            clusters="5",
            opts="",
            sector_opts="",
            planning_horizons="2030",
            configfiles="config/test/config.overnight.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    assign_carriers(n)
    assign_locations(n)

    for output in OUTPUTS:
        globals()["calculate_" + output](n).to_csv(snakemake.output[output])
