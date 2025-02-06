# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Create summary CSV files for all scenario runs including costs, capacities,
capacity factors, curtailment, energy balances, prices and other metrics.
"""

import logging

import pandas as pd
import pypsa
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)
OPT_NAME = {"Store": "e", "Line": "s", "Transformer": "s"}


def assign_carriers(n):
    if "carrier" not in n.lines:
        n.lines["carrier"] = "AC"


def assign_locations(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.unique():
            names = ifind.index[ifind == i]
            c.df.loc[names, "location"] = "" if i == -1 else names.str[:i]


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

    return pd.concat(costs)


def calculate_supply_energy(n):
    """
    Calculate the total energy supply/consuption of each component at the buses
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

                if len(items) == 0:
                    continue

                s = (-1) * c.pnl["p" + end][items].multiply(
                    n.snapshot_weightings.generators, axis=0
                ).sum().groupby(c.df.loc[items, "carrier"]).sum()
                s.index = s.index + end
                s = pd.concat([s], keys=[c.list_name])
                s = pd.concat([s], keys=[i])

                supply_energy.append(s)

    return pd.concat(supply_energy)


def calculate_metrics(n):
    metrics = {}

    metrics["line_volume_DC"] = (n.links.length * n.links.p_nom_opt)[
        n.links.carrier == "DC"
    ].sum()
    metrics["line_volume_AC"] = (n.lines.length * n.lines.s_nom_opt).sum()
    metrics["line_volume"] = metrics["line_volume_AC"] + metrics["line_volume_DC"]

    if "lv_limit" in n.global_constraints.index:
        metrics["line_volume_limit"] = n.global_constraints.at[
            "lv_limit", "constant"
        ]
        metrics["line_volume_shadow"] = n.global_constraints.at[
            "lv_limit", "mu"
        ]

    if "CO2Limit" in n.global_constraints.index:
        metrics["co2_shadow"] = n.global_constraints.at["CO2Limit", "mu"]

    if "co2_sequestration_limit" in n.global_constraints.index:
        metrics["co2_storage_shadow"] = n.global_constraints.at[
            "co2_sequestration_limit", "mu"
        ]

    return pd.Series(metrics)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "make_summary",
            opts="",
            clusters="115",
            ll="vopt",
            sector_opts="imp-CF+costs+year+2020",
            planning_horizons="2050",
            configfiles="../../config/config.20240826-z1.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)

    assign_carriers(n)
    assign_locations(n)

    costs = calculate_costs(n)
    costs.round(3).to_csv(snakemake.output.costs)

    supply_energy = calculate_supply_energy(n)
    supply_energy.round(3).to_csv(snakemake.output.supply_energy)

    metrics = calculate_metrics(n)
    metrics.loc["total costs"] = costs.sum()
    metrics.to_csv(snakemake.output.metrics)
