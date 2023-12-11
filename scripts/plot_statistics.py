#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import matplotlib.pyplot as plt
import pypsa
import seaborn as sns
from _helpers import configure_logging

sns.set_theme("paper", style="whitegrid")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_elec_statistics",
            simpl="",
            opts="Ept-12h",
            clusters="37",
            ll="v1.0",
        )
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    n.loads.carrier = "load"
    n.carriers.loc["load", ["nice_name", "color"]] = "Load", "darkred"
    colors = n.carriers.set_index("nice_name").color.where(
        lambda s: s != "", "lightgrey"
    )

    def rename_index(ds):
        specific = ds.index.map(lambda x: f"{x[1]}\n({x[0]})")
        generic = ds.index.get_level_values("carrier")
        duplicated = generic.duplicated(keep=False)
        index = specific.where(duplicated, generic)
        return ds.set_axis(index)

    def plot_static_per_carrier(ds, ax, drop_zero=True):
        if drop_zero:
            ds = ds[ds != 0]
        ds = ds.dropna()
        c = colors[ds.index.get_level_values("carrier")]
        ds = ds.pipe(rename_index)
        label = f"{ds.attrs['name']} [{ds.attrs['unit']}]"
        ds.plot.barh(color=c.values, xlabel=label, ax=ax)
        ax.grid(axis="y")

    fig, ax = plt.subplots()
    ds = n.statistics.capacity_factor().dropna()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.capacity_factor_bar)

    fig, ax = plt.subplots()
    ds = n.statistics.installed_capacity().dropna()
    ds = ds.drop("Line")
    ds = ds.drop(("Generator", "Load"))
    ds = ds / 1e3
    ds.attrs["unit"] = "GW"
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.installed_capacity_bar)

    fig, ax = plt.subplots()
    ds = n.statistics.optimal_capacity()
    ds = ds.drop("Line")
    ds = ds.drop(("Generator", "Load"))
    ds = ds / 1e3
    ds.attrs["unit"] = "GW"
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.optimal_capacity_bar)

    fig, ax = plt.subplots()
    ds = n.statistics.capex()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.capital_expenditure_bar)

    fig, ax = plt.subplots()
    ds = n.statistics.opex()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.operational_expenditure_bar)

    fig, ax = plt.subplots()
    ds = n.statistics.curtailment()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.curtailment_bar)

    fig, ax = plt.subplots()
    ds = n.statistics.supply()
    ds = ds.drop("Line")
    ds = ds / 1e6
    ds.attrs["unit"] = "TWh"
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.supply_bar)

    fig, ax = plt.subplots()
    ds = n.statistics.withdrawal()
    ds = ds.drop("Line")
    ds = ds / -1e6
    ds.attrs["unit"] = "TWh"
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.withdrawal_bar)

    fig, ax = plt.subplots()
    ds = n.statistics.market_value()
    plot_static_per_carrier(ds, ax)
    fig.savefig(snakemake.output.market_value_bar)

    # touch file
    with open(snakemake.output.barplots_touch, "a"):
        pass
