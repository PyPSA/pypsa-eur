#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import matplotlib.pyplot as plt
import pandas as pd
import pypsa
import seaborn as sns
from _helpers import configure_logging, set_scenario_config
from pypsa.statistics import get_bus_and_carrier

sns.set_theme("paper", style="whitegrid")

carrier_groups = {
    "Offshore Wind (AC)": "Offshore Wind",
    "Offshore Wind (DC)": "Offshore Wind",
    "Open-Cycle Gas": "Gas",
    "Combined-Cycle Gas": "Gas",
    "Reservoir & Dam": "Hydro",
    "Pumped Hydro Storage": "Hydro",
}


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_validation_electricity_production",
            opts="Ept",
            clusters="37c",
            ll="v1.0",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    n.loads.carrier = "load"

    historic = pd.read_csv(
        snakemake.input.electricity_production,
        index_col=0,
        header=[0, 1],
        parse_dates=True,
    )
    subset_technologies = ["Geothermal", "Nuclear", "Biomass", "Lignite", "Oil", "Coal"]
    lowercase_technologies = [
        technology.lower() if technology in subset_technologies else technology
        for technology in historic.columns.levels[1]
    ]
    historic.columns = historic.columns.set_levels(lowercase_technologies, level=1)

    colors = n.carriers.set_index("nice_name").color.where(
        lambda s: s != "", "lightgrey"
    )
    colors["Offshore Wind"] = colors["Offshore Wind (AC)"]
    colors["Gas"] = colors["Combined-Cycle Gas"]
    colors["Hydro"] = colors["Reservoir & Dam"]
    colors["Other"] = "lightgray"

    if len(historic.index) > len(n.snapshots):
        historic = historic.resample(n.snapshots.inferred_freq).mean().loc[n.snapshots]

    optimized = n.statistics.dispatch(
        groupby=get_bus_and_carrier, aggregate_time=False
    ).T
    optimized = optimized[["Generator", "StorageUnit"]].droplevel(0, axis=1)
    optimized = optimized.rename(columns=n.buses.country, level=0)
    optimized = optimized.rename(columns=carrier_groups, level=1)
    optimized = optimized.T.groupby(level=[0, 1]).sum().T

    data = pd.concat([historic, optimized], keys=["Historic", "Optimized"], axis=1)
    data.columns.names = ["Kind", "Country", "Carrier"]
    data = data.mul(n.snapshot_weightings.generators, axis=0)

    # total production per carrier
    fig, ax = plt.subplots(figsize=(6, 6))

    df = data.groupby(level=["Kind", "Carrier"], axis=1).sum().sum().unstack().T
    df = df / 1e6  # TWh
    df.plot.barh(ax=ax, xlabel="Electricity Production [TWh]", ylabel="")
    ax.grid(axis="y")
    fig.savefig(snakemake.output.production_bar, bbox_inches="tight")
    plt.close(fig)

    # highest diffs

    fig, ax = plt.subplots(figsize=(6, 10))

    df = data.sum() / 1e6  # TWh
    df = df["Optimized"] - df["Historic"]
    df = df.dropna().sort_values()
    df = pd.concat([df.iloc[:5], df.iloc[-5:]])
    c = colors[df.index.get_level_values(1)]
    df.plot.barh(
        xlabel="Optimized Production - Historic Production [TWh]", ax=ax, color=c.values
    )
    ax.set_title("Strongest Deviations")
    ax.grid(axis="y")
    fig.savefig(snakemake.output.production_deviation_bar, bbox_inches="tight")
    plt.close(fig)

    # seasonal operation

    fig, axes = plt.subplots(3, 1, figsize=(9, 9))

    df = (
        data.groupby(level=["Kind", "Carrier"], axis=1)
        .sum()
        .resample("1W")
        .mean()
        .clip(lower=0)
    )
    df = df / 1e3

    order = (
        (df["Historic"].diff().abs().sum() / df["Historic"].sum()).sort_values().index
    )
    c = colors[order]
    optimized = df["Optimized"].reindex(order, axis=1, level=1)
    historical = df["Historic"].reindex(order, axis=1, level=1)

    kwargs = dict(color=c, legend=False, ylabel="Production [GW]", xlabel="")

    optimized.plot.area(ax=axes[0], **kwargs, title="Optimized")
    historical.plot.area(ax=axes[1], **kwargs, title="Historic")

    diff = optimized - historical
    diff.clip(lower=0).plot.area(
        ax=axes[2], **kwargs, title="$\Delta$ (Optimized - Historic)"
    )
    lim = axes[2].get_ylim()[1]
    diff.clip(upper=0).plot.area(ax=axes[2], **kwargs)
    axes[2].set_ylim(bottom=-lim, top=lim)

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(
        h[::-1],
        l[::-1],
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        ncol=1,
        frameon=False,
        labelspacing=1,
    )
    fig.savefig(snakemake.output.seasonal_operation_area, bbox_inches="tight")
    plt.close(fig)

    # touch file
    with open(snakemake.output.plots_touch, "a"):
        pass
