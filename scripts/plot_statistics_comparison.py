#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from _helpers import configure_logging
from plot_summary import rename_techs

logger = logging.getLogger(__name__)
sns.set_theme("paper", style="whitegrid")
STACKED = {
    "capacity_factor": False,
    "market_value": False,
}


def rename_index(df):
    # rename index and drop duplicates
    specific = df.index.map(lambda x: f"{x[1]}({x[0]})")
    generic = df.index.get_level_values("carrier")
    duplicated = generic.duplicated(keep=False)
    index = specific.where(duplicated, generic)
    df = df.set_axis(index)
    # rename columns and drop duplicates
    opts = df.columns.get_level_values(-1).str.split("_", expand=True)
    scenario = df.columns.get_level_values(0)
    opts = [
        opts.get_level_values(level).unique()
        for level in range(opts.nlevels)
        if not opts.get_level_values(level).duplicated(keep=False).all()
    ]
    opts = [""] if not opts else opts  # in case all opts are the same
    if df.columns.nlevels > 1:
        columns = pd.MultiIndex.from_product([scenario, opts])
    # columns = pd.MultiIndex.from_arrays(columns)
    df.columns = columns.map("\n".join)
    return df


def get_carrier_colors(carriers, tech_colors):
    if not carriers.difference(tech_colors.index).empty:
        print(
            f"Missing colors for carrier: {carriers.difference(tech_colors.index).values}\n Dark grey used instead."
        )
    carrier_colors = carriers.map(lambda x: tech_colors.get(x, "#808080")).values
    return carrier_colors


def plot_static_comparison(df, ax, stacked=False, tech_colors=None):
    factor, unit = conversion[output]
    df = df.round(2)
    df = df[df != 0]
    df = df.dropna(axis=0, how="all").fillna(0)
    if df.empty:
        return
    df = df.div(float(factor)) if factor != "-" else df
    df = df.rename(index=rename_techs).groupby(["component", "carrier"]).sum()
    # sort values in descending order
    df = df.reindex(df.abs().sum(1).sort_values().index)
    carriers = df.index.get_level_values("carrier")
    carrier_colors = get_carrier_colors(carriers, tech_colors)
    df = df.pipe(rename_index).T
    df.plot.bar(
        color=carrier_colors,
        ax=ax,
        stacked=stacked,
        legend=False,
        ylabel=unit,
        linewidth=0.1,
    )
    ax.legend(
        bbox_to_anchor=(1, 1),
        loc="upper left",
    )
    ax.grid(axis="x")


def read_csv(input, output, num_run):
    try:
        # filter required csv to plot the wanted output
        files = list(filter(lambda x: output in x, input))
        # retrieves network labels from folder name
        if num_run == 1:
            network_labels = [file.split("/")[-3] for file in files]
        else:
            network_labels = [
                (file.split("/")[-6], file.split("/")[-3]) for file in files
            ]
        df = pd.concat(
            [
                pd.read_csv(f, skiprows=2).set_index(["component", "carrier"])
                for f in files
            ],
            axis=1,
            keys=network_labels,
            names=["scenario", "opts", "to_drop"],
        )
        df.columns = df.columns.droplevel("to_drop")
    except Exception as e:
        print(f"Error reading csv file for {output}: {e}")
        df = pd.DataFrame()
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_statistics_comparison",
            run="",
            country="all",
            carrier="electricity",
        )

    configure_logging(snakemake)

    plotting = snakemake.params.plotting
    tech_colors = pd.Series(plotting["tech_colors"]).groupby(rename_techs).first()
    conversion = pd.Series(snakemake.params.statistics)

    all_runs = snakemake.config["run"]["name"]
    runs = plotting["statistics"].get("scenario_comparison", all_runs)
    num_run = len(np.atleast_1d(runs))

    for output in snakemake.output.keys():
        if "touch" in output:
            with open(snakemake.output.barplots_touch, "a"):
                pass
                continue
        fig, ax = plt.subplots()
        df = read_csv(snakemake.input, output, num_run)
        if df.empty:
            ax.text(
                0.5,
                0.5,
                "No data available.",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=14,
                color="red",
            )
            fig.savefig(snakemake.output[output])
            continue

        plot_static_comparison(
            df, ax, stacked=STACKED.get(output, True), tech_colors=tech_colors
        )

        fig.savefig(snakemake.output[output], bbox_inches="tight")
