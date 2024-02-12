#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import re

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns
from _helpers import configure_logging
from plot_summary import rename_techs

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
    columns = df.columns.str.split("_", expand=True)
    columns = [
        columns.get_level_values(level).unique()
        for level in range(columns.nlevels)
        if not columns.get_level_values(level).duplicated(keep=False).all()
    ]
    columns = pd.MultiIndex.from_arrays(columns)
    df.columns = columns.map("\n".join)
    return df


def plot_static_comparison(df, ax, stacked=False):
    df = df[df != 0]
    df = df.dropna(axis=0, how="all").fillna(0)
    c = tech_colors[df.index.get_level_values("carrier").map(rename_techs)]
    df = df.pipe(rename_index).T
    df.plot.bar(color=c.values, ax=ax, stacked=stacked, legend=False)
    ax.legend(
        bbox_to_anchor=(1, 1),
        loc="upper left",
    )
    ax.grid(axis="x")


def read_csv(input, output):
    try:
        # filter required csv to plot the wanted output
        files = list(filter(lambda x: output in x, input))
        pattern = r"elec_.*?(\d{4})"
        network_labels = [re.search(pattern, f).group() for f in files]
        df = pd.concat(
            [pd.read_csv(f).set_index(["component", "carrier"]) for f in files],
            axis=1,
            keys=network_labels,
        )
        # get plot label and drop from index
        label = df.columns.get_level_values(1).unique()[0]
        df.columns = df.columns.droplevel(1)
    except Exception as e:
        print(f"Error reading csv file for {output}: {e}")
        df = pd.DataFrame()
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_statistics_comparison",
            country="all",
            carrier="electricity",
        )
    configure_logging(snakemake)

    tech_colors = pd.Series(snakemake.params.plotting["tech_colors"])

    for output in snakemake.output.keys():
        if "touch" in output:
            with open(snakemake.output.barplots_touch, "a"):
                pass
                continue
        fig, ax = plt.subplots()
        if "energy_balance" not in output:
            df = read_csv(snakemake.input, output)
        else:
            supply = read_csv(snakemake.input, "supply")
            withdrawal = read_csv(snakemake.input, "withdrawal")
            df = (
                pd.concat([supply, withdrawal.mul(-1)])
                .groupby(["component", "carrier"])
                .sum()
            )
        if df.empty:
            fig.savefig(snakemake.output[output])
            continue

        plot_static_comparison(df, ax, stacked=STACKED.get(output, True))

        fig.savefig(snakemake.output[output], bbox_inches="tight")
