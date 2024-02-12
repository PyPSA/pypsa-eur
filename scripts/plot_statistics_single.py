#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from _helpers import configure_logging
from plot_summary import rename_techs
from pypsa.statistics import get_carrier

sns.set_theme("paper", style="whitegrid")


def rename_index(ds):
    specific = ds.index.map(lambda x: f"{x[1]}\n({x[0]})")
    generic = ds.index.get_level_values("carrier")
    duplicated = generic.duplicated(keep=False)
    index = specific.where(duplicated, generic)
    return ds.set_axis(index)


def plot_static_per_carrier(ds, ax):
    ds = ds.dropna()
    c = tech_colors[ds.index.get_level_values("carrier").map(rename_techs)]
    ds = ds.pipe(rename_index)
    ds.T.plot.barh(color=c.values, ax=ax)
    ax.grid(axis="x")


def read_csv(input):
    try:
        df = pd.read_csv(input, skiprows=2)
        df = df.set_index(["component", "carrier"]).squeeze()
    except Exception as e:
        print(f"An error occurred reading file {input}: {e}")
        df = pd.DataFrame()
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_statistics_single",
            simpl="",
            ll="v1.5",
            clusters="5",
            opts="",
            sector_opts="24h-T-H-B-I-A-dist1",
            planning_horizons="2040",
            country="all",
            carrier="heat",
        )
    configure_logging(snakemake)

    tech_colors = pd.Series(snakemake.params.plotting["tech_colors"])

    for output in snakemake.output.keys():
        if "touch" in output:
            with open(snakemake.output.barplots_touch, "a"):
                pass
            continue
        fig, ax = plt.subplots()
        ds = read_csv(snakemake.input[output])
        if ds.empty:
            fig.savefig(snakemake.output[output])
            continue
        plot_static_per_carrier(ds, ax)
        fig.savefig(snakemake.output[output], bbox_inches="tight")
