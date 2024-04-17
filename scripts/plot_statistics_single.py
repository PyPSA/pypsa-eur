#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import logging
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from _helpers import configure_logging
from plot_summary import rename_techs

logger = logging.getLogger(__name__)
sns.set_theme("paper", style="whitegrid")


def rename_index(ds):
    specific = ds.index.map(lambda x: f"{x[1]}\n({x[0]})")
    generic = ds.index.get_level_values("carrier")
    duplicated = generic.duplicated(keep=False)
    index = specific.where(duplicated, generic)
    return ds.set_axis(index)


def plot_static_single(ds, ax):
    factor, unit = conversion[output]
    ds = ds.dropna()
    carriers = ds.index.get_level_values("carrier").map(rename_techs)
    if not carriers.difference(tech_colors.index).empty:
        print(
            f"Missing colors for carrier: {carriers.difference(tech_colors.index).values}\n Dark grey used instead."
        )
    c = carriers.map(lambda x: tech_colors.get(x, "#808080"))
    ds = ds.pipe(rename_index)
    ds = ds.div(float(factor)) if factor != "-" else ds
    ds.T.plot.barh(color=c.values, ax=ax, xlabel=unit)
    ax.grid(axis="y")


def read_csv(input):
    try:
        df = pd.read_csv(input, skiprows=2)
        df = df.set_index(["component", "carrier"]).iloc[:, 0]
    except Exception as e:
        print(f"An error occurred reading file {input}: {e}")
        df = pd.DataFrame()
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_statistics_single",
            run = "240219-test/normal",
            simpl="",
            ll="v1.2",
            clusters="22",
            opts="",
            sector_opts="none",
            planning_horizons="2020",
            country="DE",
            carrier="heat",
        )
    configure_logging(snakemake)

    tech_colors = pd.Series(snakemake.params.plotting["tech_colors"])
    conversion = pd.Series(snakemake.params.statistics)

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
        plot_static_single(ds, ax)
        fig.savefig(snakemake.output[output], bbox_inches="tight")
