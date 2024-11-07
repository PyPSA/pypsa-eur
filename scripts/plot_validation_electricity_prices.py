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

sns.set_theme("paper", style="whitegrid")

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_electricity_prices",
            opts="Ept-12h",
            clusters="37",
            ll="v1.0",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    n.loads.carrier = "load"

    historic = pd.read_csv(
        snakemake.input.electricity_prices,
        index_col=0,
        header=0,
        parse_dates=True,
    )

    if len(historic.index) > len(n.snapshots):
        historic = historic.resample(n.snapshots.inferred_freq).mean().loc[n.snapshots]

    optimized = n.buses_t.marginal_price.groupby(n.buses.country, axis=1).mean()

    data = pd.concat([historic, optimized], keys=["Historic", "Optimized"], axis=1)
    data.columns.names = ["Kind", "Country"]

    fig, ax = plt.subplots(figsize=(6, 6))

    df = data.mean().unstack().T
    df.plot.barh(ax=ax, xlabel="Electricity Price [€/MWh]", ylabel="")
    ax.grid(axis="y")
    fig.savefig(snakemake.output.price_bar, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots()

    df = data.groupby(level="Kind", axis=1).mean()
    df.plot(ax=ax, xlabel="", ylabel="Electricity Price [€/MWh]", alpha=0.8)
    ax.grid(axis="x")
    fig.savefig(snakemake.output.price_line, bbox_inches="tight")
    plt.close(fig)

    # touch file
    with open(snakemake.output.plots_touch, "a"):
        pass
