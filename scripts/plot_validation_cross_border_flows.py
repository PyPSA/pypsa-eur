#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import country_converter as coco
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
import seaborn as sns
from _helpers import configure_logging

sns.set_theme("paper", style="whitegrid")

cc = coco.CountryConverter()

color_country = {
    "AL": "#440154",
    "AT": "#482677",
    "BA": "#43398e",
    "BE": "#3953a4",
    "BG": "#2c728e",
    "CH": "#228b8d",
    "CZ": "#1f9d8a",
    "DE": "#29af7f",
    "DK": "#3fbc73",
    "EE": "#5ec962",
    "ES": "#84d44b",
    "FI": "#addc30",
    "FR": "#d8e219",
    "GB": "#fde725",
    "GR": "#f0f921",
    "HR": "#f1c25e",
    "HU": "#f4a784",
    "IE": "#f78f98",
    "IT": "#f87ea0",
    "LT": "#f87a9a",
    "LU": "#f57694",
    "LV": "#f3758d",
    "ME": "#f37685",
    "MK": "#f37b7c",
    "NL": "#FF6666",
    "NO": "#FF3333",
    "PL": "#eb0000",
    "PT": "#d70000",
    "RO": "#c00000",
    "RS": "#a50000",
    "SE": "#8a0000",
    "SI": "#6f0000",
    "SK": "#550000",
}


def sort_one_country(country, df):
    indices = [link for link in df.columns if country in link]
    df_country = df[indices].copy()
    for link in df_country.columns:
        if country in link[5:]:
            df_country[link] = -df_country[link]
            link_reverse = str(link[5:] + " - " + link[:2])
            df_country = df_country.rename(columns={link: link_reverse})

    return df_country.reindex(sorted(df_country.columns), axis=1)


def cross_border_time_series(countries, data):
    fig, ax = plt.subplots(2 * len(countries), 1, figsize=(15, 10 * len(countries)))
    axis = 0

    for country in countries:
        ymin = 0
        ymax = 0
        for df in data:
            df_country = sort_one_country(country, df)
            df_neg, df_pos = df_country.clip(upper=0), df_country.clip(lower=0)

            color = [color_country[link[5:]] for link in df_country.columns]

            df_pos.plot.area(
                ax=ax[axis], stacked=True, linewidth=0.0, color=color, ylim=[-1, 1]
            )

            df_neg.plot.area(
                ax=ax[axis], stacked=True, linewidth=0.0, color=color, ylim=[-1, 1]
            )
            title = "Historic" if (axis % 2) == 0 else "Optimized"
            ax[axis].set_title(
                f"{title} Import / Export for " + cc.convert(country, to="name_short")
            )

            # Custom legend elements
            legend_elements = []

            for link in df_country.columns:
                legend_elements = legend_elements + [
                    plt.fill_between(
                        [],
                        [],
                        color=color_country[link[5:]],
                        label=cc.convert(link[5:], to="name_short"),
                    )
                ]

            # Create the legend
            ax[axis].legend(handles=legend_elements, loc="upper right")

            # rescale the y axis
            neg_min = df_neg.sum(axis=1).min() * 1.2
            if neg_min < ymin:
                ymin = neg_min

            pos_max = df_pos.sum(axis=1).max() * 1.2
            if pos_max < ymax:
                ymax = pos_max

            axis = axis + 1

        for x in range(axis - 2, axis):
            ax[x].set_ylim([neg_min, pos_max])

    fig.savefig(snakemake.output.trade_time_series, bbox_inches="tight")


def cross_border_bar(countries, data):
    df_positive = pd.DataFrame()
    df_negative = pd.DataFrame()
    color = []

    for country in countries:
        order = 0
        for df in data:
            df_country = sort_one_country(country, df)
            df_neg, df_pos = df_country.clip(upper=0), df_country.clip(lower=0)

            title = "Historic" if (order % 2) == 0 else "Optimized"
            df_positive_new = pd.DataFrame(data=df_pos.sum()).T.rename(
                {0: f"{title} " + cc.convert(country, to="name_short")}
            )
            df_negative_new = pd.DataFrame(data=df_neg.sum()).T.rename(
                {0: f"{title} " + cc.convert(country, to="name_short")}
            )

            df_positive = pd.concat([df_positive_new, df_positive])
            df_negative = pd.concat([df_negative_new, df_negative])

            order = order + 1

    color = [color_country[link[5:]] for link in df_positive.columns]

    fig, ax = plt.subplots(figsize=(15, 60))

    df_positive.plot.barh(ax=ax, stacked=True, color=color, zorder=2)
    df_negative.plot.barh(ax=ax, stacked=True, color=color, zorder=2)

    plt.grid(axis="x", zorder=0)
    plt.grid(axis="y", zorder=0)

    # Custom legend elements
    legend_elements = []

    for country in list(color_country.keys()):
        legend_elements = legend_elements + [
            plt.fill_between(
                [],
                [],
                color=color_country[country],
                label=cc.convert(country, to="name_short"),
            )
        ]

    # Create the legend
    plt.legend(handles=legend_elements, loc="upper right")

    fig.savefig(snakemake.output.cross_border_bar, bbox_inches="tight")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_electricity_prices",
            simpl="",
            opts="Ept-12h",
            clusters="37",
            ll="v1.0",
        )
    configure_logging(snakemake)

    countries = snakemake.params.countries

    n = pypsa.Network(snakemake.input.network)
    n.loads.carrier = "load"

    historic = pd.read_csv(
        snakemake.input.cross_border_flows,
        index_col=0,
        header=0,
        parse_dates=True,
    )

    if len(historic.index) > len(n.snapshots):
        historic = historic.resample(n.snapshots.inferred_freq).mean().loc[n.snapshots]

    # Preparing network data to be shaped similar to ENTSOE datastructure
    optimized_links = n.links_t.p0.rename(
        columns=dict(n.links.bus0.str[:2] + " - " + n.links.bus1.str[:2])
    )
    optimized_lines = n.lines_t.p0.rename(
        columns=dict(n.lines.bus0.str[:2] + " - " + n.lines.bus1.str[:2])
    )
    optimized = pd.concat([optimized_links, optimized_lines], axis=1)

    # Drop internal country connection
    optimized.drop(
        [c for c in optimized.columns if c[:2] == c[5:]], axis=1, inplace=True
    )

    # align columns name
    for c1 in optimized.columns:
        for c2 in optimized.columns:
            if c1[:2] == c2[5:] and c2[:2] == c1[5:]:
                optimized = optimized.rename(columns={c1: c2})

    optimized = optimized.groupby(lambda x: x, axis=1).sum()

    cross_border_bar(countries, [historic, optimized])

    cross_border_time_series(countries, [historic, optimized])

    # touch file
    with open(snakemake.output.plots_touch, "a"):
        pass
