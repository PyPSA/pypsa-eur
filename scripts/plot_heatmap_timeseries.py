# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot heatmap time series (results).
"""

import logging
import pypsa
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool

logger = logging.getLogger(__name__)

THRESHOLD_MW = 1e3 # 1 GW
THRESHOLD_MWh = 100e3 # 100 GWh

MARGINAL_PRICES = [
    "AC",
    "H2",
    "NH3",
    "gas",
    "co2 stored",
    "methanol",
    "oil",
    "rural heat",
    "urban central heat",
    "urban decentral heat"
]

SKIP_UTILISATION_RATES = [
    'DC',
    'H2 pipeline',
    'electricity distribution grid',
    'gas for industry',
    'gas for industry CC',
    'gas pipeline',
    'gas pipeline new',
    'process emissions',
    'process emissions CC',
    'solid biomass for industry',
    'solid biomass for industry CC',
]


def unstack_day_hour(s, sns):
    if isinstance(sns, dict):
        sns = pd.date_range(freq="h", **sns)
    s_h = s.reindex(sns).ffill()
    grouped = s_h.groupby(s_h.index.hour).agg(list)
    index = [f"{i:02d}:00" for i in grouped.index]
    columns = pd.date_range(s_h.index[0], s_h.index[-1], freq="D")
    return pd.DataFrame(grouped.to_list(), index=index, columns=columns)


def plot_heatmap(
    df,
    vmin=None,
    vmax=None,
    center=None,
    cmap="Greens",
    label="",
    title="",
    cbar_kws={},
    tag="",
    dir="",
):
    _cbar_kws = dict(label=label, aspect=17, pad=0.015)
    _cbar_kws.update(cbar_kws)
    fig, ax = plt.subplots(figsize=(8.5, 4), constrained_layout=True)
    sns.heatmap(
        df,
        cmap=cmap,
        ax=ax,
        vmin=vmin,
        vmax=vmax,
        center=center,
        cbar_kws=_cbar_kws,
    )
    plt.ylabel("hour of the day")
    plt.xlabel("day of the year")
    plt.title(title, fontsize='large')

    ax.grid(axis='y')

    hours = list(range(0,24))
    ax.set_yticks(hours[0::2])
    ax.set_yticklabels(df.index[0::2], rotation=0)
    ax.set_yticks(hours, minor=True)

    major_ticks = [i for i, date in enumerate(df.columns) if date.day == 1]
    minor_ticks = [i for i, date in enumerate(df.columns) if date.day == 15]
    ax.set_xticks(major_ticks)
    ax.set_xticklabels([df.columns[i].strftime('%e\n%b') for i in major_ticks], rotation=0)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticklabels([df.columns[i].strftime('%e') for i in minor_ticks], rotation=0, minor=True, color='grey')

    cb = ax.collections[0].colorbar
    cb.outline.set_linewidth(0.75)

    for spine in ax.spines.values():
        spine.set_visible(True)

    title = title.lower().replace(" ", "_").replace("(", "").replace(")", "")
    fn = f"ts-heatmap-{tag}-{title}"
    plt.savefig(dir + "/" + fn + ".pdf")
    plt.savefig(dir + "/" + fn + ".png")
    plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_heatmap_timeseries",
            simpl="",
            clusters=100,
            ll="v1.5",
            opts="",
            sector_opts="Co2L0-2190SEG-T-H-B-I-S-A",
            planning_horizons=2050,
            configfiles="../../config/config.100n-seg.yaml"
        )

    plt.style.use(["bmh", snakemake.input.rc])


    n = pypsa.Network(snakemake.input.network)

    dir = snakemake.output[0]
    snapshots = snakemake.config["snapshots"]
    carriers = n.carriers

    p_nom_opt = n.generators.groupby("carrier").p_nom_opt.sum()
    data = (
        n.generators_t.p.groupby(n.generators.carrier, axis=1).sum()
        / p_nom_opt * 100 # %
    )
    data = data.loc[:, p_nom_opt > THRESHOLD_MW]

    def process_generator_utilisation(carrier, s):
        df = unstack_day_hour(s, snapshots)
        title = carriers.nice_name.get(carrier, carrier)
        label = "utilisation rate [%]"
        plot_heatmap(
            df,
            cmap="Blues",
            label=label,
            title=title,
            vmin=0,
            vmax=90,
            cbar_kws=dict(extend='max'),
            tag="utilisation_rate",
            dir=dir,
        )

    with Pool(processes=snakemake.threads) as pool:
        pool.starmap(process_generator_utilisation, data.items())

    # marginal prices
    data = n.buses_t.marginal_price.groupby(n.buses.carrier, axis=1).mean()
    data = data[data.columns.intersection(MARGINAL_PRICES)].round(2)

    def process_marginal_prices(carrier, s):
        df = unstack_day_hour(s, snapshots)
        label = "marginal price [€/t]" if "co2" in carrier.lower() else "marginal price [€/MWh]"
        plot_heatmap(
            df,
            cmap="Spectral_r",
            label=label,
            title=carrier,
            cbar_kws=dict(extend='both'),
            tag="marginal_price",
            dir=dir,
        )

    with Pool(processes=snakemake.threads) as pool:
        pool.starmap(process_marginal_prices, data.items())

    # SOCs
    e_nom_opt = n.stores.groupby("carrier").e_nom_opt.sum()
    data = (
        n.stores_t.e.groupby(n.stores.carrier, axis=1).sum()
        / e_nom_opt * 100
    )
    data = data.loc[:, e_nom_opt > THRESHOLD_MWh]

    def process_socs(carrier, s):
        df = unstack_day_hour(s, snapshots)
        label = "SOC [%]"
        title = carriers.nice_name.get(carrier, carrier)
        plot_heatmap(
            df,
            cmap="Purples",
            vmin=0,
            vmax=90,
            label=label,
            title=title,
            cbar_kws=dict(extend='max'),
            tag="soc",
            dir=dir,
        )

    with Pool(processes=snakemake.threads) as pool:
        pool.starmap(process_socs, data.items())

    # link utilisation rates
    p_nom_opt = n.links.groupby("carrier").p_nom_opt.sum()
    data = (
        n.links_t.p0.groupby(n.links.carrier, axis=1).sum()
        / p_nom_opt * 100
    )
    data = data[data.columns.difference(SKIP_UTILISATION_RATES)]
    data = data.loc[:,p_nom_opt > THRESHOLD_MW]

    def process_link_utilisation(carrier, s):
        df = unstack_day_hour(s, snapshots)
        label = "utilisation rate [%]"
        title = carriers.nice_name.get(carrier, carrier)
        plot_heatmap(
            df,
            cmap="Reds",
            vmin=0,
            vmax=100,
            label=label,
            title=title,
            tag="utilisation_rate",
            dir=dir,
        )

    with Pool(processes=snakemake.threads) as pool:
        pool.starmap(process_link_utilisation, data.items())
