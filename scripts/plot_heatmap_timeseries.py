# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Plot heatmap time series of marginal prices, utilisation rates, state of charge profiles.
"""

import logging
import os

import matplotlib.pyplot as plt
import pandas as pd
import pypsa
import seaborn as sns
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def unstack_day_hour(s: pd.Series, sns: pd.DatetimeIndex) -> pd.DataFrame:
    if isinstance(sns, dict):
        sns = pd.date_range(freq="h", **sns)
    s_h = s.reindex(sns).ffill()
    grouped = s_h.groupby(s_h.index.hour).agg(list)
    index = [f"{i:02d}:00" for i in grouped.index]
    columns = pd.date_range(s_h.index[0], s_h.index[-1], freq="D")
    return pd.DataFrame(grouped.to_list(), index=index, columns=columns)


def plot_heatmap(
    df: pd.DataFrame,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str = "Greens",
    label: str = "",
    title: str = "",
    cbar_kws: dict = {},
    fn: str | None = None,
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
        cbar_kws=_cbar_kws,
    )
    plt.ylabel("hour of the day")
    plt.xlabel("day of the year")
    plt.title(title, fontsize="large")

    ax.grid(axis="y")

    hours = list(range(0, 24))
    ax.set_yticks(hours[0::2])
    ax.set_yticklabels(df.index[0::2], rotation=0)
    ax.set_yticks(hours, minor=True)

    major_ticks = [i for i, date in enumerate(df.columns) if date.day == 1]
    minor_ticks = [i for i, date in enumerate(df.columns) if date.day == 15]
    ax.set_xticks(major_ticks)
    ax.set_xticklabels(
        [df.columns[i].strftime("%e\n%b") for i in major_ticks], rotation=0
    )
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticklabels(
        [df.columns[i].strftime("%e") for i in minor_ticks],
        rotation=0,
        minor=True,
        color="grey",
    )

    cb = ax.collections[0].colorbar
    cb.outline.set_linewidth(0)

    if fn is not None:
        plt.savefig(fn)
        plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_heatmap_timeseries",
            simpl="",
            clusters="10",
            ll="vopt",
            opts="",
            sector_opts="",
            planning_horizons=2050,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    plt.style.use(["bmh", snakemake.input.rc])

    config = snakemake.params.plotting["heatmap_timeseries"]

    output_dir = snakemake.output[0]
    os.makedirs(output_dir, exist_ok=True)

    n = pypsa.Network(snakemake.input.network)

    snapshots = snakemake.params["snapshots"]
    carriers = n.carriers

    # filter for build capacities
    optimal_capacity = (
        n.statistics.optimal_capacity(nice_names=False).groupby("carrier").sum()
    )
    built_idx = optimal_capacity.where(optimal_capacity > 100).dropna().index

    # utilisation rates
    cf = (
        n.statistics.capacity_factor(aggregate_time=False, nice_names=False)
        .dropna()
        .groupby("carrier")
        .sum()
        .mul(100)
    )
    idx = cf.index.intersection(config["utilisation_rate"]).intersection(built_idx)
    cf = cf.loc[idx]

    for carrier, s in cf.iterrows():
        logger.info(f"Plotting utilisation rate heatmap time series for {carrier}")
        df = unstack_day_hour(s, snapshots)
        label = "utilisation rate [%]"
        fn = (
            output_dir
            + "/ts-heatmap-utilisation_rate-"
            + carrier.replace(" ", "_")
            + ".pdf"
        )
        plot_heatmap(
            df,
            cmap="Greens",
            label=label,
            title=carrier,
            vmin=0,
            vmax=100,
            fn=fn,
        )

    # marginal prices
    prices = n.buses_t.marginal_price.T.groupby(n.buses.carrier).mean()
    prices = prices.loc[prices.index.intersection(config["marginal_price"])]

    for carrier, s in prices.iterrows():
        logger.info(f"Plotting marginal prices heatmap time series for {carrier}")
        df = unstack_day_hour(s, snapshots)
        label = (
            "marginal price [€/t]"
            if "co2" in carrier.lower()
            else "marginal price [€/MWh]"
        )
        fn = (
            output_dir
            + "/ts-heatmap-marginal_price-"
            + carrier.replace(" ", "_")
            + ".pdf"
        )
        plot_heatmap(
            df,
            cmap="Spectral_r",
            label=label,
            title=carrier,
            fn=fn,
        )

    # SOCs
    e_nom_opt = n.stores.groupby("carrier").e_nom_opt.sum()
    socs = n.stores_t.e.T.groupby(n.stores.carrier).sum().div(e_nom_opt, axis=0) * 100
    socs = socs.loc[socs.index.intersection(config["soc"]).intersection(built_idx)]

    for carrier, s in socs.iterrows():
        logger.info(f"Plotting SOC heatmap time series for {carrier}")
        df = unstack_day_hour(s, snapshots)
        label = "SOC [%]"
        fn = output_dir + "/ts-heatmap-soc-" + carrier.replace(" ", "_") + ".pdf"
        plot_heatmap(
            df,
            cmap="Blues",
            vmin=0,
            vmax=100,
            label=label,
            title=carrier,
            fn=fn,
        )
