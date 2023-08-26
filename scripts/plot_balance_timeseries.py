# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot balance time series.
"""

import logging
import pypsa
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from multiprocessing import Pool

logger = logging.getLogger(__name__)

THRESHOLD = 5 # GW

CARRIER_GROUPS = {
    "electricity": ["AC", "low voltage"],
    "heat": ["urban central heat", "urban decentral heat", "rural heat", "residential urban decentral heat", "residential rural heat", "services urban decentral heat", "services rural heat"],
    "hydrogen": "H2",
    "oil": "oil",
    "methanol": "methanol",
    "ammonia": "NH3",
    "biomass": ["solid biomass", "biogas"],
    "CO2 atmosphere": "co2",
    "CO2 stored": "co2 stored",
    "methane": "gas",
}

def plot_stacked_area_steplike(ax, df, colors={}):

    if isinstance(colors, pd.Series):
        colors = colors.to_dict()

    df_cum = df.cumsum(axis=1)

    previous_series = np.zeros_like(df_cum.iloc[:, 0].values)

    for col in df_cum.columns:
        ax.fill_between(
            df_cum.index,
            previous_series,
            df_cum[col],
            step='pre',
            linewidth=0,
            color=colors.get(col, 'grey'),
            label=col,
        )
        previous_series = df_cum[col].values


def plot_energy_balance_timeseries(
    df,
    time=None,
    ylim=None,
    resample=None,
    rename={},
    preferred_order=[],
    ylabel="",
    colors={},
    threshold=0,
    dir="",
):

    if time is not None:
        df = df.loc[time]

    timespan = (df.index[-1] - df.index[0])
    long_time_frame = timespan > pd.Timedelta(weeks=5)

    techs_below_threshold = df.columns[df.abs().max() < threshold].tolist()
    if techs_below_threshold:
        other = {tech: "other" for tech in techs_below_threshold}
        rename.update(other)
        colors["other"] = 'grey'
    
    if rename:
        df = df.groupby(df.columns.map(lambda a: rename.get(a, a)), axis=1).sum()

    if resample is not None:
        # upsampling to hourly resolution required to handle overlapping block
        df = df.resample("1H").ffill().resample(resample).mean()

    order = (df / df.max()).var().sort_values().index
    if preferred_order:
        order = preferred_order.intersection(order).append(
            order.difference(preferred_order)
        )
    df = df.loc[:, order]

    # fillna since plot_stacked_area_steplike cannot deal with NaNs
    pos = df.where(df > 0).fillna(0.)
    neg = df.where(df < 0).fillna(0.)

    fig, ax = plt.subplots(figsize=(10, 4), layout="constrained")

    plot_stacked_area_steplike(ax, pos, colors)
    plot_stacked_area_steplike(ax, neg, colors)

    plt.xlim((df.index[0], df.index[-1]))

    if not long_time_frame:
        # Set major ticks every Monday
        ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MONDAY))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%e\n%b'))
        # Set minor ticks every day
        ax.xaxis.set_minor_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(mdates.DateFormatter('%e'))
    else:
        # Set major ticks every first day of the month
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%e\n%b'))
        # Set minor ticks every 15th of the month
        ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=15))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter('%e'))

    ax.tick_params(axis='x', which='minor', labelcolor='grey')
    ax.grid(axis='y')

    # half the labels because pos and neg create duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    half = int(len(handles) / 2)
    fig.legend(
        handles=handles[:half],
        labels=labels[:half],
        loc="outside right upper"
    )

    ax.axhline(0, color="grey", linewidth=0.5)

    if ylim is None:
        # ensure y-axis extent is symmetric around origin in steps of 100 units
        ylim = np.ceil(max(-neg.sum(axis=1).min(), pos.sum(axis=1).max()) / 100) * 100
    plt.ylim([-ylim, ylim])

    is_kt = any(s in ylabel.lower() for s in ["co2", "steel", "hvc"])
    unit = "kt/h" if is_kt else "GW"
    plt.ylabel(f"{ylabel} balance [{unit}]")

    if not long_time_frame:
        # plot frequency of snapshots on top x-axis as ticks
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(df.index)
        ax2.grid(False)
        ax2.tick_params(axis='x', length=2)
        ax2.xaxis.set_tick_params(labelbottom=False)
        ax2.set_xticklabels([])

    if resample is None:
        resample = "native"
    fn = f"ts-balance-{ylabel.replace(' ', '_')}-{resample}"
    plt.savefig(dir + "/" + fn + ".pdf")
    plt.savefig(dir + "/" + fn + ".png")



if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_balance_timeseries",
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

    months = pd.date_range(freq="M", **snakemake.config["snapshots"]).format(
        formatter=lambda x: x.strftime("%Y-%m")
    )

    balance = n.statistics.energy_balance(aggregate_time=False)

    n.carriers.color.update(snakemake.config["plotting"]["tech_colors"])
    colors = n.carriers.color.rename(n.carriers.nice_name)

    # wrap in function for multiprocessing
    def process_group(group, carriers, balance, months, colors):
        if not isinstance(carriers, list):
            carriers = [carriers]

        mask = balance.index.get_level_values("bus_carrier").isin(carriers)
        df = balance[mask].groupby("carrier").sum().div(1e3).T

        # native resolution for each month and carrier
        for month in months:
            plot_energy_balance_timeseries(
                df,
                time=month,
                ylabel=group,
                colors=colors,
                threshold=THRESHOLD,
                dir=snakemake.output[0]
            )

        # daily resolution for each carrier
        plot_energy_balance_timeseries(
            df,
            resample="D",
            ylabel=group,
            colors=colors,
            threshold=THRESHOLD,
            dir=snakemake.output[0]
        )

    args = [(group, carriers, balance, months, colors) for group, carriers in CARRIER_GROUPS.items()]
    with Pool(processes=snakemake.threads) as pool:
        pool.starmap(process_group, args)
