# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Plot balance time series.
"""

import logging
import os
from functools import partial
from multiprocessing import Pool

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging, get_snapshots, set_scenario_config
from tqdm import tqdm

logger = logging.getLogger(__name__)


def plot_stacked_area_steplike(
    ax: plt.Axes, df: pd.DataFrame, colors: dict | pd.Series = {}
):
    """Plot stacked area chart with step-like transitions."""
    if isinstance(colors, pd.Series):
        colors = colors.to_dict()

    df_cum = df.cumsum(axis=1)
    previous_series = np.zeros_like(df_cum.iloc[:, 0].values)

    for col in df_cum.columns:
        ax.fill_between(
            df_cum.index,
            previous_series,
            df_cum[col],
            step="pre",
            linewidth=0,
            color=colors.get(col, "grey"),
            label=col,
        )
        previous_series = df_cum[col].values


def setup_time_axis(ax: plt.Axes, timespan: pd.Timedelta):
    """Configure time axis formatting based on timespan."""
    long_time_frame = timespan > pd.Timedelta(weeks=5)

    if not long_time_frame:
        ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MONDAY))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%e\n%b"))
        ax.xaxis.set_minor_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%e"))
    else:
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%e\n%b"))
        ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=15))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%e"))

    ax.tick_params(axis="x", which="minor", labelcolor="grey")


def plot_energy_balance_timeseries(
    df: pd.DataFrame,
    time: pd.DatetimeIndex | None = None,
    ylim: float | None = None,
    resample: str | None = None,
    rename: dict = {},
    preferred_order: pd.Index | list = [],
    ylabel: str = "",
    colors: dict | pd.Series = {},
    max_threshold: float = 0.0,
    mean_threshold: float = 0.0,
    directory="",
):
    """Create energy balance time series plot with positive/negative stacked areas."""
    if time is not None:
        df = df.loc[time]

    # Handle small values and renaming
    techs_below_threshold = df.columns[
        (df.abs().max() < max_threshold) & (df.abs().mean() < mean_threshold)
    ].tolist()
    if techs_below_threshold:
        rename.update({tech: "other" for tech in techs_below_threshold})
        colors["other"] = "grey"

    if rename:
        df = df.T.groupby(df.columns.map(lambda a: rename.get(a, a))).sum().T

    # Upsample to hourly resolution to handle overlapping snapshots
    if resample is not None:
        df = df.resample("1h").ffill().resample(resample).mean()

    # Sort columns by variance
    order = (df / df.max()).var().sort_values().index
    if preferred_order:
        order = preferred_order.intersection(order).append(
            order.difference(preferred_order)
        )
    df = df.loc[:, order]

    # Split into positive and negative values
    pos = df.where(df > 0).fillna(0.0)
    neg = df.where(df < 0).fillna(0.0)

    # Create figure and plot
    fig, ax = plt.subplots(figsize=(10, 4.5), layout="constrained")
    plot_stacked_area_steplike(ax, pos, colors)
    plot_stacked_area_steplike(ax, neg, colors)

    # Set x and y limits
    plt.xlim((df.index[0], df.index[-1]))
    setup_time_axis(ax, df.index[-1] - df.index[0])

    # Configure y-axis and grid
    ax.grid(axis="y")
    ax.axhline(0, color="grey", linewidth=0.5)

    if ylim is None:
        # ensure y-axis extent is symmetric around origin in steps of 50 units
        ylim = np.ceil(max(-neg.sum(axis=1).min(), pos.sum(axis=1).max()) / 50) * 50
    plt.ylim([-ylim, ylim])

    # Set labels and legend
    unit = "kt/h" if "co2" in ylabel.lower() else "GW"
    plt.ylabel(f"{ylabel} balance [{unit}]")

    # half the labels because pos and neg create duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    half = int(len(handles) / 2)
    fig.legend(
        handles=handles[:half],
        labels=labels[:half],
        loc="outside right upper",
        fontsize=6,
    )

    # Save figures
    if resample is None:
        resample = f"native-{time if time is not None else 'default'}"
    fn = f"ts-balance-{ylabel.replace(' ', '_')}-{resample}.pdf"
    plt.savefig(f"{directory}/{fn}")
    plt.close()


def process_carrier(group_item, balance, months, colors, config, output_dir):
    """Process carrier data and create plots for specific carrier group."""

    group, carriers = group_item

    if not isinstance(carriers, list):
        carriers = [carriers]

    mask = balance.index.get_level_values("bus_carrier").isin(carriers)
    df = balance[mask].groupby("carrier").sum().div(1e3).T

    kwargs = dict(
        ylabel=group,
        colors=colors,
        max_threshold=config["max_threshold"],
        mean_threshold=config["mean_threshold"],
        directory=output_dir,
    )

    # daily resolution for each carrier
    if config["annual"]:
        plot_energy_balance_timeseries(
            df, resample=config["annual_resolution"], **kwargs
        )

    # native resolution for each month and carrier
    if config["monthly"]:
        for month in months:
            plot_energy_balance_timeseries(
                df, resample=config["monthly_resolution"], time=month, **kwargs
            )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_balance_timeseries",
            simpl="",
            clusters="10",
            opts="",
            sector_opts="",
            planning_horizons=2050,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    plt.style.use(["bmh", snakemake.input.rc])

    # Load network and prepare data
    n = pypsa.Network(snakemake.input.network)
    config = snakemake.params.plotting["balance_timeseries"]
    output_dir = snakemake.output[0]
    os.makedirs(output_dir, exist_ok=True)

    # Get month ranges for plotting
    sns = snakemake.params.snapshots
    drop_leap_day = snakemake.params.drop_leap_day
    months = get_snapshots(sns, drop_leap_day, freq="ME").map(
        lambda x: x.strftime("%Y-%m")
    )

    # Calculate energy balance
    balance = n.statistics.energy_balance(aggregate_time=False, nice_names=False)

    # Get colors for carriers
    n.carriers.update({"color": snakemake.params.plotting["tech_colors"]})
    colors = n.carriers.color.copy().replace("", "grey")

    # Setup carrier groups for plotting
    groups = config["carrier_groups"]
    groups.update({c: [c] for c in config["carriers"]})

    missing_bus_carriers = set(config["carriers"]).difference(n.buses.carrier.unique())
    if missing_bus_carriers:
        logger.warning(f"Skipping missing bus carriers: {missing_bus_carriers}")
    groups = {k: v for k, v in groups.items() if k not in missing_bus_carriers}

    # Process each carrier group in parallel
    threads = snakemake.threads
    tqdm_kwargs = dict(
        ascii=False,
        unit=" carrier",
        total=len(groups),
        desc="Plotting carrier balance time series",
    )
    func = partial(
        process_carrier,
        balance=balance,
        months=months,
        colors=colors,
        config=config,
        output_dir=output_dir,
    )
    with Pool(processes=min(threads, len(groups))) as pool:
        list(tqdm(pool.imap(func, groups.items()), **tqdm_kwargs))
