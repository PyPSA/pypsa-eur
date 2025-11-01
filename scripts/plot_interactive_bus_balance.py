# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create interactive bus energy balance time series plots.

This script generates interactive HTML plots showing energy balance time series
for buses in the final network. It calculates and visualizes the contribution of
different carriers (generation, storage, loads, links) to the energy balance
at specified buses, creating stacked area charts with positive/negative values.

The plots show generation as positive values and consumption as negative values to provide an energy balance at each bus over time.
The scripts does not use `n.statistics.energy_balance` but calculates the balance directly from the time series data of generators, storage units, loads, and links connected to each bus.


Relevant Settings
-----------------

.. code:: yaml

    plotting:
        tech_colors: # Color mapping for different technologies/carriers
        balance_timeseries:
            bus_name_pattern: # Pattern to filter buses (e.g., 'DE*' for German buses)

Inputs
------
- `resources/<run_name>/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc`: Solved PyPSA network
- `config/plotting/rc.mplstyle`: Matplotlib style configuration

Outputs
-------
- `results/<run_name>/plots/balance_timeseries/`: Directory containing HTML files with interactive plots
  - `ts-balance-{bus_name}-native-{time}.html`: Interactive time series plot for each bus

Notes
-----
Uses Plotly for interactive visualization. Supports filtering buses by pattern
(shell-style wildcards). Processes multiple buses in parallel for efficiency.
"""

import fnmatch
import logging
import os
from functools import partial
from multiprocessing import Pool

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pypsa
from plotly.subplots import make_subplots
from tqdm import tqdm

from scripts._helpers import configure_logging, get_snapshots, set_scenario_config

logger = logging.getLogger(__name__)


# Add a function to convert matplotlib color codes to Plotly compatible colors
def convert_color_for_plotly(color: str) -> str:
    """Convert matplotlib color codes to Plotly-compatible color formats."""
    # Common matplotlib single letter color codes
    mpl_to_plotly = {
        "k": "black",
        "b": "blue",
        "r": "red",
        "g": "green",
        "y": "yellow",
        "c": "cyan",
        "m": "magenta",
        "w": "white",
    }

    if color in mpl_to_plotly:
        return mpl_to_plotly[color]
    return color


def get_bus_balance(n: pypsa.Network, bus_name: str) -> pd.DataFrame:
    """
    Calculate the energy balance for a specific bus.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network
    bus_name : str
        Name of the bus

    Returns
    -------
    pd.DataFrame
        DataFrame with carrier dispatch data
    """
    carriers = {}

    # Generation from generators
    generators = n.generators.index[n.generators.bus == bus_name]
    for gen in generators:
        carrier = n.generators.carrier[gen]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        carriers[carrier] += n.generators_t.p[gen]

    # Generation from storage units
    storage_units = n.storage_units.index[n.storage_units.bus == bus_name]
    for su in storage_units:
        carrier = n.storage_units.carrier[su]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        carriers[carrier] += n.storage_units_t.p[su]

    # Handle links connected to the bus via any bus connection (bus0, bus1, bus2, bus3)

    # Links with bus_name as bus0 (input power p0)
    for link in n.links.index[n.links.bus0 == bus_name]:
        carrier = n.links.carrier[link]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        carriers[carrier] += -n.links_t.p0[
            link
        ]  # negative sign for consumption at bus0

    # Links with bus_name as bus1 (output power p1)
    for link in n.links.index[n.links.bus1 == bus_name]:
        carrier = n.links.carrier[link]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        carriers[carrier] += -n.links_t.p1[
            link
        ]  # negative sign because p1 is outflow from link perspective

    # Links with bus_name as bus2 (output power p2)
    bus2_links = n.links.index[n.links.bus2 == bus_name]
    for link in bus2_links:
        carrier = n.links.carrier[link]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        if link in n.links_t.p2.columns:
            carriers[carrier] += -n.links_t.p2[
                link
            ]  # negative sign because p2 is outflow from link perspective

    # Links with bus_name as bus3 (output power p3)
    bus3_links = n.links.index[n.links.bus3 == bus_name]
    for link in bus3_links:
        carrier = n.links.carrier[link]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        if link in n.links_t.p3.columns:
            carriers[carrier] += -n.links_t.p3[
                link
            ]  # negative sign because p3 is outflow from link perspective

    # Conventional load
    loads = n.loads.index[n.loads.bus == bus_name]
    for load in loads:
        carrier = n.loads.carrier[load]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        carriers[carrier] += -n.loads_t.p[load]  # negative sign for consumption

    # Storage charging (consumption)
    stores = n.stores.index[n.stores.bus == bus_name]
    for store in stores:
        carrier = n.stores.carrier[store]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        carriers[carrier] += -n.stores_t.p[store]  # negative sign for consumption

    # Create a dataframe with all carriers
    result = pd.DataFrame(carriers, index=n.snapshots)

    # Remove columns (carriers) that have no contribution (all zeros or all NaN)
    # First fill NaN with 0 to properly check for zero columns
    result_filled = result.fillna(0)

    # Keep only columns that have at least one non-zero value
    non_zero_columns = result_filled.columns[result_filled.abs().sum() > 1e-2]
    result = result[non_zero_columns]

    # Drop rows (time periods) where all carriers are zero/NaN
    result = result.dropna(axis=0, how="all")

    return result


def prepare_all_buses_data(
    n: pypsa.Network, bus_name_pattern: str | None = None
) -> dict[str, pd.DataFrame]:
    """
    Prepare data for all buses in the network, optionally filtering by a pattern.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network
    bus_name_pattern : str, optional
        Shell-style pattern to filter bus names (e.g., 'DE*')
        If set to "NONE_BY_DEFAULT", no buses will be selected.

    Returns
    -------
    dict
        Dictionary with bus names as keys and carrier dispatch DataFrames as values
    """
    buses_data = {}

    # Check for the special "NONE_BY_DEFAULT" pattern
    if bus_name_pattern == "NONE_BY_DEFAULT":
        logger.info(
            "No buses selected for plotting due to NONE_BY_DEFAULT pattern. Set a different pattern to plot specific buses."
        )
        return buses_data

    # Get all buses from the network
    all_buses = list(n.buses.index)
    if bus_name_pattern:
        all_buses = [bus for bus in all_buses if fnmatch.fnmatch(bus, bus_name_pattern)]
    logger.info(
        f"Found {len(all_buses)} buses matching pattern '{bus_name_pattern}'"
        if bus_name_pattern
        else f"Found {len(all_buses)} total buses in network"
    )

    # Process all buses
    for bus in all_buses:
        try:
            # Get the bus balance
            bus_balance = get_bus_balance(n, bus)

            # Only add buses that have non-empty carrier data
            if len([col for col in bus_balance.columns if col != "time"]) > 0:
                buses_data[bus] = bus_balance
        except Exception as e:
            raise RuntimeError(f"Error processing bus {bus}: {e}")
    logger.info(f"Successfully processed {len(buses_data)} buses with carrier data")

    return buses_data


def plot_stacked_area_steplike(
    ax: plt.Axes, df: pd.DataFrame, colors: dict[str, str] | pd.Series = {}
) -> None:
    """
    Plot stacked area chart with step-like transitions.

    Parameters
    ----------
    ax : matplotlib.pyplot.Axes
        Matplotlib axes object to plot on.
    df : pd.DataFrame
        DataFrame with time series data to plot.
    colors : dict[str, str] | pd.Series, optional
        Color mapping for carriers.
    """
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


def setup_time_axis(ax: plt.Axes, timespan: pd.Timedelta) -> None:
    """
    Configure time axis formatting based on timespan.

    Parameters
    ----------
    ax : matplotlib.pyplot.Axes
        Matplotlib axes object to configure.
    timespan : pd.Timedelta
        Time range of the data for appropriate formatting.
    """
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
    rename: dict[str, str] = {},
    ylabel: str = "",
    colors: dict[str, str] | pd.Series = {},
    directory: str = "",
) -> None:
    """
    Create interactive energy balance time series plot with positive/negative stacked areas.

    Parameters
    ----------
    df : pd.DataFrame
        Energy balance data with carriers as columns.
    time : pd.DatetimeIndex, optional
        Time period to plot. If None, plots full time series.
    ylim : float, optional
        Y-axis limits. If None, calculated automatically.
    resample : str, optional
        Resampling frequency (e.g., 'D' for daily).
    rename : dict[str, str], optional
        Mapping to rename carriers for display.
    ylabel : str, optional
        Label for y-axis and filename.
    colors : dict[str, str] | pd.Series, optional
        Color mapping for carriers.
    directory : str, optional
        Output directory for HTML file.
    """
    if time is not None:
        df = df.loc[time]

    if rename:
        df = df.T.groupby(df.columns.map(lambda a: rename.get(a, a))).sum().T

    # Upsample to hourly resolution to handle overlapping snapshots
    if resample is not None:
        df = df.resample("1h").ffill().resample(resample).mean()

    # Sort columns alphabetically
    order = df.columns.sort_values()

    df = df.loc[:, order]

    # Split into positive and negative values
    pos = df.where(df > 0).fillna(0.0)
    neg = df.where(df < 0).fillna(0.0)

    # Convert matplotlib color codes to Plotly compatible formats
    plotly_colors = {
        carrier: convert_color_for_plotly(color) for carrier, color in colors.items()
    }

    # Create plotly figure
    fig = make_subplots(specs=[[{"secondary_y": False}]])

    # Keep track of carriers already added to avoid duplicates
    carriers_added = set()

    # Plot positive values
    pos_reset = pos.reset_index()
    date_col = pos_reset.columns[0]  # Get the actual name of the reset index column
    pos_df_melt = pos_reset.melt(
        id_vars=date_col, var_name="carrier", value_name="value"
    )

    for carrier in pos.columns:
        carrier_df = pos_df_melt[pos_df_melt["carrier"] == carrier]
        # Only add if carrier has actual positive values and hasn't been added yet
        if not carrier_df.empty and carrier_df["value"].max() > 0:
            fig.add_trace(
                go.Scatter(
                    x=carrier_df[date_col],
                    y=carrier_df["value"],
                    mode="lines",
                    name=carrier,
                    line=dict(width=0, color=plotly_colors.get(carrier, "grey")),
                    stackgroup="positive",
                    fill="tonexty",
                    hovertemplate=f"{carrier}: %{{y:.2f}}<extra></extra>",
                )
            )
            carriers_added.add(carrier)

    # Plot negative values
    neg_reset = neg.reset_index()
    date_col_neg = neg_reset.columns[0]
    neg_df_melt = neg_reset.melt(
        id_vars=date_col_neg, var_name="carrier", value_name="value"
    )

    for carrier in neg.columns:
        carrier_df = neg_df_melt[neg_df_melt["carrier"] == carrier]
        # Only add if carrier has actual negative values and hasn't been added yet
        if (
            not carrier_df.empty
            and carrier_df["value"].min() < 0
            and carrier not in carriers_added
        ):
            fig.add_trace(
                go.Scatter(
                    x=carrier_df[date_col_neg],
                    y=carrier_df["value"],
                    mode="lines",
                    name=carrier,
                    line=dict(width=0, color=plotly_colors.get(carrier, "grey")),
                    stackgroup="negative",
                    fill="tonexty",
                    hovertemplate=f"{carrier}: %{{y:.2f}}<extra></extra>",
                )
            )
            carriers_added.add(carrier)

    # Set y-axis limits
    if ylim is None:
        # ensure y-axis extent is symmetric around origin in steps of 50 units
        ylim = np.ceil(max(-neg.sum(axis=1).min(), pos.sum(axis=1).max()) / 50) * 50

    # Set layout
    unit = "kt/h" if "co2" in ylabel.lower() else "MW"
    fig.update_layout(
        title="",
        yaxis_title=f"{ylabel} balance [{unit}]",
        hovermode="x unified",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        yaxis=dict(range=[-ylim, ylim]),
        plot_bgcolor="white",
    )

    # Add a horizontal line at y=0
    fig.add_shape(
        type="line",
        x0=df.index[0],
        y0=0,
        x1=df.index[-1],
        y1=0,
        line=dict(color="grey", width=1),
    )

    # Add grid lines
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor="lightgrey")

    # Save as interactive HTML
    if resample is None:
        resample = f"native-{time if time is not None else 'default'}"
    fn = f"ts-balance-{ylabel.replace(' ', '_')}-{resample}.html"
    fig.write_html(f"{directory}/{fn}")


def process_carrier(
    bus_name: str,
    bus_data: dict[str, pd.DataFrame],
    colors: dict[str, str] | pd.Series,
    output_dir: str,
) -> None:
    """
    Process carrier data and create plots for specific bus.

    Parameters
    ----------
    bus_name : str
        Name of the bus to process.
    bus_data : dict[str, pd.DataFrame]
        Dictionary mapping bus names to energy balance DataFrames.
    colors : dict[str, str] | pd.Series
        Color mapping for carriers.
    output_dir : str
        Directory to save output files.
    """

    df = bus_data[bus_name]

    kwargs = dict(
        ylabel=bus_name,
        colors=colors,
        directory=output_dir,
    )

    plot_energy_balance_timeseries(df, resample=None, **kwargs)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

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
    # config = snakemake.params.plotting["balance_timeseries"]
    output_dir = snakemake.output[0]
    os.makedirs(output_dir, exist_ok=True)

    # Get bus name pattern from Snakemake params (preferred) or config (fallback)
    bus_name_pattern = snakemake.params.get(
        "bus_name_pattern", snakemake.params.bus_name_pattern
    )

    # Get month ranges for plotting
    sns = snakemake.params.snapshots
    drop_leap_day = snakemake.params.drop_leap_day
    months = get_snapshots(sns, drop_leap_day, freq="ME").map(
        lambda x: x.strftime("%Y-%m")
    )

    # Calculate energy balance
    bus_data = prepare_all_buses_data(n, bus_name_pattern=bus_name_pattern)
    if len(bus_data) == 0:
        logger.warning(
            f"No buses found matching the pattern {bus_name_pattern}. Exiting."
        )

    else:
        logger.info(f"Found {len(bus_data)} buses with carrier data for plotting.")
        # Get colors for carriers
        n.carriers.update({"color": snakemake.params.plotting["tech_colors"]})
        colors = n.carriers.color.copy().replace("", "grey")
        # Process each carrier group in partial
        threads = snakemake.threads
        tqdm_kwargs = dict(
            ascii=False,
            unit=" carrier",
            total=len(bus_data.keys()),
            desc="Plotting carrier balance time series",
        )
        func = partial(
            process_carrier,
            bus_data=bus_data,
            colors=colors,
            output_dir=output_dir,
        )
        with Pool(processes=min(threads, len(bus_data.keys()))) as pool:
            list(tqdm(pool.imap(func, bus_data.keys()), **tqdm_kwargs))
