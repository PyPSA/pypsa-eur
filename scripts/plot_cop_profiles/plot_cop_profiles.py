# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Plot coefficient of performance (COP) profiles for heat pumps in different regions.
Generates interactive HTML plots showing COP profiles for central heating.
"""

import logging

import numpy as np
import pandas as pd
import xarray as xr
from _helpers import configure_logging
from bokeh.io import output_file, save
from bokeh.layouts import column, row
from bokeh.models import (
    ColumnDataSource,
    CustomJS,
    Div,
    HoverTool,
    Select,
)
from bokeh.palettes import Category10
from bokeh.plotting import figure
from bokeh.transform import dodge

from scripts.definitions.heat_system_type import HeatSystemType

logger = logging.getLogger(__name__)


def prepare_cop_data(cop_profiles, heat_system_type: HeatSystemType, region_dim="name"):
    """
    Prepare COP data for plotting.
    Handles 4-dimensional data (time, name, heat_source, heat_system)
    and extracts data for the specified heat system type.

    Parameters
    ----------
    cop_profiles : xarray.Dataset or xarray.DataArray
        COP profiles dataset or data array
    heat_system : str
        Heat system type to filter for (e.g., "urban central")

    Returns
    -------
    pandas.DataFrame
        Prepared COP data for plotting
    pandas.DataFrame
        Monthly average COP data for bar chart
    list
        List of region names (nodes)
    list
        List of heat source names
    """
    # If we have a Dataset, convert to DataArray
    if isinstance(cop_profiles, xr.Dataset):
        # Assuming there's one main variable in the dataset
        var_name = list(cop_profiles.data_vars)[0]
        logger.info(f"Converting Dataset to DataArray using variable: {var_name}")
        cop_profiles = cop_profiles[var_name]

    # Filter to the specified heat system type
    try:
        # Check if the specified heat system exists in the heat_system dimension
        cop_data = cop_profiles.sel(heat_system=heat_system_type)
        logger.info(f"Selected '{heat_system_type}' heat system")
    except Exception as e:
        raise RuntimeError(f"Error selecting heat system: {e}")

    # Get the name of the region dimension
    # Capture heat source names before pivoting
    try:
        heat_sources = [val for val in cop_data.coords["heat_source"].values]
        # logger.info(f"Heat sources: {heat_sources}")
    except Exception as e:
        raise RuntimeError(f"Error retrieving heat sources: {e}")

    # Convert to pandas for plotting
    # We need to reshape data to have heat sources as columns
    try:
        # Start with a multi-index DataFrame
        df = cop_data.to_dataframe().reset_index()
        # Check for NaN values before pivoting
        cop_col = var_name

        # Pivot the table to have heat sources as columns
        pivot_df = df.pivot_table(
            index=["time", region_dim],
            columns="heat_source",
            values=cop_col,
            # Do not include NaN values in the average
            aggfunc="mean",
            fill_value=None,  # Keep NaN values as NaN
        ).reset_index()

    except Exception as e:
        logger.error(f"Unexpected data structure. Error processing COP data: {e}")

    # Calculate monthly averages for bar chart, ignoring NaN values
    if "time" in pivot_df.columns:
        # Add month column
        pivot_df["month"] = pd.to_datetime(pivot_df["time"]).dt.month

        # Group by month and region, calculate mean of each heat source
        monthly_avg = (
            pivot_df.groupby(["month", region_dim])
            .mean(numeric_only=True)
            .reset_index()
        )
        # Add month column and month names in one go
        pivot_df["month"] = pd.to_datetime(pivot_df["time"]).dt.month
        pivot_df["month_name"] = pd.to_datetime(pivot_df["time"]).dt.strftime("%b")

        # Group by month and region, calculate mean of each heat source
        monthly_avg = (
            pivot_df.groupby(["month", region_dim, "month_name"])
            .mean(numeric_only=True)
            .reset_index()
            # Sort by month
            .sort_values("month")
        )

    else:
        # Create empty DataFrame if time data is not available
        monthly_avg = pd.DataFrame(
            columns=["month", region_dim, "month_name"] + heat_sources
        )

    # Get unique region names
    regions = sorted(pivot_df[region_dim].unique())
    logger.info(f"Found {len(regions)} regions: {regions[:5]}...")

    # Add heat system information to the dataframes for reference
    pivot_df["heat_system"] = heat_system_type
    monthly_avg["heat_system"] = heat_system_type

    return pivot_df, monthly_avg, regions, heat_sources


def create_interactive_cop_plot(
    cop_df, monthly_avg_df, regions, heat_sources, region_dim: str = "name"
):
    """
    Create an interactive Bokeh plot for COP profiles with a monthly average bar chart.

    Parameters
    ----------
    cop_df : pandas.DataFrame
        DataFrame with COP data
    monthly_avg_df : pandas.DataFrame
        DataFrame with monthly average COP data
    regions : list
        List of region names
    heat_sources : list
        List of heat source names

    Returns
    -------
    bokeh.layouts.column
        Bokeh layout with interactive plots
    """
    # Use the first region as default
    default_region = regions[0] if regions else None

    if default_region is None:
        logger.error("No regions found in COP data")
        return None

    # Filter for the default region
    region_data = cop_df[cop_df[region_dim] == default_region]
    monthly_region_data = monthly_avg_df[monthly_avg_df[region_dim] == default_region]

    # Create a color palette with enough colors for all heat sources
    colors = Category10[10] * (len(heat_sources) // 10 + 1)
    colors = colors[: len(heat_sources)]

    # Create ColumnDataSources for both plots
    timeseries_source = ColumnDataSource(region_data)
    monthly_source = ColumnDataSource(monthly_region_data)

    # Create the time series figure
    p_timeseries = figure(
        width=900,
        height=500,
        title=f"Heat Pump COP Profiles for {default_region}",
        x_axis_type="datetime",
        toolbar_location="above",
        tools=[
            "pan",
            "wheel_zoom",
            "box_zoom",
            "reset",
            "save",
        ],  # Removed hover tool from here
    )

    # Add lines for each heat source to time series plot
    for i, heat_source in enumerate(heat_sources):
        if heat_source in region_data.columns:
            line = p_timeseries.line(
                x="time",
                y=heat_source,
                source=timeseries_source,
                line_width=2,
                color=colors[i],
                legend_label=heat_source,
                name=heat_source,
            )

            # Add hover tool only to the first heat source
            if i == 0:
                # Create tooltip data for time series
                timeseries_tooltips = [("Date", "@time{%F}")]
                for hs in heat_sources:
                    timeseries_tooltips.append((hs, f"@{{{hs}}}{{0.00}}"))

                # Create hover tool for time series
                timeseries_hover = HoverTool(
                    renderers=[line],  # Only attach to this line
                    tooltips=timeseries_tooltips,
                    formatters={"@time": "datetime"},
                    mode="vline",  # Still show all values at the hovered time
                    point_policy="follow_mouse",
                )
                p_timeseries.add_tools(timeseries_hover)
        else:
            logger.warning(
                f"Heat source '{heat_source}' not found in timeseries data for {default_region}"
            )

    # Style the time series plot
    p_timeseries.title.text_font_size = "14pt"
    p_timeseries.xaxis.axis_label = "Time"
    p_timeseries.yaxis.axis_label = "Coefficient of Performance (COP)"
    p_timeseries.legend.location = "top_left"
    p_timeseries.legend.click_policy = "hide"
    p_timeseries.grid.grid_line_alpha = 0.3

    # Create monthly average bar chart
    p_monthly = figure(
        width=500,
        height=500,
        title=f"Monthly Average COP for {default_region}",
        x_range=monthly_region_data["month_name"].tolist(),
        toolbar_location="above",
        tools=["pan", "wheel_zoom", "box_zoom", "reset", "save"],
    )

    # Create tooltip for monthly bar chart
    monthly_tooltips = [("Month", "@month_name")]
    for heat_source in heat_sources:
        monthly_tooltips.append((heat_source, f"@{{{heat_source}}}{{0.00}}"))

    p_monthly.add_tools(HoverTool(tooltips=monthly_tooltips))

    # Add grouped bars for each heat source
    bar_width = 0.8 / len(heat_sources)

    for i, heat_source in enumerate(heat_sources):
        if heat_source in monthly_source.data:
            # Calculate offset for grouped bars
            offset = (i - len(heat_sources) / 2 + 0.5) * bar_width

            p_monthly.vbar(
                x=dodge("month_name", offset, range=p_monthly.x_range),
                top=heat_source,
                width=bar_width,
                source=monthly_source,
                color=colors[i],
                legend_label=heat_source,
                name=heat_source,
            )
        else:
            logger.warning(
                f"Heat source '{heat_source}' not found in monthly data for {default_region}"
            )

    # Style the monthly bar chart
    p_monthly.title.text_font_size = "14pt"
    p_monthly.xaxis.axis_label = "Month"
    p_monthly.yaxis.axis_label = "Average COP"
    p_monthly.xgrid.grid_line_color = None
    p_monthly.legend.location = "top_left"
    p_monthly.legend.click_policy = "hide"

    # Rotate x-axis labels for better readability
    p_monthly.xaxis.major_label_orientation = np.pi / 4

    # Create region selector dropdown
    region_select = Select(
        title="Select Region:",
        value=default_region,
        options=regions,
        width=200,
    )

    # Create callback to update both plots when region is changed
    callback = CustomJS(
        args=dict(
            timeseries_source=timeseries_source,
            monthly_source=monthly_source,
            region_select=region_select,
            p_timeseries=p_timeseries,
            p_monthly=p_monthly,
            all_timeseries_data=ColumnDataSource(cop_df),
            all_monthly_data=ColumnDataSource(monthly_avg_df),
            region_dim=region_dim,
        ),
        code="""
        const region = region_select.value;
        const regionColumn = region_dim;  // Use the provided region dimension name

        // Update timeseries plot
        const tsData = all_timeseries_data.data;
        const tsOut = {};
        // Filter data for selected region
        const tsIndices = [];
        for (let i = 0; i < tsData[regionColumn].length; i++) {
            if (tsData[regionColumn][i] === region) {
                tsIndices.push(i);
            }
        }

        // Get all column names
        const tsColumns = Object.keys(tsData);

        // For each column, filter the data
        for (let col of tsColumns) {
            tsOut[col] = tsIndices.map(i => tsData[col][i]);
        }

        // Update the source data
        timeseries_source.data = tsOut;

        // Update monthly plot
        const monthlyData = all_monthly_data.data;
        const monthlyOut = {};

        // Filter monthly data for selected region
        const monthlyIndices = [];
        for (let i = 0; i < monthlyData[regionColumn].length; i++) {
            if (monthlyData[regionColumn][i] === region) {
                monthlyIndices.push(i);
            }
        }

        // Get all column names
        const monthlyColumns = Object.keys(monthlyData);

        // For each column, filter the data
        for (let col of monthlyColumns) {
            monthlyOut[col] = monthlyIndices.map(i => monthlyData[col][i]);
        }

        // Update the source data
        monthly_source.data = monthlyOut;

        // Update the plot titles
        p_timeseries.title.text = `Heat Pump COP Profiles for ${region}`;
        p_monthly.title.text = `Monthly Average COP for ${region}`;

        // Update x-range for bar chart to ensure proper ordering
        if (monthlyOut.month_name && monthlyOut.month_name.length > 0) {
            p_monthly.x_range.factors = monthlyOut.month_name;
        }
        """,
    )

    # Connect callback to dropdown
    region_select.js_on_change("value", callback)

    # Create info text
    info_div = Div(
        text="""
        <p> Warning: Experimental - check values!</p>
        <p>This plot shows Coefficient of Performance (COP) profiles for heat pumps in different regions.</p>
        <p>COP represents the efficiency of heat pumps - higher values mean better performance.</p>
        <p>Select a region from the dropdown to see its specific COP profiles for different heat sources.</p>
        <p>Click on legend items to hide/show specific heat sources.</p>
        <p>The bar chart shows monthly averages to help identify seasonal patterns.</p>
        """,
        width=400,
        styles={"font-size": "12px"},
    )

    # Create layout
    controls = column(region_select, info_div, width=200)
    plots = row(p_timeseries, p_monthly)
    layout = column(row(controls, plots))

    return layout


if __name__ == "__main__":
    """Generate interactive COP profile plots for all heat systems in a single file."""
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_cop_profiles",
            clusters=48,
            planning_horizons=2030,
        )

    configure_logging(snakemake)

    # Load COP profiles
    cop_profiles = xr.open_dataset(snakemake.input.cop_profiles)
    logger.info(f"Loaded COP profiles with dimensions: {cop_profiles.dims}")

    # Process each heat system and collect layouts
    all_layouts = []

    for heat_system_type in HeatSystemType:
        logger.info(f"Processing heat system: {heat_system_type.value}")

        # Prepare data for plotting
        cop_df, monthly_avg_df, regions, heat_sources = prepare_cop_data(
            cop_profiles, heat_system_type.value
        )

        if len(regions) == 0:
            logger.error(
                f"No regions found in COP data for heat system {heat_system_type.value}"
            )
            continue

        if len(heat_sources) == 0:
            logger.error(
                f"No heat sources found in COP data for heat system {heat_system_type.value}"
            )
            continue

        # Create layout
        layout = create_interactive_cop_plot(
            cop_df, monthly_avg_df, regions, heat_sources
        )

        # Create a header for this heat system
        header = Div(
            text=f"<h2>Heat Pump COP Profiles - {heat_system_type.value}</h2>",
            width=800,
            styles={"font-size": "16px", "font-weight": "bold", "margin-top": "20px"},
        )

        # Add header and layout to the collection
        all_layouts.append(header)
        all_layouts.append(layout)

    # Combine all layouts into a single column
    if all_layouts:
        combined_layout = column(all_layouts)

        # Save the combined plot to snakemake.output[0]
        output_file(snakemake.output[0], title="COP Profiles")
        save(combined_layout)
        logger.info(f"Interactive COP profile plots saved to {snakemake.output[0]}")
    else:
        logger.error("No heat system data could be processed. No output generated.")
