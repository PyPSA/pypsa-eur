# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Plot carrier dispatch time series for different buses with interactive dropdown selection.
"""

import logging
import os
from pathlib import Path

import pandas as pd
import pypsa
from bokeh.io import output_file, save

# Import layout components with aliases to avoid naming conflicts
from bokeh.layouts import column as bokeh_column
from bokeh.layouts import row as bokeh_row
from bokeh.models import (
    ColumnDataSource,
    CustomJS,
    Div,
    HoverTool,
    Select,
)
from bokeh.palettes import Category20
from bokeh.plotting import figure

logger = logging.getLogger(__name__)


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

    # Generation from links (only taking links with bus_name as bus1)
    for link in n.links.index[n.links.bus1 == bus_name]:
        carrier = n.links.carrier[link]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        carriers[carrier] += -n.links_t.p1[
            link
        ]  # negative sign because bus1 is outflow

    # Load from links (only taking links with bus_name as bus0)
    for link in n.links.index[n.links.bus0 == bus_name]:
        carrier = n.links.carrier[link]
        p0 = n.links_t.p0[link]
        if f"{carrier} load" not in carriers:
            carriers[f"{carrier} load"] = 0.0
        carriers[f"{carrier} load"] += -p0  # negative sign for consumption

    # Conventional load
    loads = n.loads.index[n.loads.bus == bus_name]
    for load in loads:
        carrier = n.loads.carrier[load]
        if f"{carrier} load" not in carriers:
            carriers[f"{carrier} load"] = 0.0
        carriers[f"{carrier} load"] += -n.loads_t.p[
            load
        ]  # negative sign for consumption

    # Storage charging (consumption)
    stores = n.stores.index[n.stores.bus == bus_name]
    for store in stores:
        carrier = n.stores.carrier[store]
        if f"{carrier} charge" not in carriers:
            carriers[f"{carrier} charge"] = 0.0
        carriers[f"{carrier} charge"] += -n.stores_t.p[
            store
        ]  # negative sign for consumption

    # Create a dataframe with all carriers
    result = pd.DataFrame(carriers, index=n.snapshots)

    # Add a datetime column for plotting
    result["time"] = pd.to_datetime(result.index)

    return result


def prepare_all_buses_data(n: pypsa.Network) -> dict:
    """
    Prepare data for all buses in the network.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network

    Returns
    -------
    dict
        Dictionary with bus names as keys and carrier dispatch DataFrames as values
    """
    buses_data = {}

    # Get all buses from the network
    all_buses = list(n.buses.index)
    logger.info(f"Found {len(all_buses)} total buses in network")

    # Process all buses
    for bus in all_buses:
        try:
            # Get the bus balance
            bus_balance = get_bus_balance(n, bus)

            # Only add buses that have non-empty carrier data
            if len([col for col in bus_balance.columns if col != "time"]) > 0:
                buses_data[bus] = bus_balance
        except Exception as e:
            raise RuntimeError(
                f"Error processing bus {bus}: {e}"
            )
    logger.info(f"Successfully processed {len(buses_data)} buses with carrier data")

    return buses_data


def create_dispatch_plot(bus_data, default_bus_name):
    """
    Create an interactive Bokeh plot for carrier dispatch.

    Parameters
    ----------
    bus_data : dict
        Dictionary with bus names as keys and carrier dispatch DataFrames as values
    default_bus_name : str
        Name of the default bus to display

    Returns
    -------
    bokeh.layouts.column
        Bokeh layout with interactive plot
    """
    # Get list of all buses
    buses = sorted(bus_data.keys())

    # If the requested bus is not in the data, use the first available bus
    if default_bus_name not in buses:
        if buses:
            default_bus_name = buses[0]
        else:
            logger.error("No buses found with dispatch data")
            return None

    # Create a ColumnDataSource for each bus
    sources = {}
    for bus, data in bus_data.items():
        sources[bus] = ColumnDataSource(data)

    # Create a dictionary to track carriers and colors for each bus
    bus_carriers = {}

    for bus, data in bus_data.items():
        current_carriers = [col for col in data.columns if col != "time"]
        # Determine positive (generation) and negative (consumption) carriers
        pos_carriers = []
        neg_carriers = []
        for carrier in current_carriers:
            if data[carrier].mean() >= 0:
                pos_carriers.append(carrier)
            else:
                neg_carriers.append(carrier)

        # Store them for this bus
        bus_carriers[bus] = {
            "all": current_carriers,
            "positive": pos_carriers,
            "negative": neg_carriers,
        }

    # Create a consistent color mapping across all buses
    # Collect all unique carriers
    all_unique_carriers = set()
    for bus_data in bus_carriers.values():
        all_unique_carriers.update(bus_data["all"])

    # Create a palette with enough colors
    color_palette = Category20[20] * ((len(all_unique_carriers) // 20) + 1)
    color_palette = color_palette[: len(all_unique_carriers)]

    # Create a mapping from carrier to color
    carrier_colors = dict(zip(sorted(all_unique_carriers), color_palette))

    # DIFFERENT APPROACH: Instead of creating one plot with multiple renderers,
    # we'll create a separate plot for each bus (initially hidden) and switch between them.
    # This ensures each plot has its own properly configured legend.

    plots = {}
    for bus, bus_source in sources.items():
        # Create a new figure for this bus
        p = figure(
            width=900,
            height=500,
            title=f"Carrier Dispatch for {bus}",
            x_axis_type="datetime",
            toolbar_location="above",
            tools=["pan", "wheel_zoom", "box_zoom", "reset", "save"],
        )

        # Get carriers for this bus
        pos_carriers = bus_carriers[bus]["positive"]
        neg_carriers = bus_carriers[bus]["negative"]

        # Create stacked areas for positive (generation) carriers
        if pos_carriers:
            pos_colors = [carrier_colors[carrier] for carrier in pos_carriers]
            p.varea_stack(
                pos_carriers,
                x="time",
                color=pos_colors,
                legend_label=pos_carriers,
                source=bus_source,
            )

        # Create stacked areas for negative (consumption) carriers
        if neg_carriers:
            neg_colors = [carrier_colors[carrier] for carrier in neg_carriers]
            p.varea_stack(
                neg_carriers,
                x="time",
                color=neg_colors,
                legend_label=neg_carriers,
                source=bus_source,
            )

        # Create hover tooltips dynamically based on the bus carriers
        tooltips = [("Time", "@time{%F %H:%M}")]
        for carrier in bus_carriers[bus]["all"]:
            tooltips.append((carrier, f"@{{{carrier}}}{{0.00}} MW"))

        hover = HoverTool(
            tooltips=tooltips,
            formatters={"@time": "datetime"},
            mode="vline",
        )
        p.add_tools(hover)

        # Style the plot
        p.title.text_font_size = "14pt"
        p.xaxis.axis_label = "Time"
        p.yaxis.axis_label = "Power [MW]"
        p.axis.axis_label_text_font_style = "normal"
        p.grid.grid_line_alpha = 0.3
        p.legend.location = "right"
        p.legend.click_policy = "hide"

        # Store the plot with its bus name
        plots[bus] = p

        # Initially hide all plots except for the default bus
        if bus != default_bus_name:
            p.visible = False

    # Create a container to hold all plots
    plot_container = bokeh_column(
        children=[plot for plot in plots.values()], sizing_mode="stretch_width"
    )

    # Create bus selector dropdown
    bus_select = Select(
        title="Select Bus:",
        value=default_bus_name,
        options=buses,
        width=400,
    )

    # Simple callback to show/hide plots based on selection
    callback = CustomJS(
        args=dict(plots=plots, bus_select=bus_select),
        code="""
        // Get the selected bus name
        const selected_bus = bus_select.value;
        
        // Show/hide plots based on selection
        Object.keys(plots).forEach(bus => {
            plots[bus].visible = (bus === selected_bus);
        });
        """,
    )

    # Connect callback to dropdown
    bus_select.js_on_change("value", callback)

    # Create info text
    info_div = Div(
        text="""
        <p>This plot shows the carrier dispatch for different buses in the network.</p>
        <p>Positive values (above the x-axis) represent generation or inflow.</p>
        <p>Negative values (below the x-axis) represent consumption or outflow.</p>
        <p>Select a bus from the dropdown to view its carrier dispatch.</p>
        <p>Click on legend items to hide/show specific carriers.</p>
        """,
        width=400,
        styles={"font-size": "12px"},
    )

    # Create layout
    controls = bokeh_column(bus_select, info_div, width=400)
    layout = bokeh_row(controls, plot_container)

    return layout


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_bus_carrier_dispatch",
            clusters=10,
            opts="",
            sector_opts="",
            planning_horizons="2030",
        )

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(snakemake.log[0], mode="w"),
        ],
    )

    # Load network
    network = pypsa.Network(snakemake.input.network)
    logger.info(f"Loaded network with {len(network.buses)} buses")

    # Apply time subset if specified in the parameters
    if "snapshots" in snakemake.params and snakemake.params["snapshots"]:
        logger.info(f"Selecting time subset: {snakemake.params['snapshots']}")
        network.set_snapshots(network.snapshots[snakemake.params["snapshots"]])

    # Prepare data for all buses
    logger.info("Preparing carrier dispatch data for all buses")
    bus_data = prepare_all_buses_data(network)
    logger.info(f"Processed carrier dispatch for {len(bus_data)} buses")

    # Create output directory if it doesn't exist
    output_path = Path(snakemake.output.html)
    os.makedirs(output_path.parent, exist_ok=True)

    # Configure the output file
    output_file(snakemake.output.html, title="Bus Carrier Dispatch")

    # Create plot (starting with the first bus for now)
    if bus_data:
        default_bus = list(bus_data.keys())[0]
        logger.info(f"Creating plot with default bus: {default_bus}")
        layout = create_dispatch_plot(bus_data, default_bus)

        # Save the plot
        save(layout)
        logger.info(
            f"Interactive carrier dispatch plot saved to {snakemake.output.html}"
        )
    else:
        logger.error("No bus data available for plotting")
