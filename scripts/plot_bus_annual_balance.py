# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Plot annual carrier energy balance for different buses with interactive dropdown selection.
This shows the annual energy balance aggregated over the year rather than time series.
"""

import logging
import os
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from bokeh.io import output_file, save
# Import layout components with clear namespace
from bokeh.layouts import column as bokeh_column, row as bokeh_row
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


def get_annual_bus_balance(n, bus_name):
    """
    Calculate the annual energy balance for a specific bus.
    
    Parameters
    ----------
    n : pypsa.Network
        PyPSA network
    bus_name : str
        Name of the bus
        
    Returns
    -------
    dict
        Dictionary with carriers as keys and annual energy values as values
    """
    carriers = {}
    
    # Generation from generators
    generators = n.generators.index[n.generators.bus == bus_name]
    for gen in generators:
        carrier = n.generators.carrier[gen]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        # Sum over all time steps to get annual generation
        carriers[carrier] += n.generators_t.p[gen].sum()
    
    # Generation from storage units
    storage_units = n.storage_units.index[n.storage_units.bus == bus_name]
    for su in storage_units:
        carrier = n.storage_units.carrier[su]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        # Sum over all time steps
        carriers[carrier] += n.storage_units_t.p[su].sum()
    
    # Generation from links (only taking links with bus_name as bus1)
    for link in n.links.index[n.links.bus1 == bus_name]:
        carrier = n.links.carrier[link]
        if carrier not in carriers:
            carriers[carrier] = 0.0
        # Negative sign because bus1 is outflow, sum over all time steps
        carriers[carrier] += -n.links_t.p1[link].sum()
    
    # Load from links (only taking links with bus_name as bus0)
    for link in n.links.index[n.links.bus0 == bus_name]:
        carrier = n.links.carrier[link]
        eff = n.links.efficiency[link]
        p0 = n.links_t.p0[link]
        if f"{carrier} load" not in carriers:
            carriers[f"{carrier} load"] = 0.0
        # Negative sign for consumption, sum over all time steps
        carriers[f"{carrier} load"] += -p0.sum()
    
    # Conventional load
    loads = n.loads.index[n.loads.bus == bus_name]
    for load in loads:
        carrier = n.loads.carrier[load]
        if f"{carrier} load" not in carriers:
            carriers[f"{carrier} load"] = 0.0
        # Negative sign for consumption, sum over all time steps
        carriers[f"{carrier} load"] += -n.loads_t.p[load].sum()
    
    # Storage charging (consumption)
    stores = n.stores.index[n.stores.bus == bus_name]
    for store in stores:
        carrier = n.stores.carrier[store]
        if f"{carrier} charge" not in carriers:
            carriers[f"{carrier} charge"] = 0.0
        # Negative sign for consumption, sum over all time steps
        carriers[f"{carrier} charge"] += -n.stores_t.p[store].sum()
    
    return carriers


def prepare_all_buses_annual_data(n):
    """
    Prepare annual energy balance data for all buses in the network.
    
    Parameters
    ----------
    n : pypsa.Network
        PyPSA network
        
    Returns
    -------
    dict
        Dictionary with bus names as keys and carrier energy DataFrames as values
    """
    buses_data = {}
    
    # Process all buses that have generators, loads, or links connected
    active_buses = set()
    active_buses.update(n.generators.bus.unique())
    active_buses.update(n.loads.bus.unique())
    active_buses.update(n.links.bus0.unique())
    active_buses.update(n.links.bus1.unique())
    active_buses.update(n.storage_units.bus.unique())
    active_buses.update(n.stores.bus.unique())
    
    # For each active bus, calculate the annual carrier energy
    for bus in active_buses:
        try:
            annual_balance = get_annual_bus_balance(n, bus)
            
            # Convert to DataFrame
            df = pd.DataFrame({
                'carrier': list(annual_balance.keys()),
                'energy_mwh': list(annual_balance.values())
            })
            
            # Add 'type' column: generation (positive) or demand (negative)
            df['type'] = df['energy_mwh'].apply(
                lambda x: 'Generation' if x >= 0 else 'Demand'
            )
            
            # Take absolute values for better visualization
            df['energy_mwh_abs'] = df['energy_mwh'].abs()
            
            # Convert to TWh for better readability
            df['energy_twh'] = df['energy_mwh'] / 1e6
            df['energy_twh_abs'] = df['energy_mwh_abs'] / 1e6
            
            buses_data[bus] = df
        except Exception as e:
            logger.warning(f"Error processing bus {bus}: {e}")
    
    return buses_data


def create_annual_balance_plot(bus_data, bus_name):
    """
    Create an interactive Bokeh plot for annual carrier energy balance.
    
    Parameters
    ----------
    bus_data : dict
        Dictionary with bus names as keys and carrier energy DataFrames as values
    bus_name : str
        Name of the default bus to display
        
    Returns
    -------
    bokeh.layouts.column
        Bokeh layout with interactive plot
    """
    # Get list of all buses
    buses = sorted(bus_data.keys())
    
    # If the requested bus is not in the data, use the first available bus
    if bus_name not in buses:
        if buses:
            bus_name = buses[0]
        else:
            logger.error("No buses found with energy data")
            return None
    
    # Create one main figure that will display the balance for all buses
    p = figure(
        width=900,
        height=500,
        title=f"Annual Energy Balance",
        x_range=['Generation', 'Demand'],
        toolbar_location="above",
        tools=["pan", "wheel_zoom", "box_zoom", "reset", "save"],
    )
    
    # Style the plot
    p.title.text_font_size = "14pt"
    p.xaxis.axis_label = "Energy Flow Type"
    p.yaxis.axis_label = "Annual Energy (TWh)"
    p.axis.axis_label_text_font_style = "normal"
    p.grid.grid_line_alpha = 0.3
    p.legend.location = "right"
    p.legend.click_policy = "hide"
    
    # Add hover tooltips
    tooltips = [
        ("Carrier", "$name"),
        ("Energy (TWh/year)", "@$name{0.00}"),
    ]
    hover = HoverTool(tooltips=tooltips)
    p.add_tools(hover)
    
    # Create data sources dictionary for all buses
    sources = {}
    renderers = {}
    
    # Create stacked data for each bus
    for bus in buses:
        df = bus_data[bus]
        
        # Separate generation and demand
        gen_df = df[df['type'] == 'Generation'].sort_values('energy_twh_abs', ascending=False)
        dem_df = df[df['type'] == 'Demand'].sort_values('energy_twh_abs', ascending=False)
        
        # Generate colors with enough for all items
        gen_palette = Category20[20] * ((len(gen_df) // 20) + 1)
        dem_palette = Category20[20] * ((len(dem_df) // 20) + 1)
        
        gen_colors = list(gen_palette[:len(gen_df)]) if len(gen_df) > 0 else []
        dem_colors = list(dem_palette[:len(dem_df)]) if len(dem_df) > 0 else []
        
        # Create data structures for stacked bars
        gen_data = {'categories': ['Generation']}
        for i, row in gen_df.iterrows():
            carrier = row['carrier']
            gen_data[carrier] = [row['energy_twh']]
        
        dem_data = {'categories': ['Demand']}
        for i, row in dem_df.iterrows():
            carrier = row['carrier']
            # Store absolute values that will be made negative later
            dem_data[carrier] = [abs(row['energy_twh'])]
        
        # Create ColumnDataSources
        gen_source = ColumnDataSource(gen_data)
        dem_source = ColumnDataSource(dem_data)
        
        sources[bus] = {
            'gen': gen_source,
            'dem': dem_source,
            'gen_carriers': gen_df['carrier'].tolist(),
            'dem_carriers': dem_df['carrier'].tolist(),
            'gen_colors': gen_colors,
            'dem_colors': dem_colors
        }
        
        # Create renderers for this bus (they'll be hidden initially except for default bus)
        bus_renderers = {'gen': [], 'dem': []}
        
        # Generation bars
        if gen_df.empty:
            pass
        else:
            r = p.vbar_stack(
                stackers=gen_df['carrier'].tolist(),
                x='categories',
                width=0.6,
                source=gen_source,
                color=gen_colors,
                legend_label=gen_df['carrier'].tolist(),
                name=f'gen_{bus}',
                visible=(bus == bus_name)  # Only visible for default bus
            )
            bus_renderers['gen'] = r
        
        # Demand bars
        if dem_df.empty:
            pass
        else:
            r = p.vbar_stack(
                stackers=dem_df['carrier'].tolist(),
                x='categories',
                width=0.6,
                source=dem_source,
                color=dem_colors,
                legend_label=dem_df['carrier'].tolist(),
                name=f'dem_{bus}',
                visible=(bus == bus_name)  # Only visible for default bus
            )
            
            # Make demand values negative
            for renderer in r:
                field = renderer.name
                if field in dem_source.data:
                    for i in range(len(dem_source.data[field])):
                        dem_source.data[field][i] = -dem_source.data[field][i]
            
            # Update the data source
            dem_source.data = dict(dem_source.data)
        
        renderers[bus] = bus_renderers
    
    # Update plot title for default bus
    p.title.text = f"Annual Energy Balance for {bus_name}"
    
    # Create bus selector dropdown
    bus_select = Select(
        title="Select Bus:",
        value=bus_name,
        options=buses,
        width=400,
    )
    
    # Create a JavaScript callback for dropdown changes
    callback = CustomJS(
        args=dict(
            p=p,
            renderers=renderers,
            buses=buses
        ),
        code="""
        // Get the selected bus from the dropdown
        const selected_bus = cb_obj.value;
        
        // Update the plot title
        p.title.text = `Annual Energy Balance for ${selected_bus}`;
        
        // Set visibility of all renderers
        for (const bus of buses) {
            const is_visible = (bus === selected_bus);
            
            // Handle generation renderers
            if (renderers[bus]['gen']) {
                for (const r of renderers[bus]['gen']) {
                    r.visible = is_visible;
                }
            }
            
            // Handle demand renderers
            if (renderers[bus]['dem']) {
                for (const r of renderers[bus]['dem']) {
                    r.visible = is_visible;
                }
            }
        }
        """,
    )
    
    # Connect callback to dropdown
    bus_select.js_on_change('value', callback)
    
    # Create info text
    info_div = Div(
        text="""
        <p>This plot shows the annual energy balance for different buses in the network.</p>
        <p>Positive values (above the x-axis) represent generation or inflow in TWh/year.</p>
        <p>Negative values (below the x-axis) represent consumption or outflow in TWh/year.</p>
        <p>Select a bus from the dropdown to view its annual energy balance.</p>
        <p>Click on legend items to hide/show specific carriers.</p>
        """,
        width=400,
        styles={"font-size": "12px"},
    )
    
    # Create layout
    controls = bokeh_column(bus_select, info_div, width=400)
    layout = bokeh_row(controls, p)
    
    return layout


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_bus_carrier_dispatch",
            clusters=48,
            opts="Co2L0-24H",
            sector_opts="T-H-B-I-A-solar+p3-dist1",
            planning_horizons=2030,
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
    
    # Prepare annual data for all buses
    logger.info("Preparing annual energy balance data for all buses")
    bus_data = prepare_all_buses_annual_data(network)
    logger.info(f"Processed annual energy balance for {len(bus_data)} buses")
    
    # Create output directory if it doesn't exist
    output_path = Path(snakemake.output.html)
    os.makedirs(output_path.parent, exist_ok=True)
    
    # Configure the output file
    output_file(snakemake.output.html, title="Bus Annual Energy Balance")
    
    # Create plot (starting with the first bus for now)
    if bus_data:
        default_bus = list(bus_data.keys())[0]
        logger.info(f"Creating plot with default bus: {default_bus}")
        layout = create_annual_balance_plot(bus_data, default_bus)
        
        # Save the plot
        save(layout)
        logger.info(f"Interactive annual energy balance plot saved to {snakemake.output.html}")
    else:
        logger.error("No bus data available for plotting")