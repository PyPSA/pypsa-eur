# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import json
import logging
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pydeck as pdk
import pypsa
from shapely import wkt
from shapely.geometry import LineString, Polygon

from scripts._helpers import configure_logging, set_scenario_config
from scripts.base_network import _get_linetype_by_voltage

logger = logging.getLogger(__name__)


DISTANCE_CRS = "EPSG:3035"
GEO_CRS = "EPSG:4326"
LINE_SIMPLIFY = 100  # meters
STATIONS_SIMPLIFY = 100  # meters
BUSES_POLYGON_SIMPLIFY = 5  # meters
BUSES_COLUMNS = [
    "bus_id",
    "voltage",
    "dc",
    "symbol",
    "under_construction",
    "tags",
    "x",
    "y",
    "country",
    "geometry",
]
LINES_COLUMNS = [
    "line_id",
    "bus0",
    "bus1",
    "voltage",
    "i_nom",
    "circuits",
    "s_nom",
    "r",
    "x",
    "b",
    "length",
    "underground",
    "under_construction",
    "type",
    "tags",
    "geometry",
]
LINKS_COLUMNS = [
    "link_id",
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "length",
    "underground",
    "under_construction",
    "tags",
    "geometry",
]
TRANSFORMERS_COLUMNS = [
    "transformer_id",
    "bus0",
    "bus1",
    "voltage_bus0",
    "voltage_bus1",
    "s_nom",
    "geometry",
]
CONVERTERS_COLUMNS = [
    "converter_id",
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "geometry",
]
ROUNDING_DECIMALS = {
    "s_nom": 3,
    "r": 6,
    "x": 6,
    "b": 8,
}


def export_clean_csv(
    df: pd.DataFrame,
    columns: list[str],
    output_file: str | Path,
    rename_idx: str,
) -> pd.DataFrame:
    """
    Export a cleaned DataFrame to a CSV file.

    This function renames the DataFrame index, applies column renaming,
    selects specified columns, converts boolean values to 't'/'f', and
    exports the result to a CSV file.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to be exported.
    columns : list[str]
        Column names to include in the exported CSV file.
    output_file : str | Path
        Path to the output CSV file.
    rename_idx : str
        Name to use for the renamed index column.

    Returns
    -------
    df : pd.DataFrame
        The cleaned DataFrame that was exported.

    """

    rename_dict = {
        "v_nom": "voltage",
        "num_parallel": "circuits",
    }

    df = (
        df.rename_axis(index=rename_idx)
        .reset_index()
        .rename(columns=rename_dict)
        .loc[:, columns]
        .replace({True: "t", False: "f"})
    )

    df.to_csv(output_file, index=False, quotechar="'")

    return df


def linestring_to_coords(linestring: LineString) -> list[list[float]]:
    """
    Convert shapely LineString to list of [lon, lat] coordinates.

    Parameters
    ----------
    linestring : LineString
        A shapely LineString geometry object.

    Returns
    -------
    list[list[float]]
        List of [longitude, latitude] coordinate pairs.

    """
    return [[coord[0], coord[1]] for coord in linestring.coords]


def polygon_to_coords(polygon: Polygon) -> list[list[float]]:
    """
    Convert shapely Polygon to coordinate list for PyDeck.

    Parameters
    ----------
    polygon : Polygon
        A shapely Polygon geometry object.

    Returns
    -------
    list[list[float]]
        List of [longitude, latitude] coordinate pairs from exterior ring.

    """
    return [[coord[0], coord[1]] for coord in polygon.exterior.coords]


def get_line_colors(voltages: pd.Series) -> list:
    """
    Get RGB colors for transmission lines based on voltage level ranges (vectorized).

    Uses discrete color ranges: yellow (low) → orange → red (high) → magenta (ultra-high).
    DC lines are colored neon green regardless of voltage.

    Parameters
    ----------
    voltages : pd.Series
        Voltage levels in kV.

    Returns
    -------
    list
        List of RGBA colors as [R, G, B, A] arrays.

    """
    voltages = voltages.to_numpy()
    n = len(voltages)
    colors = np.zeros((n, 4), dtype=int)
    colors[:, 3] = 200  # Alpha channel

    # Define voltage ranges and their colors [R, G, B]
    # Format: (min_voltage, max_voltage, [R, G, B])
    voltage_ranges = [
        (0, 110, [255, 255, 0]),  # < 110 kV: Bright yellow
        (110, 150, [255, 235, 0]),  # 110-150 kV: Light yellow
        (151, 240, [255, 200, 0]),  # 151-240 kV: Gold/Orange-yellow
        (241, 280, [255, 165, 0]),  # 241-280 kV: Orange
        (281, 330, [255, 120, 0]),  # 281-330 kV: Dark orange
        (331, 350, [255, 69, 0]),  # 331-350 kV: Orange-red
        (351, 420, [255, 40, 0]),  # 351-420 kV: Red-orange
        (421, 520, [220, 20, 60]),  # 421-520 kV: Crimson
        (521, np.inf, [255, 0, 255]),  # > 520 kV: Magenta
    ]

    # Apply colors based on voltage ranges
    for v_min, v_max, rgb in voltage_ranges:
        if v_max == np.inf:
            mask = voltages >= v_min
        else:
            mask = (voltages >= v_min) & (voltages <= v_max)
        colors[mask, :3] = rgb

    return colors.tolist()


def inject_custom_controls(deck: pdk.Deck, release_version: str) -> str:  # noqa: W291, W293
    """
    Inject custom layer visibility controls into a pydeck object.

    Adds interactive checkboxes to toggle visibility of different network
    component layers and voltage filtering with tag-based selection.

    Parameters
    ----------
    deck : pdk.Deck
        pydeck Deck instance to export with custom controls.
    release_version : str
        Release version string to include in HTML title.

    Returns
    -------
    str
        Updated html file in string format with injected controls.

    """
    html = deck.to_html(as_string=True)

    # Expose deck instance to window for JavaScript access
    html = html.replace(
        "});\n\n  </script>",
        "});\nwindow.deck = deckInstance;\n\n  </script>",
    )

    # HTML/JavaScript for collapsible layer visibility controls
    # ruff: noqa
    controls = """
        <style>
        .voltage-tag {
            display: inline-block;
            background: #6c757d;
            color: white;
            padding: 4px 8px;
            margin: 2px;
            border-radius: 12px;
            font-size: 12px;
            cursor: pointer;
        }
        .voltage-tag:hover {
            background: #0056b3;
        }
        .voltage-tag.selected {
            background: #28a745;
        }
        .selected-voltage-tag {
            display: inline-block;
            background: #28a745;
            color: white;
            padding: 4px 8px;
            margin: 2px;
            border-radius: 12px;
            font-size: 12px;
        }
        .selected-voltage-tag .remove {
            margin-left: 6px;
            cursor: pointer;
            font-weight: bold;
        }
        #voltage-suggestions {
            max-height: 150px;
            overflow-y: auto;
            border: 1px solid #ddd;
            border-radius: 3px;
            padding: 6px;
            margin-top: 4px;
            background: white;
        }
        </style>
        <button id="menu-toggle" 
            style="position:absolute;top:10px;left:10px;
                   background:rgba(255,255,255,0.9);border:none;
                   padding:10px;cursor:pointer;font-size:20px;
                   border-radius:4px;z-index:1001;box-shadow:0 2px 4px rgba(0,0,0,0.2)">
            ☰
        </button>
        <div id="layer-controls"
            style="position:absolute;top:10px;left:10px;
                   background:rgba(255,255,255,0.95);padding:12px;
                   font-family:sans-serif;z-index:1000;
                   border-radius:4px;box-shadow:0 2px 8px rgba(0,0,0,0.3);
                   display:none;margin-top:50px;min-width:250px;max-width:300px">
            
            <!-- Voltage Filter -->
            <div style="margin-bottom:12px;padding-bottom:12px;border-bottom:1px solid #ddd">
                <label style="display:block;margin-bottom:4px;font-size:13px">
                    Filter by voltage (kV):
                </label>
                <div id="selected-voltages" style="min-height:20px;margin-bottom:4px"></div>
                <input type="text" id="voltage-filter" placeholder="Type to search..."
                    style="width:100%;padding:6px;border:1px solid #ccc;border-radius:3px;
                           font-size:13px;box-sizing:border-box">
                <div id="voltage-suggestions"></div>
                <button onclick="clearVoltageFilter()"
                    style="width:100%;margin-top:8px;padding:6px;background:#6c757d;
                           color:white;border:none;border-radius:3px;cursor:pointer;
                           font-size:13px">
                    Clear filter
                </button>
            </div>
            
            <!-- Layer Toggles -->
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Lines', this.checked)"
                        style="margin-right:8px">
                    Lines
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Links', this.checked)"
                        style="margin-right:8px">
                    Links
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Converters', this.checked)"
                        style="margin-right:8px">
                    Converters
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Transformers', this.checked)"
                        style="margin-right:8px">
                    Transformers
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Buses', this.checked)"
                        style="margin-right:8px">
                    Buses
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Buses (Polygons)', this.checked)"
                        style="margin-right:8px">
                    Buses (Polygons)
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Stations', this.checked)"
                        style="margin-right:8px">
                    Stations
                </label>
            </div>
        </div>
        <script>
        // Store original data and available voltages
        let originalLayerData = {};
        let availableVoltages = new Set();
        let selectedVoltages = new Set();
        
        // Extract all unique voltages from dataset
        function initializeVoltages() {
            const deck = window.deck;
            if (!deck) return;
            
            deck.props.layers.forEach(layer => {
                originalLayerData[layer.id] = layer.props.data;
                
                layer.props.data.forEach(item => {
                    if (item.voltage !== undefined) {
                        availableVoltages.add(item.voltage);
                    }
                    if (item.voltage_bus0 !== undefined) {
                        availableVoltages.add(item.voltage_bus0);
                    }
                    if (item.voltage_bus1 !== undefined) {
                        availableVoltages.add(item.voltage_bus1);
                    }
                });
            });
            
            // Sort voltages in descending order
            availableVoltages = new Set([...availableVoltages].sort((a, b) => b - a));
            updateSuggestions();
        }
        
        function updateSuggestions(searchTerm = '') {
            const container = document.getElementById('voltage-suggestions');
            container.innerHTML = '';
            
            const filteredVoltages = [...availableVoltages].filter(v => 
                String(v).includes(searchTerm) && !selectedVoltages.has(v)
            );
            
            if (filteredVoltages.length === 0 && searchTerm === '') {
                container.innerHTML = '<div style="color:#999;font-size:12px">No more voltages available</div>';
                return;
            }
            
            if (filteredVoltages.length === 0) {
                container.innerHTML = '<div style="color:#999;font-size:12px">No matches</div>';
                return;
            }
            
            filteredVoltages.forEach(voltage => {
                const tag = document.createElement('span');
                tag.className = 'voltage-tag';
                tag.textContent = voltage;
                tag.onclick = (e) => {
                    e.stopPropagation();  // Prevent menu from closing
                    addVoltage(voltage);
                };
                container.appendChild(tag);
            });
        }
        
        function addVoltage(voltage) {
            selectedVoltages.add(voltage);
            updateSelectedDisplay();
            updateSuggestions();
            applyVoltageFilter();
            document.getElementById('voltage-filter').value = '';
        }
        
        function removeVoltage(voltage, event) {
            if (event) event.stopPropagation();  // Prevent menu from closing
            selectedVoltages.delete(voltage);
            updateSelectedDisplay();
            updateSuggestions();
            applyVoltageFilter();
        }
        
        function updateSelectedDisplay() {
            const container = document.getElementById('selected-voltages');
            container.innerHTML = '';
            
            [...selectedVoltages].sort((a, b) => b - a).forEach(voltage => {
                const tag = document.createElement('span');
                tag.className = 'selected-voltage-tag';
                const removeBtn = document.createElement('span');
                removeBtn.className = 'remove';
                removeBtn.textContent = '×';
                removeBtn.onclick = (e) => removeVoltage(voltage, e);
                
                tag.textContent = voltage;
                tag.appendChild(removeBtn);
                container.appendChild(tag);
            });
        }
        
        function setLayerVisibility(layerId, visible) {
            const deck = window.deck;
            if (!deck) return;
            const layers = deck.props.layers.map(l =>
                l.id === layerId ? l.clone({ visible }) : l
            );
            deck.setProps({ layers });
        }
        
        function applyVoltageFilter() {
            const deck = window.deck;
            if (!deck) return;
            
            if (selectedVoltages.size === 0) {
                // No filter - show all data
                const layers = deck.props.layers.map(layer => {
                    if (originalLayerData[layer.id]) {
                        return layer.clone({ data: originalLayerData[layer.id] });
                    }
                    return layer;
                });
                deck.setProps({ layers });
                return;
            }
            
            const voltages = [...selectedVoltages];
            
            const layers = deck.props.layers.map(layer => {
                const originalData = originalLayerData[layer.id];
                const filteredData = originalData.filter(item => {
                    if (item.voltage !== undefined) {
                        return voltages.includes(item.voltage);
                    }
                    if (item.voltage_bus0 !== undefined || item.voltage_bus1 !== undefined) {
                        return voltages.includes(item.voltage_bus0) || 
                               voltages.includes(item.voltage_bus1);
                    }
                    return true;
                });
                
                return layer.clone({ data: filteredData });
            });
            
            deck.setProps({ layers });
        }
        
        function clearVoltageFilter() {
            selectedVoltages.clear();
            updateSelectedDisplay();
            updateSuggestions();
            applyVoltageFilter();
            document.getElementById('voltage-filter').value = '';
        }
        
        // Search input handler
        document.addEventListener('DOMContentLoaded', function() {
            initializeVoltages();
            
            const input = document.getElementById('voltage-filter');
            if (input) {
                input.addEventListener('input', function(e) {
                    updateSuggestions(e.target.value);
                });
            }
        });
        
        // Toggle menu visibility
        document.getElementById('menu-toggle').addEventListener('click', function() {
            const menu = document.getElementById('layer-controls');
            menu.style.display = menu.style.display === 'none' ? 'block' : 'none';
        });
        
        // Close menu when clicking outside
        document.addEventListener('click', function(event) {
            const menu = document.getElementById('layer-controls');
            const toggle = document.getElementById('menu-toggle');
            if (menu.style.display === 'block' && 
                !menu.contains(event.target) && 
                !toggle.contains(event.target)) {
                menu.style.display = 'none';
            }
        });
        </script>
    """

    # Inject controls before deck container
    html = html.replace(
        '<div id="deck-container">',
        controls + '<div id="deck-container">',
    )

    # Inject custom title
    html = html.replace(
        "<title>pydeck</title>",
        "<title>PyPSA network (OSM " + release_version + ")</title>",
    )

    return html


def compress_html(html: str) -> str:
    """
    Compress HTML by removing unnecessary whitespace and comments.

    Reduces file size by minifying HTML, CSS, and JavaScript content
    while preserving functionality. Aggressively compresses JSON data.

    Parameters
    ----------
    html : str
        HTML string to compress.

    Returns
    -------
    str
        Compressed HTML string.

    """
    # Remove HTML comments
    html = re.sub(r"<!--.*?-->", "", html, flags=re.DOTALL)

    # Minify CSS only
    def minify_css(match):
        css = match.group(1)
        css = re.sub(r"\s+", " ", css)
        css = re.sub(r"\s*([{};:,])\s*", r"\1", css)
        css = css.strip()
        return f"<style>{css}</style>"

    html = re.sub(r"<style>(.*?)</style>", minify_css, html, flags=re.DOTALL)

    # Compress JSON data in the main script tag
    def compress_json_in_script(match):
        script_content = match.group(0)  # Get the entire match including <script> tags

        # Find the jsonInput declaration
        json_match = re.search(r"const jsonInput = ({.*?});", script_content, re.DOTALL)
        if json_match:
            json_str = json_match.group(1)

            # Parse and re-serialize without whitespace
            try:
                json_obj = json.loads(json_str)
                compact_json = json.dumps(json_obj, separators=(",", ":"))
                script_content = script_content.replace(json_str, compact_json)
            except json.JSONDecodeError:
                pass  # If parsing fails, leave as is

        return script_content

    # Apply to the main script (not the controls script)
    html = re.sub(
        r"<script>\s*const container = document\.getElementById.*?</script>",
        compress_json_in_script,
        html,
        flags=re.DOTALL,
    )

    # Remove extra whitespace in HTML (but NOT in script tags)
    parts = re.split(r"(<script>.*?</script>)", html, flags=re.DOTALL | re.IGNORECASE)

    for i in range(len(parts)):
        if not parts[i].startswith("<script>") and not parts[i].lower().startswith(
            "<script>"
        ):
            parts[i] = re.sub(r">\s+<", "><", parts[i])
            parts[i] = re.sub(r"\n\s+", "\n", parts[i])

    html = "".join(parts)

    return html


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_osm_network_release",
            configfiles=["config/config.osm-release.yaml"],
        )

    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    # Params and import
    line_types = snakemake.params.line_types
    release_version = snakemake.params.release_version
    n = pypsa.Network(snakemake.input.base_network)

    #############
    ### Buses ###
    #############

    logger.info("Cleaning and formatting bus data for release.")
    buses = n.buses.copy()
    buses["dc"] = buses.pop("carrier").map({"DC": "t", "AC": "f"})
    buses.v_nom = buses.v_nom.astype(int)
    buses.sort_index(inplace=True)

    logger.info(f"Exporting {len(buses)} buses to %s", snakemake.output.buses)
    buses = export_clean_csv(buses, BUSES_COLUMNS, snakemake.output.buses, "bus_id")

    #############
    ### Lines ###
    #############

    logger.info("Cleaning and formatting line data for release.")
    n.lines.loc[:, "type"] = n.lines.v_nom.apply(
        lambda x: _get_linetype_by_voltage(x, line_types)
    )
    n.calculate_dependent_values()  # Calculate dependent variables (r, x)
    lines = n.lines.copy()
    lines["i_nom"] = (
        (lines.s_nom / lines.v_nom / lines.num_parallel).div(np.sqrt(3)).round(3)
    )  # kA
    lines.v_nom = lines.v_nom.astype(int)
    lines.num_parallel = lines.num_parallel.astype(int)
    lines.length = lines.length * 1e3
    lines[list(ROUNDING_DECIMALS.keys())] = lines[list(ROUNDING_DECIMALS.keys())].round(
        ROUNDING_DECIMALS
    )
    lines.sort_index(inplace=True)

    logger.info(f"Exporting {len(lines)} lines to %s", snakemake.output.lines)
    lines = export_clean_csv(lines, LINES_COLUMNS, snakemake.output.lines, "line_id")

    ##########################
    ### Links + Converters ###
    ##########################

    logger.info("Cleaning and formatting link and converter data for release.")
    links = n.links.copy()
    links.voltage = links.voltage.astype(int)
    links.p_nom = links.p_nom.astype(int)
    links.length = links.length * 1e3
    links.sort_index(inplace=True)

    # Boolean that specifies if link element is a converter
    is_converter = links.index.str.startswith("conv") == True
    links_dc = links[~is_converter].copy()
    converters = links[is_converter].copy()

    logger.info(
        f"Exporting {len(links_dc)} links to %s",
        snakemake.output.links,
    )
    links_dc = export_clean_csv(
        links_dc, LINKS_COLUMNS, snakemake.output.links, "link_id"
    )

    logger.info(
        f"Exporting {len(converters)} converters to %s",
        snakemake.output.converters,
    )
    converters = export_clean_csv(
        converters,
        CONVERTERS_COLUMNS,
        snakemake.output.converters,
        "converter_id",
    )

    ####################
    ### Transformers ###
    ####################

    logger.info("Cleaning and formatting transformer data for release.")
    transformers = n.transformers.copy()
    transformers.voltage_bus0 = transformers.voltage_bus0.astype(int)
    transformers.voltage_bus1 = transformers.voltage_bus1.astype(int)
    transformers.s_nom = transformers.s_nom.astype(int)
    transformers.sort_index(inplace=True)

    logger.info(
        f"Exporting {len(transformers)} transformers to %s",
        snakemake.output.transformers,
    )
    transformers = export_clean_csv(
        transformers,
        TRANSFORMERS_COLUMNS,
        snakemake.output.transformers,
        "transformer_id",
    )

    ######################################
    ### Pydeck map with layer controls ###
    ######################################

    buses["geometry"] = buses["geometry"].apply(wkt.loads)
    buses = gpd.GeoDataFrame(buses, geometry="geometry", crs=GEO_CRS)

    # Import stations
    stations_polygon = gpd.read_file(snakemake.input.stations_polygon).to_crs(GEO_CRS)

    # Only keep stations_polygon that contain buses points
    stations_polygon = gpd.sjoin(
        stations_polygon, buses, how="left", predicate="contains"
    )
    stations_polygon = stations_polygon[stations_polygon.index_right.notnull()]
    stations_polygon = stations_polygon.drop_duplicates(subset=["station_id"])
    stations_polygon = stations_polygon[["station_id", "geometry"]]

    # Simplify
    stations_polygon = stations_polygon.to_crs(DISTANCE_CRS)
    stations_polygon["geometry"] = stations_polygon["geometry"].simplify(
        STATIONS_SIMPLIFY
    )
    stations_polygon = stations_polygon.to_crs(GEO_CRS)

    stations_polygon["poly"] = stations_polygon["geometry"].apply(polygon_to_coords)

    # Import bus polygons
    buses_polygon = gpd.read_file(snakemake.input.buses_polygon).to_crs(GEO_CRS)

    # Only keep bus polygons that contain buses points
    buses_polygon = gpd.sjoin(
        buses_polygon,
        buses,
        how="left",
        predicate="contains",
    )
    buses_polygon = buses_polygon[buses_polygon.index_right.notnull()]
    buses_polygon = buses_polygon.drop_duplicates(subset=["bus_id_left"])
    buses_polygon = buses_polygon[["geometry"]]

    # Simplify
    buses_polygon = buses_polygon.to_crs(DISTANCE_CRS)
    buses_polygon["geometry"] = buses_polygon["geometry"].simplify(
        BUSES_POLYGON_SIMPLIFY
    )
    buses_polygon = buses_polygon.to_crs(GEO_CRS)

    buses_polygon["poly"] = buses_polygon["geometry"].apply(polygon_to_coords)

    # Prepare geometries for pydeck
    # Lines
    lines["geometry"] = lines["geometry"].apply(wkt.loads)
    lines = gpd.GeoDataFrame(lines, geometry="geometry", crs=GEO_CRS)
    lines = lines.to_crs(DISTANCE_CRS)
    lines["geometry"] = lines["geometry"].simplify(LINE_SIMPLIFY)
    lines = lines.to_crs(GEO_CRS)
    lines["path"] = lines["geometry"].apply(linestring_to_coords)
    lines["color"] = get_line_colors(lines["voltage"])

    # Links
    links_dc["path"] = links_dc["geometry"].apply(
        lambda wkt_str: linestring_to_coords(wkt.loads(wkt_str))
    )

    # Converters
    converters["path"] = converters["geometry"].apply(
        lambda wkt_str: linestring_to_coords(wkt.loads(wkt_str))
    )

    # Transformers
    transformers["path"] = transformers["geometry"].apply(
        lambda wkt_str: linestring_to_coords(wkt.loads(wkt_str))
    )

    logger.info("Creating interactive map with pydeck.")

    # Stations PolygonLayer
    stations_polygon_layer = pdk.Layer(
        "PolygonLayer",
        data=stations_polygon.drop(columns=["geometry"]),
        get_polygon="poly",
        get_fill_color=[0, 0, 255, 100],
        pickable=True,
        auto_highlight=True,
        id="Stations",
        parameters={"depthTest": False},
    )

    # Buses PolygonLayer
    buses_polygon_layer = pdk.Layer(
        "PolygonLayer",
        data=buses_polygon.drop(columns=["geometry"]),
        get_polygon="poly",
        get_fill_color=[255, 0, 255, 100],
        get_line_color=[255, 255, 255],
        extruded=True,
        wireframe=True,
        get_elevation=200,
        pickable=True,
        auto_highlight=True,
        id="Buses (Polygons)",
        parameters={"depthTest": False},
    )

    # Lines PathLayer
    lines_layer = pdk.Layer(
        "PathLayer",
        data=lines.drop(columns=["geometry"]),
        get_path="path",
        get_color="color",
        width_scale=1,
        width_min_pixels=2,
        pickable=True,
        auto_highlight=True,
        parameters={"depthTest": False},
        id="Lines",
    )

    # Links PathLayer
    links_layer = pdk.Layer(
        "PathLayer",
        data=links_dc.drop(columns=["geometry"]),
        get_path="path",
        get_color=[0, 255, 0, 160],
        width_scale=1,
        width_min_pixels=2,
        pickable=True,
        auto_highlight=True,
        parameters={"depthTest": False},
        id="Links",
    )

    # Converters PathLayer
    converters_layer = pdk.Layer(
        "PathLayer",
        data=converters.drop(columns=["geometry"]),
        get_path="path",
        get_color=[255, 100, 100, 160],
        width_scale=1,
        width_min_pixels=2,
        pickable=True,
        auto_highlight=True,
        parameters={"depthTest": False},
        id="Converters",
    )

    # Transformers PathLayer
    transformers_layer = pdk.Layer(
        "PathLayer",
        data=transformers.drop(columns=["geometry"]),
        get_path="path",
        get_color=[255, 255, 0, 160],
        width_scale=1,
        width_min_pixels=2,
        pickable=True,
        auto_highlight=True,
        parameters={"depthTest": False},
        id="Transformers",
    )

    # Buses ScatterplotLayer
    buses_layer = pdk.Layer(
        "ColumnLayer",
        data=buses,
        get_position=["x", "y"],
        get_fill_color=[255, 0, 155, 160],
        radius=20,
        get_elevation=195,
        pickable=True,
        auto_highlight=True,
        parameters={"depthTest": False},
        id="Buses",
    )

    # Calculate center point for initial view
    all_coords = [coord for path in lines["path"] for coord in path]
    center_lon = sum(c[0] for c in all_coords) / len(all_coords)
    center_lat = sum(c[1] for c in all_coords) / len(all_coords)

    # Create the deck
    map = pdk.Deck(
        layers=[
            stations_polygon_layer,
            buses_polygon_layer,
            lines_layer,
            links_layer,
            converters_layer,
            transformers_layer,
            buses_layer,
        ],
        initial_view_state=pdk.ViewState(
            latitude=center_lat,
            longitude=center_lon,
            zoom=5,
            pitch=30,
        ),
    )

    logger.info("Injecting custom layer controls into map HTML.")
    map_ctrl = inject_custom_controls(map, release_version)

    # logger.info("Compressing map HTML to reduce file size.")
    map_ctrl = compress_html(map_ctrl)

    # Export to HTML
    logger.info("Exporting interactive map to %s", snakemake.output.map)
    with open(snakemake.output.map, "w") as f:
        f.write(map_ctrl)
