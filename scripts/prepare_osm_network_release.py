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
LINE_SIMPLIFY = 30  # meters
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
    export: bool = True,
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
    export : bool, optional
        Whether to perform store the cleaned CSV in the output directory, by default True.

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

    if export:
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
    colors[:, 3] = 150  # Alpha channel

    # Format: (min_voltage, max_voltage, [R, G, B])
    voltage_ranges = [
        (0, 110, [255, 255, 0]),  # < 110 kV: Bright yellow
        (110, 150, [255, 235, 0]),  # 110-150 kV: Light yellow
        (151, 240, [0, 250, 0]),  # 151-240 kV: Green
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


def inject_custom_controls(deck: pdk.Deck, release_version: str) -> str:
    """
    Inject interactive controls and URL-based state management into PyDeck HTML.

    Adds layer toggles, voltage/text filtering, theme switching, click-based
    tooltips with OSM links, and URL hash synchronization for map state.

    Parameters
    ----------
    deck : pdk.Deck
        PyDeck Deck instance to export with custom controls.
    release_version : str
        Release version string to include in HTML title.

    Returns
    -------
    str
        HTML string with injected interactive controls and state management.

    """
    html = deck.to_html(as_string=True)

    html = html.replace(
        "});\n\n  </script>",
        "});\nwindow.deck = deckInstance;\n\n  </script>",
    )

    controls = """
        <meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no">
        <style>
        html, body {
            margin: 0;
            padding: 0;
            width: 100%;
            height: 100%;
            overflow: hidden;
        }
        #deck-container {
            position: fixed !important;
            inset: 0 !important;
            width: 100% !important;
            height: 100% !important;
        }
        #deck-container canvas {
            width: 100% !important;
            height: 100% !important;
        }
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
        .map-btn {
            position: absolute;
            top: 10px;
            background: rgba(255,255,255,0.7);
            border: none;
            cursor: pointer;
            border-radius: 4px;
            z-index: 10001;
            box-shadow: 0 2px 4px rgba(0,0,0,0.2);
            line-height: 1;
            height: 40px;
            width: 40px;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: background 0.15s, box-shadow 0.15s, transform 0.1s;
        }
        .map-btn:hover {
            background: rgba(255,255,255,0.95);
            box-shadow: 0 3px 8px rgba(0,0,0,0.3);
        }
        .map-btn:active {
            background: rgba(220,220,220,0.9);
            box-shadow: 0 1px 2px rgba(0,0,0,0.2);
            transform: scale(0.93);
        }
        .ctrl-btn {
            width: 100%;
            margin-top: 8px;
            padding: 6px;
            background: #6c757d;
            color: white;
            border: none;
            border-radius: 3px;
            cursor: pointer;
            font-size: 13px;
            transition: background 0.15s, box-shadow 0.15s, transform 0.1s;
        }
        .ctrl-btn:hover {
            background: #5a6268;
            box-shadow: 0 2px 6px rgba(0,0,0,0.25);
        }
        .ctrl-btn:active {
            background: #4e555b;
            box-shadow: 0 1px 2px rgba(0,0,0,0.2);
            transform: scale(0.97);
        }
        </style>
        <button id="menu-toggle" class="map-btn"
            style="left:10px;font-size:16px;padding:10px 12px">
            ☰
        </button>
        <button id="theme-toggle" class="map-btn"
            style="left:60px;font-size:16px;padding:10px 12px">
            ◐
        </button>
        <button id="zoom-in" class="map-btn"
            onclick="zoomIn()"
            style="left:110px;font-size:18px;padding:10px 12px">
            +
        </button>
        <button id="zoom-out" class="map-btn"
            onclick="zoomOut()"
            style="left:160px;font-size:18px;padding:10px 12px">
            −
        </button>
        <div id="layer-controls"
            style="position:absolute;top:10px;left:10px;
                background:rgba(255,255,255,0.75);padding:12px;
                font-family:sans-serif;z-index:10001;
                border-radius:4px;box-shadow:0 2px 8px rgba(0,0,0,0.3);
                display:none;margin-top:50px;min-width:250px;max-width:300px">

            <div style="margin-bottom:12px;padding-bottom:12px;border-bottom:1px solid #ddd">
                <label style="display:block;margin-bottom:4px;font-size:13px">
                    Search all fields:
                </label>
                <input type="text" id="text-search-filter" placeholder="term1 & term2 | term3"
                    style="width:100%;padding:6px;border:1px solid #ccc;border-radius:3px;
                        font-size:13px;box-sizing:border-box">
                <div style="font-size:11px;color:#666;margin-top:4px">
                    Use & for AND, | for OR
                </div>
                <button class="ctrl-btn" onclick="clearTextSearchFilter()">
                    Clear search
                </button>
            </div>

            <div style="margin-bottom:12px;padding-bottom:12px;border-bottom:1px solid #ddd">
                <label style="display:block;margin-bottom:4px;font-size:13px">
                    Filter by voltage (kV):
                </label>
                <div id="selected-voltages" style="min-height:20px;margin-bottom:4px"></div>
                <input type="text" id="voltage-filter" placeholder="Type to search..."
                    style="width:100%;padding:6px;border:1px solid #ccc;border-radius:3px;
                        font-size:13px;box-sizing:border-box">
                <div id="voltage-suggestions"></div>
                <button class="ctrl-btn" onclick="clearVoltageFilter()">
                    Clear filter
                </button>
            </div>

            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center;font-size:12px">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Lines', this.checked)"
                        style="margin-right:8px">
                    Lines
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center;font-size:12px">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Links', this.checked)"
                        style="margin-right:8px">
                    Links
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center;font-size:12px">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Converters', this.checked)"
                        style="margin-right:8px">
                    Converters
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center;font-size:12px">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Transformers', this.checked)"
                        style="margin-right:8px">
                    Transformers
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center;font-size:12px">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Buses', this.checked)"
                        style="margin-right:8px">
                    Buses
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center;font-size:12px">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Buses (Polygons)', this.checked)"
                        style="margin-right:8px">
                    Buses (Polygons)
                </label>
            </div>
            <div style="margin:6px 0">
                <label style="cursor:pointer;display:flex;align-items:center;font-size:12px">
                    <input type="checkbox" checked
                        onchange="setLayerVisibility('Stations', this.checked)"
                        style="margin-right:8px">
                    Stations
                </label>
            </div>
        </div>
        <script>
        let originalLayerData = {};
        let availableVoltages = new Set();
        let selectedVoltages = new Set();
        let currentTextSearch = '';
        let isDarkMode = true;
        let hashUpdateTimeout = null;
        let isMobile = /iPhone|iPad|iPod|Android/i.test(navigator.userAgent);

        const MAP_STYLES = {
            light: 'https://basemaps.cartocdn.com/gl/positron-gl-style/style.json',
            dark: 'https://basemaps.cartocdn.com/gl/dark-matter-gl-style/style.json'
        };

        function parseHash() {
            const hash = window.location.hash.substring(1);
            if (!hash) return null;

            const parts = hash.split('/');
            if (parts.length >= 4) {
                return {
                    theme: parts[0],
                    zoom: parseFloat(parts[1]),
                    latitude: parseFloat(parts[2]),
                    longitude: parseFloat(parts[3])
                };
            }
            return null;
        }

        function updateHash(viewState, immediate = false) {
            if (hashUpdateTimeout) {
                clearTimeout(hashUpdateTimeout);
            }

            const doUpdate = () => {
                const theme = isDarkMode ? 'dark' : 'light';
                const hash = `#${theme}/${viewState.zoom.toFixed(2)}/${viewState.latitude.toFixed(4)}/${viewState.longitude.toFixed(4)}`;
                window.history.replaceState(null, '', hash);
            };

            if (immediate) {
                doUpdate();
            } else {
                hashUpdateTimeout = setTimeout(doUpdate, 150);
            }
        }

        function applyHashToMap() {
            const hashParams = parseHash();
            if (!hashParams || !window.deck) return;

            window.deck.setProps({
                initialViewState: {
                    latitude: hashParams.latitude,
                    longitude: hashParams.longitude,
                    zoom: hashParams.zoom,
                    pitch: 30,
                    transitionDuration: 500,
                }
            });

            if (hashParams.theme && hashParams.theme !== (isDarkMode ? 'dark' : 'light')) {
                isDarkMode = hashParams.theme === 'dark';
                document.getElementById('theme-toggle').textContent = isDarkMode ? '◐' : '○';
                window.deck.setProps({
                    mapStyle: isDarkMode ? MAP_STYLES.dark : MAP_STYLES.light
                });
            }
        }

        function toggleTheme() {
            isDarkMode = !isDarkMode;
            const deck = window.deck;
            if (!deck) return;

            document.getElementById('theme-toggle').textContent = isDarkMode ? '◐' : '○';
            deck.setProps({
                mapStyle: isDarkMode ? MAP_STYLES.dark : MAP_STYLES.light
            });

            const viewState = deck.viewManager.getViewports()[0];
            if (viewState) updateHash(viewState, true);
        }

        function zoomIn() {
            const deck = window.deck;
            if (!deck) return;
            const vp = deck.viewManager.getViewports()[0];
            if (!vp) return;
            deck.setProps({
                initialViewState: {
                    latitude: vp.latitude,
                    longitude: vp.longitude,
                    zoom: Math.min(vp.zoom + 1, 20),
                    pitch: vp.pitch || 30,
                    transitionDuration: 300,
                }
            });
        }

        function zoomOut() {
            const deck = window.deck;
            if (!deck) return;
            const vp = deck.viewManager.getViewports()[0];
            if (!vp) return;
            deck.setProps({
                initialViewState: {
                    latitude: vp.latitude,
                    longitude: vp.longitude,
                    zoom: Math.max(vp.zoom - 1, 1),
                    pitch: vp.pitch || 30,
                    transitionDuration: 300,
                }
            });
        }

        const CIRCUIT_OFFSET_METERS = 0.0003;
        const CIRCUIT_ZOOM_THRESHOLD = 11.5;
        const TAPER_DEGREES = 0.00025;  // ~20 meters at European latitudes
        let circuitsExpanded = false;
        let hoveredLineId = null;

        function pathLength(path) {
            let total = 0;
            for (let i = 1; i < path.length; i++) {
                const dx = path[i][0] - path[i-1][0];
                const dy = path[i][1] - path[i-1][1];
                total += Math.sqrt(dx*dx + dy*dy);
            }
            return total;
        }

        function interpCoord(a, b, t) {
            return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t];
        }

        // Returns [densePath, totalLength] so callers avoid recomputing pathLength.
        function insertTaperPoints(path, taper, total) {
            if (total <= taper * 2) return [path, total];
            const result = [];
            let cumDist = 0;
            const startT = taper, endT = total - taper;
            let startIns = false, endIns = false;
            for (let i = 0; i < path.length; i++) {
                if (i > 0) {
                    const dx = path[i][0] - path[i-1][0];
                    const dy = path[i][1] - path[i-1][1];
                    const seg = Math.sqrt(dx*dx + dy*dy);
                    if (!startIns && cumDist + seg >= startT) {
                        result.push(interpCoord(path[i-1], path[i], (startT - cumDist) / seg));
                        startIns = true;
                    }
                    if (!endIns && cumDist + seg >= endT) {
                        result.push(interpCoord(path[i-1], path[i], (endT - cumDist) / seg));
                        endIns = true;
                    }
                    cumDist += seg;
                    // Early exit: both taper points inserted and we've passed endT
                    if (startIns && endIns && cumDist >= endT) {
                        result.push(path[i]);
                        for (let j = i + 1; j < path.length; j++) result.push(path[j]);
                        return [result, total];
                    }
                }
                result.push(path[i]);
            }
            return [result, total];
        }

        function miterOffset(path, shiftDeg) {
            // Offset in metric space so all segment directions get equal
            // perpendicular distance regardless of lat/lon distortion.
            const n = path.length;
            if (n < 2) return path;
            const latRad = path[0][1] * Math.PI / 180;
            const Mlat = 111320.0;
            const Mlon = Mlat * Math.cos(latRad);
            const shiftM = shiftDeg * Mlat;
            const ox0 = path[0][0], oy0 = path[0][1];

            // Compute metric segment normals in one pass (no intermediate pm[] array)
            const ns = new Array(n - 1);
            let px = 0, py = 0;
            for (let i = 0; i < n - 1; i++) {
                const qx = (path[i+1][0] - ox0) * Mlon;
                const qy = (path[i+1][1] - oy0) * Mlat;
                const dx = qx - px, dy = qy - py;
                const ln = Math.sqrt(dx*dx + dy*dy) || 1;
                ns[i] = [-dy/ln, dx/ln];
                px = qx; py = qy;
            }

            // Apply miter and convert back to degrees
            const result = new Array(n);
            for (let i = 0; i < n; i++) {
                let nx, ny;
                if (i === 0) {
                    [nx, ny] = ns[0];
                } else if (i === n - 1) {
                    [nx, ny] = ns[n - 2];
                } else {
                    const [n1x, n1y] = ns[i - 1];
                    const [n2x, n2y] = ns[i];
                    let bx = n1x + n2x, by = n1y + n2y;
                    const bl = Math.sqrt(bx*bx + by*by);
                    if (bl < 1e-10) {
                        nx = n1x; ny = n1y;
                    } else {
                        bx /= bl; by /= bl;
                        const dot = bx*n1x + by*n1y;
                        const scale = Math.min(Math.abs(dot) > 1e-10 ? 1.0/dot : 4.0, 4.0);
                        nx = bx * scale; ny = by * scale;
                    }
                }
                const pmx = (path[i][0] - ox0) * Mlon + nx * shiftM;
                const pmy = (path[i][1] - oy0) * Mlat + ny * shiftM;
                result[i] = [ox0 + pmx / Mlon, oy0 + pmy / Mlat];
            }
            return result;
        }

        function offsetPath(path, offsetIndex, totalCircuits) {
            if (totalCircuits <= 1) return path;
            const center = (totalCircuits - 1) / 2;
            const shift = (offsetIndex - center) * CIRCUIT_OFFSET_METERS;
            const total = pathLength(path);
            const taper = Math.min(TAPER_DEGREES, total * 0.45);
            const [dense, denseTotal] = insertTaperPoints(path, taper, total);
            const mitered = miterOffset(dense, shift);
            let cumDist = 0;
            return mitered.map((coord, i) => {
                if (i > 0) {
                    const dx = dense[i][0] - dense[i-1][0];
                    const dy = dense[i][1] - dense[i-1][1];
                    cumDist += Math.sqrt(dx*dx + dy*dy);
                }
                const t = taper > 0
                    ? Math.min(cumDist, denseTotal - cumDist, taper) / taper
                    : 1;
                const ox = coord[0] - dense[i][0];
                const oy = coord[1] - dense[i][1];
                return [dense[i][0] + ox*t, dense[i][1] + oy*t];
            });
        }

        function buildCircuitIndex(data) {
            if (!data || !data.length || data[0].circuits === undefined) return null;
            const index = [];
            data.forEach((item, dataIdx) => {
                const circuits = parseInt(item.circuits) || 1;
                for (let i = 0; i < circuits; i++) {
                    index.push({ dataIdx, circuitIdx: i, circuits });
                }
            });
            return index;
        }

        function buildLinesView(data) {
            const index = buildCircuitIndex(data);
            if (!index) return data;
            return index.map(({ dataIdx, circuitIdx, circuits }) => {
                const item = data[dataIdx];
                return {
                    ...item,
                    path: offsetPath(item.path, circuitIdx, circuits),
                };
            });
        }

        function updateCircuitExpansion(zoom) {
            const deck = window.deck;
            if (!deck) return;

            const shouldExpand = zoom >= CIRCUIT_ZOOM_THRESHOLD;
            if (shouldExpand === circuitsExpanded) return;

            circuitsExpanded = shouldExpand;
            const originalData = originalLayerData['Lines'];
            if (!originalData) return;

            deck.setProps({
                layers: deck.props.layers.map(l =>
                    l.id === 'Lines'
                        ? l.clone({ data: shouldExpand ? buildLinesView(originalData) : originalData })
                        : l
                )
            });
        }

        function updateHoveredLine(lineId) {
            if (lineId === hoveredLineId) return;
            hoveredLineId = lineId;

            const deck = window.deck;
            if (!deck) return;

            const hoverColor = isDarkMode ? [255, 255, 255, 255] : [255, 20, 147, 255];

            deck.setProps({
                layers: deck.props.layers.map(l => {
                    if (l.id !== 'Lines') return l;
                    return l.clone({
                        getColor: d => d.line_id === hoveredLineId ? hoverColor : d.color,
                        updateTriggers: { getColor: hoveredLineId }
                    });
                })
            });
        }


        function initializeVoltages() {
            const deck = window.deck;
            if (!deck) return;

            deck.props.layers.forEach(layer => {
                originalLayerData[layer.id] = layer.props.data;
                layer.props.data.forEach(item => {
                    if (item.voltage !== undefined) availableVoltages.add(item.voltage);
                    if (item.voltage_bus0 !== undefined) availableVoltages.add(item.voltage_bus0);
                    if (item.voltage_bus1 !== undefined) availableVoltages.add(item.voltage_bus1);
                });
            });

            availableVoltages = new Set([...availableVoltages].sort((a, b) => b - a));
            updateSuggestions();
        }

        function updateSuggestions(searchTerm = '') {
            const container = document.getElementById('voltage-suggestions');
            const filteredVoltages = [...availableVoltages].filter(v =>
                String(v).includes(searchTerm) && !selectedVoltages.has(v)
            );

            if (filteredVoltages.length === 0) {
                container.innerHTML = `<div style="color:#999;font-size:12px">${searchTerm === '' ? 'No more voltages available' : 'No matches'}</div>`;
                return;
            }

            container.innerHTML = '';
            filteredVoltages.forEach(voltage => {
                const tag = document.createElement('span');
                tag.className = 'voltage-tag';
                tag.textContent = voltage;
                tag.onclick = (e) => {
                    e.stopPropagation();
                    addVoltage(voltage);
                };
                container.appendChild(tag);
            });
        }

        function addVoltage(voltage) {
            selectedVoltages.add(voltage);
            updateSelectedDisplay();
            updateSuggestions();
            applyAllFilters();
            document.getElementById('voltage-filter').value = '';
        }

        function removeVoltage(voltage, event) {
            if (event) event.stopPropagation();
            selectedVoltages.delete(voltage);
            updateSelectedDisplay();
            updateSuggestions();
            applyAllFilters();
        }

        function updateSelectedDisplay() {
            const container = document.getElementById('selected-voltages');
            container.innerHTML = '';

            [...selectedVoltages].sort((a, b) => b - a).forEach(voltage => {
                const tag = document.createElement('span');
                tag.className = 'selected-voltage-tag';
                tag.textContent = voltage;

                const removeBtn = document.createElement('span');
                removeBtn.className = 'remove';
                removeBtn.textContent = '×';
                removeBtn.onclick = (e) => removeVoltage(voltage, e);

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

        function matchesTextSearch(item, searchTerm) {
            if (!searchTerm) return true;

            const orGroups = searchTerm.split('|').map(s => s.trim());
            return orGroups.some(orGroup => {
                const andTerms = orGroup.split('&').map(s => s.trim()).filter(s => s.length > 0);
                return andTerms.every(term => {
                    const lowerTerm = term.toLowerCase();
                    for (const [key, value] of Object.entries(item)) {
                        if (key === 'geometry' || key === 'path' || key === 'poly' || key === 'color') continue;
                        if (String(value).toLowerCase().includes(lowerTerm)) return true;
                    }
                    return false;
                });
            });
        }

        function applyAllFilters() {
            const deck = window.deck;
            if (!deck) return;

            const hasVoltageFilter = selectedVoltages.size > 0;
            const hasTextFilter = currentTextSearch.length > 0;
            const voltages = [...selectedVoltages];

            const layers = deck.props.layers.map(layer => {
                const originalData = originalLayerData[layer.id];

                const filteredData = (hasVoltageFilter || hasTextFilter)
                    ? originalData.filter(item => {
                        let passesVoltageFilter = !hasVoltageFilter;
                        if (hasVoltageFilter) {
                            if (item.voltage !== undefined) {
                                passesVoltageFilter = voltages.includes(item.voltage);
                            } else if (item.voltage_bus0 !== undefined || item.voltage_bus1 !== undefined) {
                                passesVoltageFilter = voltages.includes(item.voltage_bus0) ||
                                    voltages.includes(item.voltage_bus1);
                            } else {
                                passesVoltageFilter = true;
                            }
                        }
                        return passesVoltageFilter && matchesTextSearch(item, currentTextSearch);
                    })
                    : originalData;

                const data = (layer.id === 'Lines' && circuitsExpanded) ? buildLinesView(filteredData) : filteredData;
                return layer.clone({ data });
            });

            deck.setProps({ layers });
        }

        function clearVoltageFilter() {
            selectedVoltages.clear();
            updateSelectedDisplay();
            updateSuggestions();
            applyAllFilters();
            document.getElementById('voltage-filter').value = '';
        }

        function clearTextSearchFilter() {
            currentTextSearch = '';
            document.getElementById('text-search-filter').value = '';
            applyAllFilters();
        }

        document.addEventListener('DOMContentLoaded', function() {
            initializeVoltages();
            applyHashToMap();

            window.addEventListener('hashchange', function() {
                applyHashToMap();
            });

            const voltageInput = document.getElementById('voltage-filter');
            if (voltageInput) {
                voltageInput.addEventListener('input', function(e) {
                    updateSuggestions(e.target.value);
                });
            }

            const textSearchInput = document.getElementById('text-search-filter');
            if (textSearchInput) {
                textSearchInput.addEventListener('input', function(e) {
                    currentTextSearch = e.target.value;
                    applyAllFilters();
                });
            }

            let currentTooltip = null;
            let currentTooltipCoords = null;
            let animationFrameId = null;

            const updateTooltipPosition = () => {
                if (!currentTooltip || !currentTooltipCoords) {
                    if (animationFrameId) {
                        cancelAnimationFrame(animationFrameId);
                        animationFrameId = null;
                    }
                    return;
                }

                const deck = window.deck;
                if (deck && deck.viewManager) {
                    try {
                        const viewport = deck.viewManager.getViewports()[0];
                        if (viewport && viewport.project) {
                            const screenCoords = viewport.project(currentTooltipCoords);
                            if (screenCoords && currentTooltip) {
                                currentTooltip.style.left = (screenCoords[0] + 5) + 'px';
                                currentTooltip.style.top = (screenCoords[1] + 5) + 'px';
                            }
                        }
                    } catch (e) {}
                }

                animationFrameId = requestAnimationFrame(updateTooltipPosition);
            };

            if (window.deck) {
                const originalOnViewStateChange = window.deck.props.onViewStateChange;

                window.deck.setProps({
                    getTooltip: null,
                    pickingRadius: isMobile ? 20 : 10,
                    onHover: ({object}) => {
                        updateHoveredLine(object ? object.line_id : null);
                    },
                    onViewStateChange: ({viewState}) => {
                        if (originalOnViewStateChange) {
                            originalOnViewStateChange({viewState});
                        }
                        updateHash(viewState);
                        updateCircuitExpansion(viewState.zoom);
                        return viewState;
                    }
                });

                const deckContainer = document.getElementById('deck-container');
                if (deckContainer) {
                    let touchStartInfo = null;
                    const TAP_MAX_MOVE = 15;
                    const TAP_MAX_DURATION = 300;

                    if (isMobile) {
                        deckContainer.addEventListener('touchstart', function(event) {
                            if (event.touches.length > 1) {
                                touchStartInfo = null;
                                return;
                            }
                            touchStartInfo = {
                                x: event.touches[0].clientX,
                                y: event.touches[0].clientY,
                                time: Date.now(),
                                wasSingleTouch: true
                            };
                        }, { passive: true });

                        deckContainer.addEventListener('touchmove', function(event) {
                            if (event.touches.length > 1 || !touchStartInfo) {
                                touchStartInfo = null;
                            }
                        }, { passive: true });
                    }

                    const eventType = isMobile ? 'touchend' : 'click';
                    deckContainer.addEventListener(eventType, function(event) {
                        if (isMobile) {
                            if (!touchStartInfo || !touchStartInfo.wasSingleTouch) {
                                touchStartInfo = null;
                                return;
                            }
                            const touch = event.changedTouches ? event.changedTouches[0] : null;
                            if (!touch) return;
                            const dx = touch.clientX - touchStartInfo.x;
                            const dy = touch.clientY - touchStartInfo.y;
                            const dist = Math.sqrt(dx * dx + dy * dy);
                            const duration = Date.now() - touchStartInfo.time;
                            touchStartInfo = null;
                            if (dist > TAP_MAX_MOVE || duration > TAP_MAX_DURATION) return;
                        }

                        const clientX = isMobile && event.changedTouches ? event.changedTouches[0].clientX : event.clientX;
                        const clientY = isMobile && event.changedTouches ? event.changedTouches[0].clientY : event.clientY;

                        if (currentTooltip) {
                            currentTooltip.remove();
                            currentTooltip = null;
                            currentTooltipCoords = null;
                            if (animationFrameId) {
                                cancelAnimationFrame(animationFrameId);
                                animationFrameId = null;
                            }
                        }

                        const deck = window.deck;
                        const pickInfo = deck.pickObject({
                            x: clientX,
                            y: clientY,
                            radius: isMobile ? 20 : 10
                        });

                        if (pickInfo && pickInfo.object) {
                            let html = '<table style="border-collapse:collapse;font-size:12px;line-height:1.4">';
                            let hasContent = false;

                            for (const [key, value] of Object.entries(pickInfo.object)) {
                                if (key === 'geometry' || key === 'path' || key === 'poly' || key === 'color') continue;

                                let displayValue = value;
                                if (key === 'tags' && value) {
                                    const ids = String(value).split(';');
                                    const links = ids.map(id => {
                                        const trimmedId = id.trim();
                                        return `<a href="https://openstreetmap.org/${trimmedId}" target="_blank" rel="noopener noreferrer" style="color:#4a9eff;text-decoration:underline">${trimmedId}</a>`;
                                    }).join('; ');
                                    displayValue = links;
                                }

                                html += `<tr>
                                    <td style="padding:3px 8px 3px 0;font-weight:600;vertical-align:top;color:#aaa;white-space:nowrap">${key}</td>
                                    <td style="padding:3px 0;vertical-align:top;color:#fff;word-break:break-all;max-width:300px">${displayValue}</td>
                                </tr>`;
                                hasContent = true;
                            }
                            html += '</table>';

                            if (!hasContent) return;

                            currentTooltipCoords = pickInfo.coordinate;

                            const tooltip = document.createElement('div');
                            tooltip.innerHTML = html;
                            tooltip.style.cssText = 'position:absolute;background:rgba(0,0,0,0.8);color:white;padding:8px 28px 8px 12px;border-radius:4px;z-index:10000;pointer-events:auto;max-width:400px;box-shadow:0 2px 8px rgba(0,0,0,0.3);font-family:sans-serif;user-select:text;cursor:text;will-change:transform';

                            const closeBtn = document.createElement('div');
                            closeBtn.innerHTML = '×';
                            closeBtn.style.cssText = 'position:absolute;top:4px;right:8px;cursor:pointer;font-size:18px;font-weight:bold;color:#ccc;line-height:1';
                            closeBtn.onclick = function(e) {
                                e.stopPropagation();
                                tooltip.remove();
                                currentTooltip = null;
                                currentTooltipCoords = null;
                                if (animationFrameId) {
                                    cancelAnimationFrame(animationFrameId);
                                    animationFrameId = null;
                                }
                            };
                            closeBtn.onmouseover = () => closeBtn.style.color = '#fff';
                            closeBtn.onmouseout = () => closeBtn.style.color = '#ccc';

                            tooltip.appendChild(closeBtn);
                            currentTooltip = tooltip;
                            document.body.appendChild(tooltip);
                            updateTooltipPosition();
                        }
                    });
                }
            }
        });

        document.getElementById('menu-toggle').addEventListener('click', function(event) {
            event.stopPropagation();
            const menu = document.getElementById('layer-controls');
            menu.style.display = menu.style.display === 'none' ? 'block' : 'none';
        });

        document.getElementById('theme-toggle').addEventListener('click', function(event) {
            event.stopPropagation();
            toggleTheme();
        });

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

    html = html.replace(
        '<div id="deck-container">',
        controls + '<div id="deck-container">',
    )

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

        snakemake = mock_snakemake("prepare_osm_network_release")

    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    # Params and import
    line_types = snakemake.params.line_types
    release_version = snakemake.params.release_version
    n = pypsa.Network(snakemake.input.base_network)
    include_polygons = snakemake.params.include_polygons
    export = snakemake.params.export

    out_buses = snakemake.output.buses if export else ""
    out_lines = snakemake.output.lines if export else ""
    out_links = snakemake.output.links if export else ""
    out_converters = snakemake.output.converters if export else ""
    out_transformers = snakemake.output.transformers if export else ""

    #############
    ### Buses ###
    #############

    logger.info("Cleaning and formatting bus data for release.")
    buses = n.buses.copy()
    buses["dc"] = buses.pop("carrier").map({"DC": "t", "AC": "f"})
    buses.v_nom = buses.v_nom.astype(int)
    buses.sort_index(inplace=True)

    buses = export_clean_csv(buses, BUSES_COLUMNS, out_buses, "bus_id", export)

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

    # Legacy: Old tag column (<= 0.6) contained voltage and circuit info, so we need to clean it up for the release. Not relevant for newer versions.
    lines["tags"] = lines["tags"].apply(
        lambda x: ";".join(set(tag.split("-")[0] for tag in x.split(";")))
    )

    lines = export_clean_csv(lines, LINES_COLUMNS, out_lines, "line_id", export)

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

    links_dc = export_clean_csv(
        links_dc,
        LINKS_COLUMNS,
        out_links,
        "link_id",
        export,
    )

    converters = export_clean_csv(
        converters,
        CONVERTERS_COLUMNS,
        out_converters,
        "converter_id",
        export,
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

    transformers = export_clean_csv(
        transformers,
        TRANSFORMERS_COLUMNS,
        out_transformers,
        "transformer_id",
        export,
    )

    ######################################
    ### Pydeck map with layer controls ###
    ######################################

    buses["geometry"] = buses["geometry"].apply(wkt.loads)
    buses = gpd.GeoDataFrame(buses, geometry="geometry", crs=GEO_CRS)

    if include_polygons:
        # Import stations
        stations_polygon = gpd.read_file(snakemake.input.stations_polygon).to_crs(
            GEO_CRS
        )

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
        buses_polygon = buses_polygon[["tags", "geometry"]]

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

    # Lines PathLayer
    lines_layer = pdk.Layer(
        "PathLayer",
        data=lines.drop(columns=["geometry"]),
        get_path="path",
        get_color="color",
        width_scale=1,
        width_min_pixels=2,
        pickable=True,
        auto_highlight=False,
        parameters={"depthTest": False},
        id="Lines",
    )

    # Links PathLayer
    links_layer = pdk.Layer(
        "PathLayer",
        data=links_dc.drop(columns=["geometry"]),
        get_path="path",
        get_color=[0, 191, 255, 200],
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
        get_elevation=10,
        pickable=True,
        auto_highlight=True,
        parameters={"depthTest": False},
        id="Buses",
    )

    if include_polygons:
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

    # Calculate center point for initial view
    all_coords = [coord for path in lines["path"] for coord in path]
    center_lon = sum(c[0] for c in all_coords) / len(all_coords)
    center_lat = sum(c[1] for c in all_coords) / len(all_coords)

    layers = [
        lines_layer,
        links_layer,
        converters_layer,
        transformers_layer,
        buses_layer,
    ]

    if include_polygons:
        layers = [stations_polygon_layer, buses_polygon_layer] + layers

    # Create the deck with onViewStateChange callback
    map = pdk.Deck(
        layers=layers,
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
