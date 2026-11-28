# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Extracts the planned transmission lines and links of the ENTSO-E TYNDP 2026 as project tables for `data/transmission_projects/tyndp2026`.

The investments, their status, commissioning year and element type come from
the workbook of the TYNDP 2026 project platform, and their start and end
points from the routes drawn on its project map. Investments drawn as an area
or a closed loop, which have no fixed route yet, run between the two most
distant points of the drawing. The TYNDP 2026 publishes no technical data, so
voltage, capacity, rating and route length are taken from the TYNDP 2024 through
the investment number, which is stable across editions. The remaining gaps are
filled from the investment description, the transfer capacity of the project
and default values. Links that repeat a project of the German network
development plan are left out. The tables are committed to the repository, so
this rule only needs to run to update them.

References
----------
- [ENTSO-E TYNDP 2026 project platform](https://tyndp2026.entsoe.eu/projects-map/transmission)
"""

import logging
import re

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import shapely
from scipy.spatial.distance import pdist, squareform

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

PROJECT_URL = "https://tyndp2026.entsoe.eu/projects-map/Transmission/"
DISTANCE_CRS = 3035
# share of the project length by which the end points of duplicates may differ
MAX_OFFSET = 0.25

COLUMNS = {
    "This investment belongs to project number…": "project",
    "Investment Name": "name",
    "Description of the investment": "description",
    "Commissioning Year": "commissioning",
    "Status": "status",
    "Type of Element": "element",
    "Total route length (km)": "length",
    "Nominal voltage (kV)": "v_nom",
    "Capacity of the investment (MW)": "p_nom",
    "MVA": "s_nom",
    "Type of project - Is it an upgrade of existing infrastructure?": "upgrade",
    "Transfer capacity increase A-B (MW)": "transfer_capacity",
}
TECHNICAL = ["length", "v_nom", "p_nom", "s_nom"]
STATUS = {
    "Under Consideration": "under_consideration",
    "Planned But Not Yet Permitting": "planned_not_yet_permitted",
    "In Permitting": "in_permitting",
    "Under Construction": "under_construction",
    "Commissioned": "completed",
    "Completed": "completed",
}
SECOND_CIRCUIT = r"(?:second|2nd) circuit"
UPGRADE = (
    rf"upgrad|uprat|reconduct|reconstruct|refurbish|moderni[sz]|HTLS|{SECOND_CIRCUIT}"
)
REPLACEMENT = r"repl[ea]+c|end of (?:their|its) technical life"
TWO_CIRCUITS = rf"double|two|2\s?x|{SECOND_CIRCUIT}"
# new 380 kV lines that replace a 220 kV line according to their promoters
VOLTAGE_UPGRADES = [630, 1473, 1768, 2001]
# investments that are entered by hand in `data/transmission_projects/manual`
MANUAL = [1458]

OUTPUT = ["project_status", "length", "build_year", "underground", "tags"]
LINE_OUTPUT = ["v_nom", "num_parallel", "type"]
LINK_OUTPUT = ["p_nom"]
COORDS = ["x0", "y0", "x1", "y1"]


def read_sheet(path: str, sheet: str) -> pd.DataFrame:
    """Read a sheet of a TYNDP project portfolio with short column names."""
    df = pd.read_excel(path, sheet_name=sheet, header=1, index_col=0)
    return df.rename(columns=lambda c: COLUMNS.get(c.strip(), c))


def read_investments(portfolio: str, technical: str) -> pd.DataFrame:
    """Read the transmission lines and cables with technical data of an older portfolio."""
    df = read_sheet(portfolio, "Trans.Investments")
    projects = read_sheet(portfolio, "Trans.Projects")[["upgrade", "transfer_capacity"]]
    technical_data = read_sheet(technical, "Trans.Investments")[TECHNICAL]
    df = df.join(projects, on="project").join(
        technical_data.apply(pd.to_numeric, errors="coerce").replace(0, np.nan)
    )
    is_branch = df.element.str.contains("TransmissionLine|TransmissionCable")
    return df[is_branch].drop(MANUAL).sort_index()


def endpoints(route: shapely.Geometry) -> shapely.LineString | None:
    """Reduce a drawn route to a straight line between its two end points."""
    if route.geom_type == "Point":
        return None
    coords = shapely.get_coordinates(route)
    ends = coords[[0, -1]]
    is_open_line = route.geom_type == "LineString" and (
        np.linalg.norm(ends[0] - ends[1]) > 0.1 * route.length
    )
    if not is_open_line:
        coords = shapely.get_coordinates(route.convex_hull)
        distances = squareform(pdist(coords))
        ends = coords[list(np.unravel_index(distances.argmax(), distances.shape))]
    return shapely.LineString(ends)


def read_routes(path: str) -> gpd.GeoSeries:
    """Read the drawn routes of the investments as straight lines in metres."""
    routes = gpd.read_file(path).query("type == 'Transmission'")
    routes = routes.set_index("investmentId").geometry
    is_valid = routes.bounds.eval("-180 < minx and maxx < 180")
    lines = routes[is_valid].to_crs(DISTANCE_CRS).map(endpoints).dropna()
    return gpd.GeoSeries(lines, crs=DISTANCE_CRS)


def hstore(tags: pd.DataFrame) -> pd.Series:
    """Format columns as key-value tags like those of the base network."""
    tags = tags.astype(str).apply(lambda col: col.str.replace('"', "'"))
    return tags.apply(
        lambda row: ", ".join(f'"{k}"=>"{v}"' for k, v in row.items()), axis=1
    )


def remove_duplicates(df: pd.DataFrame, other: pd.DataFrame) -> pd.DataFrame:
    """Remove the projects of which both end points lie close to those of another project."""
    ends = df[COORDS].values.reshape(-1, 1, 2, 2)
    other_ends = other[COORDS].values.reshape(1, -1, 2, 2)
    offsets = [
        np.linalg.norm(ends - other_ends[:, :, order], axis=3).max(axis=2)
        for order in ([0, 1], [1, 0])
    ]
    length = np.linalg.norm(ends[:, :, 0] - ends[:, :, 1], axis=2)
    is_duplicate = (np.minimum(*offsets) < MAX_OFFSET * length).any(axis=1)
    logger.info(f"Removing duplicates: {', '.join(df.index[is_duplicate])}")
    return df[~is_duplicate]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("extract_tyndp_projects")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    df = read_investments(snakemake.input.portfolio, snakemake.input.technical)
    routes = read_routes(snakemake.input.routes)
    routes = routes.reindex(df.index).dropna()
    df = df.loc[routes.index]
    df[COORDS] = [np.ravel(route.coords) for route in routes.to_crs(4326)]
    df[COORDS] = df[COORDS].round(4)

    text = df.name + ". " + df.description.fillna("")
    distance = routes.length / 1e3
    df["length"] = df.length.where(df.length.between(distance, 2 * distance))
    df["project_status"] = df.status.map(STATUS)
    df["build_year"] = df.commissioning.str[-4:].astype(int)
    df["underground"] = df.element.str.contains("Cable")
    tags = {
        "url": PROJECT_URL + df.project.astype(str),
        "name": df.name.str.strip(),
        "tyndp_project": df.project,
        "tyndp_investment": df.index,
    }
    df["tags"] = hstore(pd.DataFrame(tags))

    voltages = np.array(snakemake.params.voltages)
    described_voltage = text.str.extract(r"(\d{3})\s?-?kV", flags=re.IGNORECASE)[0]
    v_nom = df.v_nom.fillna(described_voltage.astype(float))
    is_modelled = ~(v_nom < 0.9 * voltages.min())
    nearest = np.abs(v_nom.fillna(380).values[:, None] - voltages).argmin(axis=1)
    df["v_nom"] = voltages[nearest]
    df["type"] = df.v_nom.map(snakemake.params.line_types)
    i_nom = df.type.map(pypsa.Network().line_types.i_nom)
    rated_circuits = df.s_nom.div(np.sqrt(3) * df.v_nom * i_nom).round().fillna(1)
    described_circuits = text.str.contains(TWO_CIRCUITS, case=False) + 1
    df["num_parallel"] = np.maximum(rated_circuits, described_circuits).astype(int)

    # largest transfer capacity of the project across its borders
    numbers = df.transfer_capacity.astype(str).str.extractall(r"(\d+\.?\d*)")[0]
    df["p_nom"] = df.p_nom.fillna(numbers.astype(float).groupby(level=0).max())

    is_upgrade = (
        df.upgrade.eq("Upgrade")
        | df.upgrade.str.startswith("Mostly") & text.str.contains(UPGRADE, case=False)
        | text.str.contains(REPLACEMENT, case=False)
        | df.index.isin(VOLTAGE_UPGRADES)
    )
    is_dc = df.element.str.contains("DC")
    is_line = ~is_dc & is_modelled
    df.index = ("TYNDP2026_" + df.index.astype(str)).rename(None)
    tables = {
        "new_lines": df[is_line.values & ~is_upgrade.values],
        "upgraded_lines": df[is_line.values & is_upgrade.values],
        "new_links": df[is_dc.values & ~is_upgrade.values],
        "upgraded_links": df[is_dc.values & is_upgrade.values],
    }
    for kind in ["lines", "links"]:
        nep = pd.read_csv(snakemake.input[f"nep_{kind}"], index_col=0)
        tables[f"new_{kind}"] = remove_duplicates(tables[f"new_{kind}"], nep)
    for name, table in tables.items():
        columns = LINK_OUTPUT if "links" in name else LINE_OUTPUT
        table[OUTPUT + columns + COORDS].to_csv(snakemake.output[name])
        logger.info(f"Extracted {len(table)} {name.replace('_', ' ')}.")
