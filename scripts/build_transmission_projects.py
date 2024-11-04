# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Gets the transmission projects defined in the config file, concatenates and
deduplicates them. Projects are later included in :mod:`add_electricity.py`.

Relevant Settings
-----------------

.. code:: yaml

transmission_projects:
  include:
    #tyndp: true # For later, when other TYNDP projects are combined with new version
    nep: true
  status:
  - confirmed
  - in_permitting
  - under_construction
    #- under_consideration
  link_under_construction: zero

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`transmission_projects`
Inputs
------

- ``networks/base_network.nc``:  Base network topology for the electricity grid. This is processed in :mod:`base_network.py`.
- ``data/transmission_projects/"project_name"/``: Takes the transmission projects from the subfolder of data/transmission_projects. The subfolder name is the project name.
- ``offshore_shapes.geojson``: Shapefile containing the offshore regions. Used to determine if a new bus should be added for a new line or link.
- ``europe_shape.geojson``: Shapefile containing the shape of Europe. Used to determine if a project is within the considered countries.

Outputs
-------

- ``transmission_projects/new_lines.csv``: New project lines to be added to the network. This includes new lines and upgraded lines.
- ``transmission_projects/new_links.csv``: New project links to be added to the network. This includes new links and upgraded links.
- ``transmission_projects/adjust_lines.csv``: For lines which are upgraded, the decommissioning year of the existing line is adjusted to the build year of the upgraded line.
- ``transmission_projects/adjust_links.csv``: For links which are upgraded, the decommissioning year of the existing link is adjusted to the build year of the upgraded link.
- ``transmission_projects/new_buses.csv``: For some links, we have to add new buses (e.g. North Sea Wind Power Hub).
"""

import logging
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import shapely
from _helpers import configure_logging, set_scenario_config
from pypsa.descriptors import nominal_attrs
from scipy import spatial
from shapely.geometry import LineString, Point

logger = logging.getLogger(__name__)


def add_new_buses(n, new_ports):
    # Add new buses for the ports which do not have an existing bus close by. If there are multiple ports at the same location, only one bus is added.
    duplicated = new_ports.duplicated(subset=["x", "y"], keep="first")
    to_add = new_ports[~duplicated]
    added_buses = n.add(
        "Bus",
        names=to_add.index,
        suffix=" bus",
        x=to_add.x,
        y=to_add.y,
        v_nom=380,
        under_construction=True,
        symbol="substation",
        substation_off=True,
        substation_lv=False,
        carrier="AC",
    )
    new_buses = n.buses.loc[added_buses].copy().dropna(axis=1, how="all")
    new_ports.loc[to_add.index, "neighbor"] = added_buses
    new_ports["neighbor"] = new_ports.groupby(["x", "y"])["neighbor"].transform("first")
    return new_ports, new_buses


def find_country_for_bus(bus, shapes):
    """
    Find the country of a bus based on its coordinates and the provided
    shapefile.

    Shapefile must contain a column "country" with the country names.
    """
    point = Point(bus.x, bus.y)
    country = shapes.loc[shapes.contains(point), "country"]
    return country.values[0]


def connect_new_lines(
    lines,
    n,
    new_buses_df,
    offshore_shapes=None,
    distance_upper_bound=np.inf,
    bus_carrier="AC",
):
    """
    Find the closest existing bus to the port of each line.

    If closest bus is further away than distance_upper_bound and is
    inside an offshore region, a new bus is created. and the line is
    connected to it.
    """
    bus_carrier = np.atleast_1d(bus_carrier)
    buses = n.buses.query("carrier in @bus_carrier").copy()
    bus_tree = spatial.KDTree(buses[["x", "y"]])

    for port in [0, 1]:
        lines_port = lines["geometry"].apply(
            lambda x: pd.Series(
                get_bus_coords_from_port(x, port=port), index=["x", "y"]
            )
        )
        distances, indices = bus_tree.query(lines_port)
        # Series of lines with closest bus in the existing network and whether they match the distance criterion
        lines_port["neighbor"] = buses.iloc[indices].index
        lines_port["match_distance"] = distances < distance_upper_bound
        # For buses which are not close to any existing bus, only add a new bus if the line is going offshore (e.g. North Sea Wind Power Hub)
        if not lines_port.match_distance.all() and offshore_shapes.union_all():
            potential_new_buses = lines_port[~lines_port.match_distance]
            is_offshore = potential_new_buses.apply(
                lambda x: offshore_shapes.union_all().contains(Point(x.x, x.y)), axis=1
            )
            new_buses = potential_new_buses[is_offshore]
            if not new_buses.empty:
                new_port, new_buses = add_new_buses(n, new_buses)
                new_buses["country"] = new_buses.apply(
                    lambda bus: find_country_for_bus(bus, offshore_shapes), axis=1
                )
                lines_port.loc[new_port.index, "match_distance"] = True
                lines_port.loc[new_port.index, "neighbor"] = new_port.neighbor
                new_buses_df = pd.concat([new_buses_df, new_buses])

        if not lines_port.match_distance.all():
            logging.warning(
                "Could not find bus close enough to connect the the following lines:\n"
                + str(lines_port[~lines_port.match_distance].index.to_list())
                + "\n Lines will be ignored."
            )
            lines.drop(lines_port[~lines_port.match_distance].index, inplace=True)
            lines_port = lines_port[lines_port.match_distance]

        lines.loc[lines_port.index, f"bus{port}"] = lines_port["neighbor"]

    lines = lines.assign(under_construction=True)

    return lines, new_buses_df


def get_branch_coords_from_geometry(linestring, reversed=False):
    """
    Reduces a linestring to its start and end points. Used to simplify the
    linestring which can have more than two points.

    Parameters:
    linestring: Shapely linestring
    reversed (bool, optional): If True, returns the end and start points instead of the start and end points.
                               Defaults to False.

    Returns:
    numpy.ndarray: Flattened array of start and end coordinates.
    """
    coords = np.asarray(linestring.coords)
    ind = [0, -1] if not reversed else [-1, 0]
    start_end_coords = coords[ind]
    return start_end_coords.flatten()


def get_branch_coords_from_buses(line):
    """
    Gets line string for branch component in an pypsa network.

    Parameters:
    linestring: shapely linestring
    reversed (bool, optional): If True, returns the end and start points instead of the start and end points.
                               Defaults to False.

    Returns:
    numpy.ndarray: Flattened array of start and end coordinates.
    """
    start_coords = n.buses.loc[line.bus0, ["x", "y"]].values
    end_coords = n.buses.loc[line.bus1, ["x", "y"]].values
    return np.array([start_coords, end_coords]).flatten()


def get_bus_coords_from_port(linestring, port=0):
    """
    Extracts the coordinates of a specified port from a given linestring.

    Parameters:
    linestring: The shapely linestring.
    port (int): The index of the port to extract coordinates from. Default is 0.

    Returns:
    tuple: The coordinates of the specified port as a tuple (x, y).
    """
    coords = np.asarray(linestring.coords)
    ind = [0, -1]
    coords = coords[ind]
    coords = coords[port]
    return coords


def find_closest_lines(lines, new_lines, distance_upper_bound=0.1, type="new"):
    """
    Find the closest lines in the existing set of lines to a set of new lines.

    Parameters:
    lines (pandas.DataFrame): DataFrame of the existing lines.
    new_lines (pandas.DataFrame): DataFrame with column geometry containing the new lines.
    distance_upper_bound (float, optional): Maximum distance to consider a line as a match. Defaults to 0.1 which corresponds to approximately 15 km.

    Returns:
    pandas.Series: Series containing with index the new lines and values providing closest existing line.
    """

    # get coordinates of start and end points of all lines, for new lines we need to check both directions
    treelines = lines.apply(get_branch_coords_from_buses, axis=1)
    querylines = pd.concat(
        [
            new_lines["geometry"].apply(get_branch_coords_from_geometry),
            new_lines["geometry"].apply(get_branch_coords_from_geometry, reversed=True),
        ]
    )
    treelines = np.vstack(treelines)
    querylines = np.vstack(querylines)
    tree = spatial.KDTree(treelines)
    dist, ind = tree.query(querylines, distance_upper_bound=distance_upper_bound)
    found_b = ind < len(lines)
    # since the new lines are checked in both directions, we need to find the correct index of the new line
    found_i = np.arange(len(querylines))[found_b] % len(new_lines)
    # create a DataFrame with the distances, new line and its closest existing line
    line_map = pd.DataFrame(
        dict(D=dist[found_b], existing_line=lines.index[ind[found_b] % len(lines)]),
        index=new_lines.index[found_i].rename("new_lines"),
    )
    if type == "new":
        if len(found_i) != 0:
            # compare if attribute of new line and existing line is similar
            attr = "p_nom" if "p_nom" in lines else "v_nom"
            # potential duplicates
            duplicated = line_map["existing_line"]
            # only if lines are similar in terms of p_nom or v_nom they are kept as duplicates
            to_keep = is_similar(
                new_lines.loc[duplicated.index, attr],
                duplicated.map(lines[attr]),
                percentage=10,
            )
            line_map = line_map[to_keep]
            if not line_map.empty:
                logger.warning(
                    "Found new lines similar to existing lines:\n"
                    + str(line_map["existing_line"].to_dict())
                    + "\n Lines are assumed to be duplicated and will be ignored."
                )
    elif type == "upgraded":
        if len(found_i) < len(new_lines):
            not_found = new_lines.index.difference(line_map.index)
            logger.warning(
                "Could not find upgraded lines close enough to existing lines:\n"
                + str(not_found.to_list())
                + "\n Lines will be ignored."
            )
    # only keep the closer line of the new line pair (since lines are checked in both directions)
    line_map = line_map.sort_values(by="D")[
        lambda ds: ~ds.index.duplicated(keep="first")
    ].sort_index()["existing_line"]
    return line_map


def adjust_decommissioning(upgraded_lines, line_map):
    """
    Adjust the decommissioning year of the existing lines to the built year of
    the upgraded lines.
    """
    to_update = pd.DataFrame(index=line_map)
    to_update["build_year"] = (
        1990  # dummy build_year to make existing lines decommissioned when upgraded lines are built
    )
    to_update["lifetime"] = (
        upgraded_lines.rename(line_map)["build_year"] - 1990
    )  # set lifetime to the difference between build year of upgraded line and existing line
    return to_update


def get_upgraded_lines(branch_component, n, upgraded_lines, line_map):
    """
    Get upgraded lines by merging information of existing line and upgraded
    line.
    """
    # get first the information of the existing lines which will be upgraded
    lines_to_add = n.df(branch_component).loc[line_map].copy()
    # get columns of upgraded lines which are not in existing lines
    new_columns = upgraded_lines.columns.difference(lines_to_add.columns)
    # rename upgraded lines to match existing lines
    upgraded_lines = upgraded_lines.rename(line_map)
    # set the same index names to be able to merge
    upgraded_lines.index.name = lines_to_add.index.name
    # merge upgraded lines with existing lines
    lines_to_add.update(upgraded_lines)
    # add column which was added in upgraded lines
    lines_to_add = pd.concat([lines_to_add, upgraded_lines[new_columns]], axis=1)
    # only consider columns of original upgraded lines and bus0 and bus1
    lines_to_add = lines_to_add.loc[:, ["bus0", "bus1", *upgraded_lines.columns]]
    # set capacity of upgraded lines to capacity of existing lines
    lines_to_add[nominal_attrs[branch_component]] = n.df(branch_component).loc[
        line_map, nominal_attrs[branch_component]
    ]
    # change index of new lines to avoid duplicates
    lines_to_add.index = lines_to_add.index.astype(str) + "_upgraded"
    return lines_to_add


def get_project_files(path, skip=[]):
    path = Path(path)
    lines = {}
    files = [
        p
        for p in path.iterdir()
        if p.is_file()
        and p.suffix == ".csv"
        and not any(substring in p.name for substring in skip)
    ]
    if not files:
        logger.warning(f"No projects found for {path.parent.name}")
        return lines
    for file in files:
        df = pd.read_csv(file, index_col=0)
        df["geometry"] = df.apply(
            lambda x: LineString([[x.x0, x.y0], [x.x1, x.y1]]), axis=1
        )
        df.drop(columns=["x0", "y0", "x1", "y1"], inplace=True)
        lines[file.stem] = df
    return lines


def remove_projects_outside_countries(lines, europe_shape):
    """
    Remove projects which are not in the considered countries.
    """
    europe_shape_prepped = shapely.prepared.prep(europe_shape)
    is_within_covered_countries = lines["geometry"].apply(
        lambda x: europe_shape_prepped.contains(x)
    )

    if not is_within_covered_countries.all():
        logger.warning(
            "Project lines outside of the covered area (skipping): "
            + ", ".join(str(i) for i in lines.loc[~is_within_covered_countries].index)
        )

    lines = lines.loc[is_within_covered_countries]
    return lines


def is_similar(ds1, ds2, percentage=10):
    """
    Check if values in series ds2 are within a specified percentage of series
    ds1.

    Returns:
    - A boolean series where True indicates ds2 values are within the percentage range of ds2.
    """
    lower_bound = ds1 * (1 - percentage / 100)
    upper_bound = ds1 * (1 + percentage / 100)
    return np.logical_and(ds2 >= lower_bound, ds2 <= upper_bound)


def set_underwater_fraction(new_links, offshore_shapes):
    new_links_gds = gpd.GeoSeries(new_links["geometry"])
    new_links["underwater_fraction"] = (
        new_links_gds.intersection(offshore_shapes.union_all()).length
        / new_links_gds.length
    ).round(2)


def add_projects(
    n,
    new_lines_df,
    new_links_df,
    adjust_lines_df,
    adjust_links_df,
    new_buses_df,
    europe_shape,
    offshore_shapes,
    path,
    plan,
    status=["confirmed", "under construction"],
    skip=[],
):
    lines_dict = get_project_files(path, skip=skip)
    for key, lines in lines_dict.items():
        logging.info(f"Processing {key.replace('_', ' ')} projects from {plan}.")
        lines = remove_projects_outside_countries(lines, europe_shape)
        if isinstance(status, dict):
            status = status[plan]
        lines = lines.loc[lines.project_status.isin(status)]
        if lines.empty:
            continue
        if key == "new_lines":
            new_lines, new_buses_df = connect_new_lines(
                lines, n, new_buses_df, bus_carrier="AC"
            )
            duplicate_lines = find_closest_lines(
                n.lines, new_lines, distance_upper_bound=0.10, type="new"
            )
            new_lines = new_lines.drop(duplicate_lines.index, errors="ignore")
            new_lines_df = pd.concat([new_lines_df, new_lines])
            # add new lines to network to be able to find added duplicates
            n.add("Line", new_lines.index, **new_lines)
        elif key == "new_links":
            new_links, new_buses_df = connect_new_lines(
                lines,
                n,
                new_buses_df,
                offshore_shapes=offshore_shapes,
                distance_upper_bound=0.4,
                bus_carrier=["AC", "DC"],
            )
            duplicate_links = find_closest_lines(
                n.links, new_links, distance_upper_bound=0.10, type="new"
            )
            new_links = new_links.drop(duplicate_links.index, errors="ignore")
            set_underwater_fraction(new_links, offshore_shapes)
            new_links_df = pd.concat([new_links_df, new_links])
            # add new links to network to be able to find added duplicates
            n.add("Link", new_links.index, **new_links)
        elif key == "upgraded_lines":
            line_map = find_closest_lines(
                n.lines, lines, distance_upper_bound=0.30, type="upgraded"
            )
            upgraded_lines = lines.loc[line_map.index]
            lines_to_adjust = adjust_decommissioning(upgraded_lines, line_map)
            adjust_lines_df = pd.concat([adjust_lines_df, lines_to_adjust])
            upgraded_lines = get_upgraded_lines("Line", n, upgraded_lines, line_map)
            new_lines_df = pd.concat([new_lines_df, upgraded_lines])
        elif key == "upgraded_links":
            line_map = find_closest_lines(
                n.links.query("carrier=='DC'"),
                lines,
                distance_upper_bound=0.30,
                type="upgraded",
            )
            upgraded_links = lines.loc[line_map.index]
            links_to_adjust = adjust_decommissioning(upgraded_links, line_map)
            adjust_links_df = pd.concat([adjust_links_df, links_to_adjust])
            upgraded_links = get_upgraded_lines("Link", n, upgraded_links, line_map)
            new_links_df = pd.concat([new_links_df, upgraded_links])
            set_underwater_fraction(new_links_df, offshore_shapes)
        else:
            logger.warning(f"Unknown project type {key}")
            continue
    return new_lines_df, new_links_df, adjust_lines_df, adjust_links_df, new_buses_df


def fill_length_from_geometry(line, line_factor=1.2):
    if not pd.isna(line.length):
        return line.length
    length = gpd.GeoSeries(line["geometry"], crs=4326).to_crs(3035).length.values[0]
    length = length / 1000 * line_factor
    return round(length, 1)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_transmission_projects")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    line_factor = snakemake.params.line_factor
    s_max_pu = snakemake.params.s_max_pu

    n = pypsa.Network(snakemake.input.base_network)

    new_lines_df = pd.DataFrame()
    new_links_df = pd.DataFrame()
    adjust_lines_df = pd.DataFrame()
    adjust_links_df = pd.DataFrame()
    new_buses_df = pd.DataFrame()

    europe_shape = gpd.read_file(snakemake.input.europe_shape).loc[0, "geometry"]
    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes).rename(
        {"name": "country"}, axis=1
    )

    transmission_projects = snakemake.params.transmission_projects
    projects = [
        project
        for project, include in transmission_projects["include"].items()
        if include
    ]
    paths = snakemake.input.transmission_projects
    for project in projects:
        path = list(filter(lambda path: project in path, paths))[0]
        new_lines_df, new_links_df, adjust_lines_df, adjust_links_df, new_buses_df = (
            add_projects(
                n,
                new_lines_df,
                new_links_df,
                adjust_lines_df,
                adjust_links_df,
                new_buses_df,
                europe_shape,
                offshore_shapes,
                path=path,
                plan=project,
                status=transmission_projects["status"],
                skip=transmission_projects["skip"],
            )
        )
    if not new_lines_df.empty:
        line_type = "Al/St 240/40 4-bundle 380.0"
        # Add new line type for new lines
        new_lines_df["type"] = new_lines_df["type"].fillna(line_type)
        new_lines_df["num_parallel"] = new_lines_df["num_parallel"].fillna(2)
        if "underground" in new_lines_df.columns:
            new_lines_df["underground"] = (
                new_lines_df["underground"].astype("bool").fillna(False)
            )
        # Add carrier types of lines
        new_lines_df["carrier"] = "AC"
        # Fill empty length values with length calculated from geometry
        new_lines_df["length"] = new_lines_df.apply(
            fill_length_from_geometry, args=(line_factor,), axis=1
        )
        # get s_nom from line type
        new_lines_df["s_nom"] = (
            np.sqrt(3)
            * n.line_types.loc[new_lines_df["type"], "i_nom"].values
            * new_lines_df["v_nom"]
            * new_lines_df["num_parallel"]
        ).round(2)
        # set s_max_pu
        new_lines_df["s_max_pu"] = s_max_pu
    if not new_links_df.empty:
        # Add carrier types of lines and links
        new_links_df["carrier"] = "DC"
        # Fill empty length values with length calculated from geometry
        new_links_df["length"] = new_links_df.apply(
            fill_length_from_geometry, args=(line_factor,), axis=1
        )
        # Whether to keep existing link capacity or set to zero
        not_upgraded = ~new_links_df.index.str.contains("upgraded")
        if transmission_projects["new_link_capacity"] == "keep":
            new_links_df.loc[not_upgraded, "p_nom"] = new_links_df["p_nom"].fillna(0)
        elif transmission_projects["new_link_capacity"] == "zero":
            new_links_df.loc[not_upgraded, "p_nom"] = 0
    # export csv files for new buses, lines, links and adjusted lines and links
    new_lines_df.to_csv(snakemake.output.new_lines)
    new_links_df.to_csv(snakemake.output.new_links)
    adjust_lines_df.to_csv(snakemake.output.adjust_lines)
    adjust_links_df.to_csv(snakemake.output.adjust_links)
    new_buses_df.to_csv(snakemake.output.new_buses)
