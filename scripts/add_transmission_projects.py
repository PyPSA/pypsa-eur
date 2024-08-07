# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Gets the transmission projects defined in the config file and adds them to the
network. Projects are then later included in ` add_electricity.py`.

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
  line_under_construction: keep
  link_under_construction: zero

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`transmission_projects`
Inputs
------

- ``networks/base_network.nc``:  Extract from the geographical vector data of the online `ENTSO-E Interactive Map <https://www.entsoe.eu/data/map/>`_ by the `GridKit <https://github.com/martacki/gridkit>`_ toolkit dating back to March 2022.
- ``data/"project_name"/``: confer :ref:`links`

Outputs
-------

- ``project_lines.csv``

Description
-----------
Gets the transmission projects from e.g. TYNDP or other network development plan and brings them in a network compaitbile formal which then can be imported in `add_electricity`.
"""
import logging
import os

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


def _add_new_buses(n, new_ports, port):
    # Add new buses for the ports which do not have an existing bus close by. If there are multiple ports at the same location, only one bus is added.
    duplicated = new_ports.duplicated(subset=["x", "y"], keep="first")
    to_add = new_ports[~duplicated]
    added_buses = n.madd(
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


def _find_country_for_bus(bus, shapes):
    """
    Find the country of a bus based on its coordinates and the provided
    shapefile.

    Shapefile must contain a column "country" with the country names.
    """
    point = Point(bus.x, bus.y)
    country = shapes.loc[shapes.contains(point), "country"]
    return country.values[0]


def _connect_new_lines(
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
        lines_port = lines.geometry.apply(
            lambda x: pd.Series(
                _get_bus_coords_from_port(x, port=port), index=["x", "y"]
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
                new_port, new_buses = _add_new_buses(n, new_buses, port)
                new_buses.loc[:, "country"] = new_buses.apply(
                    lambda bus: _find_country_for_bus(bus, offshore_shapes), axis=1
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

        lines.loc[lines_port.index, f"bus{port}"] = lines_port.loc[:, "neighbor"]

    lines = lines.assign(under_construction=True)

    return lines, new_buses_df


def _get_branch_coords_from_geometry(linestring, reversed=False):
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


def _get_branch_coords_from_buses(line):
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


def _get_bus_coords_from_port(linestring, port=0):
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


def _find_closest_lines(lines, new_lines, distance_upper_bound=0.1, type="new"):
    """
    Find the closest lines in the existing set of lines to a set of new lines.

    Parameters:
    lines (pandas.DataFrame): DataFrame with column geometry containing the existing lines.
    new_lines (pandas.DataFrame): DataFrame with column geometry containing the new lines.
    distance_upper_bound (float, optional): Maximum distance to consider a line as a match. Defaults to 0.1 which corresponds to approximately 15 km.

    Returns:
    pandas.Series: Series containing with index the new lines and values providing closest existing line.
    """

    # get coordinates of start and end points of all lines, for new lines we need to check both directions
    lines = lines[~lines.geometry.isna()]
    treelines = lines.apply(_get_branch_coords_from_buses, axis=1)
    querylines = pd.concat(
        [
            new_lines.geometry.apply(_get_branch_coords_from_geometry),
            new_lines.geometry.apply(_get_branch_coords_from_geometry, reversed=True),
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
            found = line_map.index
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


def _adjust_decommissioning(upgraded_lines, line_map):
    """
    Adjust the decommissioning year of the existing lines to the built year of
    the upgraded lines.
    """
    to_update = pd.DataFrame(index=line_map)
    to_update.loc[:, "build_year"] = (
        1990  # dummy build_year to make existing lines decommissioned when upgraded lines are built
    )
    to_update.loc[:, "lifetime"] = (
        upgraded_lines.rename(line_map)["build_year"] - 1990
    )  # set lifetime to the difference between build year of upgraded line and existing line
    return to_update


def _get_upgraded_lines(branch_componnent, n, upgraded_lines, line_map):
    """
    Get upgraded lines by merging information of existing line and upgraded
    line.
    """
    # get first the information of the existing lines which will be upgraded
    lines_to_add = n.df(branch_componnent).loc[line_map].copy()
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
    lines_to_add[nominal_attrs[branch_componnent]] = n.df(branch_componnent).loc[
        line_map, nominal_attrs[branch_componnent]
    ]
    # change index of new lines to avoid duplicates
    lines_to_add.index = lines_to_add.index.astype(str) + "_upgraded"
    return lines_to_add


def _get_project_files(plan="tyndp"):
    path = f"data/transmission_projects/{plan}/"
    lines = dict()
    if os.path.exists(path):
        files = os.listdir(path)
    else:
        logger.warning(f"No projects found for {plan}")
        files = []
    if files:
        for file in files:
            if file.endswith(".csv"):
                name = file.split(".")[0]
                df = pd.read_csv(path + file, index_col=0)
                df["geometry"] = df.apply(
                    lambda x: LineString([[x.x0, x.y0], [x.x1, x.y1]]), axis=1
                )
                df.drop(columns=["x0", "y0", "x1", "y1"], inplace=True)
                lines[name] = df
    return lines


def _remove_projects_outside_countries(lines, europe_shape):
    """
    Remove projects which are not in the considered countries.
    """
    europe_shape_prepped = shapely.prepared.prep(europe_shape.buffer(1))
    is_within_covered_countries = lines.geometry.apply(
        lambda x: europe_shape_prepped.contains(x)
    )

    if not is_within_covered_countries.all():
        logger.warning(
            "Project lines outside of the covered area (skipping): "
            + ", ".join(str(i) for i in lines.loc[~is_within_covered_countries].index)
        )

    lines = lines.loc[is_within_covered_countries]
    return lines


def _is_similar(ds1, ds2, percentage=10):
    """
    Check if values in series ds2 are within a specified percentage of series
    ds1.

    Returns:
    - A boolean series where True indicates ds2 values are within the percentage range of ds2.
    """
    lower_bound = ds1 * (1 - percentage / 100)
    upper_bound = ds1 * (1 + percentage / 100)
    return np.logical_and(ds2 >= lower_bound, ds2 <= upper_bound)


def _set_underwater_fraction(new_links, offshore_shapes):
    new_links_gds = gpd.GeoSeries(new_links.geometry)
    new_links.loc[:, "underwater_fraction"] = (
        new_links_gds.intersection(offshore_shapes.union_all()).length
        / new_links_gds.length
    ).round(2)


def _add_projects(
    n,
    new_lines_df,
    new_links_df,
    adjust_lines_df,
    adjust_links_df,
    new_buses_df,
    europe_shape,
    offshore_shapes,
    plan="tyndp",
    status=["confirmed", "under construction"],
):
    lines_dict = _get_project_files(plan)
    logging.info(f"Adding {len(lines_dict)} projects from {plan} to the network.")
    for key, lines in lines_dict.items():
        logging.info(f"Adding {key.replace('_', ' ')} to the network.")
        lines = _remove_projects_outside_countries(lines, europe_shape)
        if isinstance(status, dict):
            status = status[plan]
        lines = lines.loc[lines.project_status.isin(status)]
        if lines.empty:
            continue
        if "new_lines" in key:
            new_lines, new_buses_df = _connect_new_lines(
                lines, n, new_buses_df, bus_carrier="AC"
            )
            duplicate_lines = _find_closest_lines(
                n.lines, new_lines, distance_upper_bound=0.10, type="new"
            )
            # TODO: think about using build_year instead of v_nom
            # ignore duplicates where v_nom is not within a tolerance of 10%
            to_ignore = _is_similar(
                new_lines.loc[duplicate_lines.index, "v_nom"],
                duplicate_lines.map(n.lines["v_nom"]),
            )
            duplicate_lines = duplicate_lines[~to_ignore]
            new_lines = new_lines.drop(duplicate_lines.index, errors="ignore")
            new_lines_df = pd.concat([new_lines_df, new_lines])
            # add new lines to network to be able to find added duplicates
            n.import_components_from_dataframe(new_lines, "Line")
        elif "new_links" in key:
            new_links, new_buses_df = _connect_new_lines(
                lines,
                n,
                new_buses_df,
                offshore_shapes=offshore_shapes,
                distance_upper_bound=0.4,
                bus_carrier=["AC", "DC"],
            )
            duplicate_links = _find_closest_lines(
                n.links, new_links, distance_upper_bound=0.10, type="new"
            )
            # TODO: think about using build_year instead of p_nom
            # ignore duplicates where p_nom is not within a tolerance of 10%
            to_ignore = _is_similar(
                new_links.loc[duplicate_links.index, "p_nom"],
                duplicate_links.map(n.links["p_nom"]),
            )
            duplicate_links = duplicate_links[~to_ignore]
            new_links = new_links.drop(duplicate_links.index, errors="ignore")
            _set_underwater_fraction(new_links, offshore_shapes)
            new_links_df = pd.concat([new_links_df, new_links])
            # add new links to network to be able to find added duplicates
            n.import_components_from_dataframe(new_links, "Link")
        elif "upgraded_lines" in key:
            line_map = _find_closest_lines(
                n.lines, lines, distance_upper_bound=0.30, type="upgraded"
            )
            upgraded_lines = lines.loc[line_map.index]
            lines_to_adjust = _adjust_decommissioning(upgraded_lines, line_map)
            adjust_lines_df = pd.concat([adjust_lines_df, lines_to_adjust])
            upgraded_lines = _get_upgraded_lines("Line", n, upgraded_lines, line_map)
            new_lines_df = pd.concat([new_lines_df, upgraded_lines])
        elif "upgraded_links" in key:
            line_map = _find_closest_lines(
                n.links.query("carrier=='DC'"),
                lines,
                distance_upper_bound=0.30,
                type="upgraded",
            )
            upgraded_links = lines.loc[line_map.index]
            links_to_adjust = _adjust_decommissioning(upgraded_links, line_map)
            adjust_links_df = pd.concat([adjust_links_df, links_to_adjust])
            upgraded_links = _get_upgraded_lines("Link", n, upgraded_links, line_map)
            new_links_df = pd.concat([new_links_df, upgraded_links])
            _set_underwater_fraction(new_links_df, offshore_shapes)
        else:
            logger.warning(f"Unknown project type {key}")
            continue
    return new_lines_df, new_links_df, adjust_lines_df, adjust_links_df, new_buses_df


def fill_length_from_geometry(line, line_factor=1.2):
    length = gpd.GeoSeries(line.geometry, crs=4326).to_crs(3035).length.values[0]
    length = length / 1000 * line_factor
    return round(length, 1)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("add_transmission_projects", run="TYNDP")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    config = snakemake.config
    line_factor = config["lines"]["length_factor"]

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
    for project, include in transmission_projects["include"].items():
        if include:
            (
                new_lines_df,
                new_links_df,
                adjust_lines_df,
                adjust_links_df,
                new_buses_df,
            ) = _add_projects(
                n,
                new_lines_df,
                new_links_df,
                adjust_lines_df,
                adjust_links_df,
                new_buses_df,
                europe_shape,
                offshore_shapes,
                plan=project,
                status=transmission_projects["status"],
            )
    if not new_lines_df.empty:
        line_type = "Al/St 240/40 4-bundle 380.0"
        # Add new line type for new lines
        new_lines_df.loc[:, "type"] = "Al/St 240/40 4-bundle 380.0"
        new_lines_df.loc[:, "num_parallel"] = 2
        (
            new_lines_df["underground"].astype("bool").fillna(False, inplace=True)
            if "underground" in new_lines_df.columns
            else None
        )
        # Add carrier types of lines
        new_lines_df.loc[:, "carrier"] = "AC"
        # Fill empty length values with length calculated from geometry
        new_lines_df["length"] = new_lines_df.apply(
            lambda x: (
                fill_length_from_geometry(x, line_factor)
                if pd.isna(x.length)
                else x.length
            ),
            axis=1,
        )
        # get s_nom from line type
        s_nom = (
            np.sqrt(3)
            * n.line_types.loc[line_type].i_nom
            * new_lines_df["v_nom"]
            * new_lines_df["num_parallel"]
        )
        new_lines_df["s_nom"] = s_nom
    if not new_links_df.empty:
        # Add carrier types of lines and links
        new_links_df.loc[:, "carrier"] = "DC"
        # Fill empty length values with length calculated from geometry
        new_links_df["length"] = new_links_df.apply(
            lambda x: (
                fill_length_from_geometry(x, line_factor)
                if pd.isna(x.length)
                else x.length
            ),
            axis=1,
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
