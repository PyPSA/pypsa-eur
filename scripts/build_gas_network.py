# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Preprocess gas network based on data from bthe SciGRID_gas project
(https://www.gas.scigrid.de/).
"""

import logging

logger = logging.getLogger(__name__)

import geopandas as gpd
import pandas as pd
from pypsa.geo import haversine_pts
from shapely.geometry import Point


def diameter_to_capacity(pipe_diameter_mm):
    """
    Calculate pipe capacity in MW based on diameter in mm.

    20 inch (500 mm)  50 bar -> 1.5   GW CH4 pipe capacity (LHV) 24 inch
    (600 mm)  50 bar -> 5     GW CH4 pipe capacity (LHV) 36 inch (900
    mm)  50 bar -> 11.25 GW CH4 pipe capacity (LHV) 48 inch (1200 mm) 80
    bar -> 21.7  GW CH4 pipe capacity (LHV)

    Based on p.15 of
    https://gasforclimate2050.eu/wp-content/uploads/2020/07/2020_European-Hydrogen-Backbone_Report.pdf
    """
    m1 = (5000 - 1500) / (600 - 500)
    m2 = (11250 - 5000) / (900 - 600)
    a1 = -16000
    a2 = -7500
    if pipe_diameter_mm < 500:
        # slopes definitions
        m0 = (1500 - 0) / (500 - 0)
        # intercept
        a0 = 0
        return a0 + m0 * pipe_diameter_mm
    elif pipe_diameter_mm < 600:
        return a1 + m1 * pipe_diameter_mm
    elif pipe_diameter_mm < 900:
        return a2 + m2 * pipe_diameter_mm
    else:
        m3 = (21700 - 11250) / (1200 - 900)

        a3 = -20100

        return a3 + m3 * pipe_diameter_mm


def load_dataset(fn):
    df = gpd.read_file(fn)
    param = df.param.apply(pd.Series)
    method = df.method.apply(pd.Series)[["diameter_mm", "max_cap_M_m3_per_d"]]
    method.columns = method.columns + "_method"
    df = pd.concat([df, param, method], axis=1)
    to_drop = ["param", "uncertainty", "method", "tags"]
    to_drop = df.columns.intersection(to_drop)
    df.drop(to_drop, axis=1, inplace=True)
    return df


def prepare_dataset(
    df,
    length_factor=1.5,
    correction_threshold_length=4,
    correction_threshold_p_nom=8,
    bidirectional_below=10,
):
    # extract start and end from LineString
    df["point0"] = df.geometry.apply(lambda x: Point(x.coords[0]))
    df["point1"] = df.geometry.apply(lambda x: Point(x.coords[-1]))

    conversion_factor = 437.5  # MCM/day to MWh/h
    df["p_nom"] = df.max_cap_M_m3_per_d * conversion_factor

    # for inferred diameters, assume 500 mm rather than 900 mm (more conservative)
    df.loc[df.diameter_mm_method != "raw", "diameter_mm"] = 500.0

    keep = [
        "name",
        "diameter_mm",
        "is_H_gas",
        "is_bothDirection",
        "length_km",
        "p_nom",
        "max_pressure_bar",
        "start_year",
        "point0",
        "point1",
        "geometry",
    ]
    to_rename = {
        "is_bothDirection": "bidirectional",
        "is_H_gas": "H_gas",
        "start_year": "build_year",
        "length_km": "length",
    }
    df = df[keep].rename(columns=to_rename)

    df.bidirectional = df.bidirectional.astype(bool)
    df.H_gas = df.H_gas.astype(bool)

    # short lines below 10 km are assumed to be bidirectional
    short_lines = df["length"] < bidirectional_below
    df.loc[short_lines, "bidirectional"] = True

    # correct all capacities that deviate correction_threshold factor
    # to diameter-based capacities, unless they are NordStream pipelines
    # also all capacities below 0.5 GW are now diameter-based capacities
    df["p_nom_diameter"] = df.diameter_mm.apply(diameter_to_capacity)
    ratio = df.p_nom / df.p_nom_diameter
    not_nordstream = df.max_pressure_bar < 220
    df.p_nom.update(
        df.p_nom_diameter.where(
            (df.p_nom <= 500)
            | ((ratio > correction_threshold_p_nom) & not_nordstream)
            | ((ratio < 1 / correction_threshold_p_nom) & not_nordstream)
        )
    )

    # lines which have way too discrepant line lengths
    # get assigned haversine length * length factor
    df["length_haversine"] = df.apply(
        lambda p: length_factor
        * haversine_pts([p.point0.x, p.point0.y], [p.point1.x, p.point1.y]),
        axis=1,
    )
    ratio = df.eval("length / length_haversine")
    df["length"].update(
        df.length_haversine.where(
            (df["length"] < 20)
            | (ratio > correction_threshold_length)
            | (ratio < 1 / correction_threshold_length)
        )
    )

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_gas_network")

    logging.basicConfig(level=snakemake.config["logging"]["level"])

    gas_network = load_dataset(snakemake.input.gas_network)

    gas_network = prepare_dataset(gas_network)

    gas_network.to_csv(snakemake.output.cleaned_gas_network)
