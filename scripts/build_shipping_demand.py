# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build regional demand for international navigation based on outflow volume of
ports.
"""

import json

import geopandas as gpd
import pandas as pd

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_shipping_demand_per_node",
            simpl="",
            clusters=48,
        )

    scope = gpd.read_file(snakemake.input.scope).geometry[0]
    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    demand = pd.read_csv(snakemake.input.demand, index_col=0)[
        "total international navigation"
    ]

    # read port data into GeoDataFrame
    with open(snakemake.input.ports, "r", encoding="latin_1") as f:
        ports = json.load(f)
    ports = pd.json_normalize(ports, "features", sep="_")
    coordinates = ports.geometry_coordinates
    geometry = gpd.points_from_xy(coordinates.str[0], coordinates.str[1])
    ports = gpd.GeoDataFrame(ports, geometry=geometry, crs=4326)

    # filter global port data by European ports
    european_ports = ports[ports.within(scope)]

    # assign ports to nearest region
    p = european_ports.to_crs(3857)
    r = regions.to_crs(3857)
    outflows = (
        p.sjoin_nearest(r).groupby("index_right").properties_outflows.sum().div(1e3)
    )

    # calculate fraction of each country's port outflows
    countries = outflows.index.str[:2]
    outflows_per_country = outflows.groupby(countries).sum()
    fraction = outflows / countries.map(outflows_per_country)

    # distribute per-country demands to nodes based on these fractions
    nodal_demand = demand.loc[countries].fillna(0.0)
    nodal_demand.index = fraction.index
    nodal_demand = nodal_demand.multiply(fraction, axis=0)
    nodal_demand = nodal_demand.reindex(regions.index, fill_value=0)

    # export nodal international navigation demands
    nodal_demand.to_csv(snakemake.output[0])
