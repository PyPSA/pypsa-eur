# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 @LukasFranken, The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This rule extracts potential and cost for electricity generation through
enhanced geothermal systems.

For this, we use data from "From hot rock to useful energy..." by Aghahosseini, Breyer (2020)
'https://www.sciencedirect.com/science/article/pii/S0306261920312551'
Note that we input data used here is not the same as in the paper, but was passed on by the authors.

The data provides a lon-lat gridded map of Europe (1° x 1°), with each grid cell assigned
a heat potential (in GWh) and a cost (in EUR/MW).

This scripts overlays that map with the network's regions, and builds a csv with CAPEX, OPEX and p_nom_max
"""

import logging

logger = logging.getLogger(__name__)

import json

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Polygon


def prepare_egs_data(egs_file):
    with open(egs_file) as f:
        jsondata = json.load(f)

    def point_to_square(p, lon_extent=1.0, lat_extent=1.0):
        try:
            x, y = p.coords.xy[0][0], p.coords.xy[1][0]
        except IndexError:
            return p

        return Polygon(
            [
                [x - lon_extent / 2, y - lat_extent / 2],
                [x - lon_extent / 2, y + lat_extent / 2],
                [x + lon_extent / 2, y + lat_extent / 2],
                [x + lon_extent / 2, y - lat_extent / 2],
            ]
        )

    years = [2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050]
    lcoes = ["LCOE50", "LCOE100", "LCOE150"]

    egs_data = dict()

    for year in years:
        df = pd.DataFrame(columns=["Lon", "Lat", "CAPEX", "HeatSust", "PowerSust"])

        for lcoe in lcoes:
            for country_data in jsondata[lcoe]:
                try:
                    country_df = pd.DataFrame(
                        columns=df.columns,
                        index=range(len(country_data[0][years.index(year)]["Lon"])),
                    )
                except TypeError:
                    country_df = pd.DataFrame(columns=df.columns, index=range(0))

                for col in df.columns:
                    country_df[col] = country_data[0][years.index(year)][col]

                if country_df.dropna().empty:
                    continue
                elif df.empty:
                    df = country_df.dropna()
                else:
                    df = pd.concat((df, country_df.dropna()), ignore_index=True)

        gdf = gpd.GeoDataFrame(
            df.drop(columns=["Lon", "Lat"]), geometry=gpd.points_from_xy(df.Lon, df.Lat)
        ).reset_index(drop=True)

        gdf["geometry"] = gdf.geometry.apply(lambda geom: point_to_square(geom))
        egs_data[year] = gdf

    return egs_data


def get_capacity_factors(network_regions_file, air_temperatures_file):
    """
    Performance of EGS is higher for lower temperatures, due to more efficient
    air cooling Data from Ricks et al.: The Role of Flexible Geothermal Power
    in Decarbonized Elec Systems.
    """

    # these values are taken from the paper's
    # Supplementary Figure 20 from https://zenodo.org/records/7093330
    # and relate deviations of the ambient temperature from the year-average
    # ambient temperature to EGS capacity factors.
    delta_T = [-15, -10, -5, 0, 5, 10, 15, 20]
    cf = [1.17, 1.13, 1.07, 1, 0.925, 0.84, 0.75, 0.65]

    x = np.linspace(-15, 20, 200)
    y = np.interp(x, delta_T, cf)

    upper_x = np.linspace(20, 25, 50)
    m_upper = (y[-1] - y[-2]) / (x[-1] - x[-2])
    upper_y = upper_x * m_upper - x[-1] * m_upper + y[-1]

    lower_x = np.linspace(-20, -15, 50)
    m_lower = (y[1] - y[0]) / (x[1] - x[0])
    lower_y = lower_x * m_lower - x[0] * m_lower + y[0]

    x = np.hstack((lower_x, x, upper_x))
    y = np.hstack((lower_y, y, upper_y))

    network_regions = gpd.read_file(network_regions_file).set_crs(epsg=4326)
    index = network_regions["name"]

    air_temp = xr.open_dataset(air_temperatures_file)

    snapshots = pd.date_range(freq="h", **snakemake.params.snapshots)
    capacity_factors = pd.DataFrame(index=snapshots)

    # bespoke computation of capacity factors for each bus.
    # Considering the respective temperatures, we compute
    # the deviation from the average temperature and relate it
    # to capacity factors based on the data from above.
    for bus in index:
        temp = air_temp.sel(name=bus).to_dataframe()["temperature"]
        capacity_factors[bus] = np.interp((temp - temp.mean()).values, x, y)

    return capacity_factors


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_egs_potentials",
            simpl="",
            clusters=37,
        )

    egs_config = snakemake.params["sector"]["enhanced_geothermal"]
    costs_config = snakemake.params["costs"]

    sustainability_factor = egs_config["sustainability_factor"]
    # the share of heat that is replenished from the earth's core.
    # we are not constraining ourselves to the sustainable share, but
    # inversely apply it to our underlying data, which refers to the
    # sustainable heat. Source: Relative magnitude of sustainable heat vs
    # nonsustainable heat in the paper "From hot rock to useful energy..."

    egs_data = prepare_egs_data(snakemake.input.egs_cost)

    if egs_config["optimism"]:
        egs_data = egs_data[(year := costs_config["year"])]
        logger.info(
            f"EGS optimism! Building EGS potentials with costs estimated for {year}."
        )

    else:
        egs_data = egs_data[(default_year := 2020)]
        logger.info(
            f"No EGS optimism! Building EGS potentials with {default_year} costs."
        )

    egs_data = egs_data.loc[egs_data["PowerSust"] > 0].reset_index(drop=True)
    egs_regions = egs_data.geometry

    network_regions = (
        gpd.read_file(snakemake.input.regions)
        .set_index("name", drop=True)
        .set_crs(epsg=4326)
    )

    overlap_matrix = pd.DataFrame(
        index=network_regions.index,
        columns=egs_data.index,
    )

    for name, polygon in network_regions.geometry.items():
        overlap_matrix.loc[name] = (
            egs_regions.intersection(polygon).area
        ) / egs_regions.area

    overlap_matrix.to_csv(snakemake.output["egs_overlap"])

    # consider not only replenished heat
    egs_data["p_nom_max"] = egs_data["PowerSust"] / sustainability_factor

    egs_data[["p_nom_max", "CAPEX"]].to_csv(snakemake.output["egs_potentials"])

    capacity_factors = get_capacity_factors(
        snakemake.input["regions"],
        snakemake.input["air_temperature"],
    )

    capacity_factors.to_csv(snakemake.output["egs_capacity_factors"])
