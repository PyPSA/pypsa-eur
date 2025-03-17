# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
This rule extracts potential and cost for electricity generation through
enhanced geothermal systems.

For this, we use data from "From hot rock to useful energy..." by Aghahosseini, Breyer (2020)
'https://doi.org/10.1016/j.apenergy.2020.115769'
Note that we input data used here is not the same as in the paper, but was passed on by the authors.

The data provides a lon-lat gridded map of Europe (1° x 1°), with each grid cell assigned
a heat potential (in GWh) and a cost (in EUR/MW).

This scripts overlays that map with the network's regions, and builds a csv with CAPEX, OPEX and p_nom_max
"""

import json
import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from _helpers import configure_logging, set_scenario_config
from shapely.geometry import Polygon

logger = logging.getLogger(__name__)


def prepare_egs_data(egs_file):
    """
    Processes the original .json file EGS data to a more human-readable format.
    """
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


def prepare_capex(prepared_data):
    """
    The source paper provides only data for year and regions where LCOE <
    100Euro/MWh. However, this implementations starts with the costs for 2020
    for all regions and then adjusts the costs according to the user's chosen
    setting in the config file.

    As such, for regions where cost data is available only from, say,
    2035, we need to reverse-engineer the costs for 2020. This is done
    in the following (unfortunately verbose) function.
    """

    default_year = 2020

    # obtains all available CAPEX data
    capex_df = pd.DataFrame(columns=prepared_data.keys())

    for year in capex_df.columns:
        year_data = prepared_data[year].groupby("geometry").mean().reset_index()

        for g in year_data.geometry:
            if g not in year_data.geometry.tolist():
                # weird but apparently necessary
                continue

            capex_df.loc[g, year] = year_data.loc[
                year_data.geometry == g, "CAPEX"
            ].values[0]

    capex_df = capex_df.loc[:, default_year:]

    # fill up missing values assuming cost reduction factors similar to existing values
    for sooner, later in zip(capex_df.columns[::-1][1:], capex_df.columns[::-1]):
        missing_mask = capex_df[sooner].isna()
        cr_factor = (
            capex_df.loc[~missing_mask, later] / capex_df.loc[~missing_mask, sooner]
        )

        capex_df.loc[missing_mask, sooner] = (
            capex_df.loc[missing_mask, later] / cr_factor.mean()
        )

    # harmonice capacity and CAPEX
    p_nom_max = prepared_data[2050].groupby("geometry")["PowerSust"].mean()
    p_nom_max = p_nom_max.loc[p_nom_max > 0]

    capex_df = capex_df.loc[p_nom_max.index]

    data = (
        pd.concat((capex_df[default_year], p_nom_max), axis=1)
        .reset_index()
        .rename(columns={2020: "CAPEX"})
    )
    return gpd.GeoDataFrame(data, geometry=data.geometry)


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
            clusters=37,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    egs_config = snakemake.params["sector"]["enhanced_geothermal"]
    costs_config = snakemake.params["costs"]

    egs_data = prepare_egs_data(snakemake.input.egs_cost)
    egs_data = prepare_capex(egs_data)

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

    # the share of heat that is replenished from the earth's core.
    # we are not constraining ourselves to the sustainable share, but
    # inversely apply it to our underlying data, which refers to the
    # sustainable heat. Source: Relative magnitude of sustainable heat vs
    # nonsustainable heat in the paper "From hot rock to useful energy..."
    sustainability_factor = egs_config["sustainability_factor"]
    egs_data["p_nom_max"] = egs_data["PowerSust"] / sustainability_factor

    egs_data[["p_nom_max", "CAPEX"]].to_csv(snakemake.output["egs_potentials"])

    capacity_factors = get_capacity_factors(
        snakemake.input["regions"],
        snakemake.input["air_temperature"],
    )

    capacity_factors.to_csv(snakemake.output["egs_capacity_factors"])
