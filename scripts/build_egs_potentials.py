# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 @LukasFranken, The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This rule extracts potential and cost for electricity generation through enhanced geothermal systems

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
import pandas as pd
import geopandas as gpd

from shapely.geometry import Polygon


def prepare_egs_data(egs_file):
    
    with open(egs_file) as f:
        jsondata = json.load(f)

    def point_to_square(p, lon_extent=1., lat_extent=1.):

        try:
            x, y = p.coords.xy[0][0], p.coords.xy[1][0]
        except IndexError:
            return p
        
        return Polygon([
            [x-lon_extent/2, y-lat_extent/2],
            [x-lon_extent/2, y+lat_extent/2],
            [x+lon_extent/2, y+lat_extent/2],
            [x+lon_extent/2, y-lat_extent/2],
            ])

    years = [2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050]
    lcoes = ["LCOE50", "LCOE100", "LCOE150"]

    egs_data = dict()

    for year in years:
        df = pd.DataFrame(columns=["Lon", "Lat", "CAPEX", "HeatSust", "PowerSust"])

        for lcoe in lcoes:

            for country_data in jsondata[lcoe]:
                try:
                    country_df = pd.DataFrame(columns=df.columns, 
                                              index=range(len(country_data[0][years.index(year)]["Lon"])))
                except TypeError:
                    country_df = pd.DataFrame(columns=df.columns, index=range(0))

                for col in df.columns:
                    country_df[col] = country_data[0][years.index(year)][col]

                df = pd.concat((df, country_df.dropna()), axis=0, ignore_index=True)

        gdf = gpd.GeoDataFrame(
            df.drop(
                columns=["Lon", "Lat"]
                ),
            geometry=gpd.points_from_xy(df.Lon, df.Lat)
            ).reset_index(drop=True)
    
        gdf["geometry"] = gdf.geometry.apply(lambda geom: point_to_square(geom))
        egs_data[year] = gdf
    
    return egs_data



if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_egs_potentials",
            simpl="",
            clusters=37,
        )
    
    sustainability_factor = 0.0025 # factor sustainable p_nom vs p_nom
    config = snakemake.config

    egs_data = prepare_egs_data(snakemake.input.egs_cost)

    if config["sector"]["enhanced_geothermal_optimism"]:
    
        egs_data = egs_data[(year := config["costs"]["year"])]
        logger.info(f"EGS optimism! Builing EGS potentials with costs estimated for {year}.")

    else:
        egs_data = egs_data[(default_year := 2020)]
        logger.info(f"No EGS optimism! Building EGS potentials with {default_year} costs.")

    egs_data.index = egs_data.geometry.astype(str)
    egs_shapes = egs_data.geometry
    
    network_shapes = (
        gpd.read_file(snakemake.input.shapes)
        .set_index("name", drop=True)
        .set_crs(epsg=4326)
    )

    overlap_matrix = (
        pd.DataFrame(
            index=network_shapes.index,
            columns=(egs_shapes := egs_data.geometry).astype(str).values)
    )

    for name, polygon in network_shapes.geometry.items():
        overlap_matrix.loc[name] = (
            egs_shapes
            .intersection(polygon).area
        ) / egs_shapes.area

    overlap_matrix.to_csv(snakemake.output["egs_overlap"])

    egs_data["p_nom_max"] = egs_data["PowerSust"] / sustainability_factor    
    egs_data[["p_nom_max", "CAPEX"]].to_csv(snakemake.output["egs_potentials"])
