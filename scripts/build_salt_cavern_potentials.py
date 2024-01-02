# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build salt cavern potentials for hydrogen storage.

Technical Potential of Salt Caverns for Hydrogen Storage in Europe CC-BY
4.0
https://doi.org/10.20944/preprints201910.0187.v1
https://doi.org/10.1016/j.ijhydene.2019.12.161

Figure 6. Distribution of potential salt cavern sites across Europe with their corresponding
energy densities (cavern storage potential divided by the volume).

Figure 7. Total cavern storage potential in European countries
classified as onshore, offshore and within 50 km of shore.

The regional distribution is taken from the map (Figure 6) and scaled to the
capacities from the bar chart split by nearshore (<50km from sea),
onshore (>50km from sea), offshore (Figure 7).
"""


import geopandas as gpd
import pandas as pd


def concat_gdf(gdf_list, crs="EPSG:4326"):
    """
    Concatenate multiple geopandas dataframes with common coordinate reference
    system (crs).
    """
    return gpd.GeoDataFrame(pd.concat(gdf_list), crs=crs)


def load_bus_regions(onshore_path, offshore_path):
    """
    Load pypsa-eur on- and offshore regions and concat.
    """
    bus_regions_offshore = gpd.read_file(offshore_path)
    bus_regions_onshore = gpd.read_file(onshore_path)
    bus_regions = concat_gdf([bus_regions_offshore, bus_regions_onshore])
    bus_regions = bus_regions.dissolve(by="name", aggfunc="sum")

    return bus_regions


def area(gdf):
    """
    Returns area of GeoDataFrame geometries in square kilometers.
    """
    return gdf.to_crs(epsg=3035).area.div(1e6)


def salt_cavern_potential_by_region(caverns, regions):
    # calculate area of caverns shapes
    caverns["area_caverns"] = area(caverns)

    overlay = gpd.overlay(regions.reset_index(), caverns, keep_geom_type=True)

    # calculate share of cavern area inside region
    overlay["share"] = area(overlay) / overlay["area_caverns"]

    overlay["e_nom"] = overlay.eval(
        "capacity_per_area * share * area_caverns / 1000"
    )  # TWh

    return overlay.groupby(["name", "storage_type"]).e_nom.sum().unstack("storage_type")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_salt_cavern_potentials", simpl="", clusters="37"
        )

    fn_onshore = snakemake.input.regions_onshore
    fn_offshore = snakemake.input.regions_offshore

    regions = load_bus_regions(fn_onshore, fn_offshore)

    caverns = gpd.read_file(snakemake.input.salt_caverns)  # GWh/sqkm

    caverns_regions = salt_cavern_potential_by_region(caverns, regions)

    caverns_regions.to_csv(snakemake.output.h2_cavern_potential)
