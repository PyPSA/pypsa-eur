# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Create land elibility analysis for Ukraine and Moldova with different datasets.
"""

import functools
import logging
import time

import atlite
import fiona
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from _helpers import configure_logging, set_scenario_config
from atlite.gis import shape_availability
from rasterio.plot import show

logger = logging.getLogger(__name__)


def get_wdpa_layer_name(wdpa_fn, layer_substring):
    """
    Get layername from file "wdpa_fn" whose name contains "layer_substring".
    """
    l = fiona.listlayers(wdpa_fn)
    return [_ for _ in l if layer_substring in _][0]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "determine_availability_matrix_MD_UA", technology="solar"
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nprocesses = None  # snakemake.config["atlite"].get("nprocesses")
    noprogress = not snakemake.config["atlite"].get("show_progress", True)
    config = snakemake.config["renewable"][snakemake.wildcards.technology]

    cutout = atlite.Cutout(snakemake.input.cutout)
    regions = (
        gpd.read_file(snakemake.input.regions).set_index("name").rename_axis("bus")
    )
    buses = regions.index

    excluder = atlite.ExclusionContainer(crs=3035, res=100)

    corine = config.get("corine", {})
    if "grid_codes" in corine:
        # Land cover codes to emulate CORINE results
        if snakemake.wildcards.technology == "solar":
            codes = [20, 30, 40, 50, 60, 90, 100]
        elif snakemake.wildcards.technology == "onwind":
            codes = [20, 30, 40, 60, 100]
        elif snakemake.wildcards.technology == "offwind-ac":
            codes = [80, 200]
        elif snakemake.wildcards.technology == "offwind-dc":
            codes = [80, 200]
        else:
            assert False, "technology not supported"

        excluder.add_raster(
            snakemake.input.copernicus, codes=codes, invert=True, crs="EPSG:4326"
        )
    if "distance" in corine and corine.get("distance", 0.0) > 0.0:
        # Land cover codes to emulate CORINE results
        if snakemake.wildcards.technology == "onwind":
            codes = [50]
        else:
            assert False, "technology not supported"

        buffer = corine["distance"]
        excluder.add_raster(
            snakemake.input.copernicus, codes=codes, buffer=buffer, crs="EPSG:4326"
        )

    if config["natura"]:
        wdpa_fn = (
            snakemake.input.wdpa_marine
            if "offwind" in snakemake.wildcards.technology
            else snakemake.input.wdpa
        )
        layer = get_wdpa_layer_name(wdpa_fn, "polygons")
        wdpa = gpd.read_file(
            wdpa_fn,
            bbox=regions.geometry,
            layer=layer,
        ).to_crs(3035)
        if not wdpa.empty:
            excluder.add_geometry(wdpa.geometry)

        layer = get_wdpa_layer_name(wdpa_fn, "points")
        wdpa_pts = gpd.read_file(
            wdpa_fn,
            bbox=regions.geometry,
            layer=layer,
        ).to_crs(3035)
        wdpa_pts = wdpa_pts[wdpa_pts["REP_AREA"] > 1]
        wdpa_pts["buffer_radius"] = np.sqrt(wdpa_pts["REP_AREA"] / np.pi) * 1000
        wdpa_pts = wdpa_pts.set_geometry(
            wdpa_pts["geometry"].buffer(wdpa_pts["buffer_radius"])
        )
        if not wdpa_pts.empty:
            excluder.add_geometry(wdpa_pts.geometry)

    if "max_depth" in config:
        # lambda not supported for atlite + multiprocessing
        # use named function np.greater with partially frozen argument instead
        # and exclude areas where: -max_depth > grid cell depth
        func = functools.partial(np.greater, -config["max_depth"])
        excluder.add_raster(snakemake.input.gebco, codes=func, crs=4236, nodata=-1000)

    if "min_shore_distance" in config:
        buffer = config["min_shore_distance"]
        excluder.add_geometry(snakemake.input.country_shapes, buffer=buffer)

    if "max_shore_distance" in config:
        buffer = config["max_shore_distance"]
        excluder.add_geometry(
            snakemake.input.country_shapes, buffer=buffer, invert=True
        )

    if "ship_threshold" in config:
        shipping_threshold = config["ship_threshold"] * 8760 * 6
        func = functools.partial(np.less, shipping_threshold)
        excluder.add_raster(
            snakemake.input.ship_density, codes=func, crs=4326, allow_no_overlap=True
        )

    kwargs = dict(nprocesses=nprocesses, disable_progressbar=noprogress)
    if noprogress:
        logger.info("Calculate landuse availabilities...")
        start = time.time()
        availability = cutout.availabilitymatrix(regions, excluder, **kwargs)
        duration = time.time() - start
        logger.info(f"Completed availability calculation ({duration:2.2f}s)")
    else:
        availability = cutout.availabilitymatrix(regions, excluder, **kwargs)

    regions_geometry = regions.to_crs(3035).geometry
    band, transform = shape_availability(regions_geometry, excluder)
    fig, ax = plt.subplots(figsize=(4, 8))
    gpd.GeoSeries(regions_geometry.unary_union).plot(ax=ax, color="none")
    show(band, transform=transform, cmap="Greens", ax=ax)
    plt.axis("off")
    plt.savefig(snakemake.output.availability_map, bbox_inches="tight", dpi=500)

    # Limit results only to buses for UA and MD
    buses = regions.loc[regions["country"].isin(["UA", "MD"])].index.values
    availability = availability.sel(bus=buses)

    # Save and plot for verification
    availability.to_netcdf(snakemake.output.availability_matrix)
