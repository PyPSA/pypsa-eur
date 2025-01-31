# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create land elibility analysis for Ukraine and Moldova with different datasets.
"""

import functools
import logging
import os
import time
from tempfile import NamedTemporaryFile

import atlite
import fiona
import geopandas as gpd
import numpy as np
from _helpers import configure_logging, set_scenario_config

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
            "determine_availability_matrix_MD_UA", clusters=100, technology="solar"
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    noprogress = not snakemake.config["atlite"].get("show_progress", True)
    config = snakemake.params["renewable"][snakemake.wildcards.technology]

    cutout = atlite.Cutout(snakemake.input.cutout)
    regions = (
        gpd.read_file(snakemake.input.regions).set_index("name").rename_axis("bus")
    )
    # Limit to "UA" and "MD" regions
    buses = regions.filter(regex="(UA|MD)", axis=0).index.values
    regions = regions.loc[buses]

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

        # temporary file needed for parallelization
        with NamedTemporaryFile(suffix=".geojson", delete=False) as f:
            plg_tmp_fn = f.name
        if not wdpa.empty:
            wdpa[["geometry"]].to_file(plg_tmp_fn)
            while not os.path.exists(plg_tmp_fn):
                time.sleep(1)
            excluder.add_geometry(plg_tmp_fn)

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

        # temporary file needed for parallelization
        with NamedTemporaryFile(suffix=".geojson", delete=False) as f:
            pts_tmp_fn = f.name
        if not wdpa_pts.empty:
            wdpa_pts[["geometry"]].to_file(pts_tmp_fn)
            while not os.path.exists(pts_tmp_fn):
                time.sleep(1)
            excluder.add_geometry(pts_tmp_fn)

    if config.get("max_depth"):
        # lambda not supported for atlite + multiprocessing
        # use named function np.greater with partially frozen argument instead
        # and exclude areas where: -max_depth > grid cell depth
        func = functools.partial(np.greater, -config["max_depth"])
        excluder.add_raster(snakemake.input.gebco, codes=func, crs=4236, nodata=-1000)

    if config.get("min_shore_distance"):
        buffer = config["min_shore_distance"]
        excluder.add_geometry(snakemake.input.country_shapes, buffer=buffer)

    if config.get("max_shore_distance"):
        buffer = config["max_shore_distance"]
        excluder.add_geometry(
            snakemake.input.country_shapes, buffer=buffer, invert=True
        )

    if config.get("ship_threshold"):
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

    for fn in [pts_tmp_fn, plg_tmp_fn]:
        if os.path.exists(fn):
            os.remove(fn)

    availability = availability.sel(bus=buses)

    # Save and plot for verification
    availability.to_netcdf(snakemake.output.availability_matrix)
