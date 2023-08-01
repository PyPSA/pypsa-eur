# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import functools
import logging
import time

import atlite
import geopandas as gpd
import numpy as np
from _helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "determine_availability_matrix_MD_UA", technology="solar"
        )
    configure_logging(snakemake)

    nprocesses = snakemake.config["atlite"].get("nprocesses")
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
            codes = [
                20,
                30,
                40,
                60,
                100,
                111,
                112,
                113,
                114,
                115,
                116,
                121,
                122,
                123,
                124,
                125,
                126,
            ]
        elif snakemake.wildcards.technology == "offshore-ac":
            codes = [80, 200]
        elif snakemake.wildcards.technology == "offshore-dc":
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

    kwargs = dict(nprocesses=nprocesses, disable_progressbar=noprogress)
    if noprogress:
        logger.info("Calculate landuse availabilities...")
        start = time.time()
        availability = cutout.availabilitymatrix(regions, excluder, **kwargs)
        duration = time.time() - start
        logger.info(f"Completed availability calculation ({duration:2.2f}s)")
    else:
        availability = cutout.availabilitymatrix(regions, excluder, **kwargs)

    # Limit results only to buses for UA and MD
    buses = regions.loc[regions["country"].isin(["UA", "MD"])].index.values
    availability = availability.sel(bus=buses)

    # Save and plot for verification
    availability.to_netcdf(snakemake.output.availability_matrix)
