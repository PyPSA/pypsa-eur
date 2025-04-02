# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
The script performs a land eligibility analysis of what share of land is
availability for developing the selected technology at each cutout grid cell.
The script uses the `atlite <https://github.com/pypsa/atlite>`_ library and
several GIS datasets like the CORINE land use data, LUISA land use data,
Natura2000 nature reserves, GEBCO bathymetry data, and shipping lanes.


Inputs
------

- ``data/bundle/corine/g250_clc06_V18_5.tif``: `CORINE Land Cover (CLC)
  <https://land.copernicus.eu/pan-european/corine-land-cover>`_ inventory on `44
  classes <https://wiki.openstreetmap.org/wiki/Corine_Land_Cover#Tagging>`_ of
  land use (e.g. forests, arable land, industrial, urban areas) at 100m
  resolution.

    .. image:: img/corine.png
        :scale: 33 %

- ``data/LUISA_basemap_020321_50m.tif``: `LUISA Base Map
  <https://publications.jrc.ec.europa.eu/repository/handle/JRC124621>`_ land
  coverage dataset at 50m resolution similar to CORINE. For codes in relation to
  CORINE land cover, see `Annex 1 of the technical documentation
  <https://publications.jrc.ec.europa.eu/repository/bitstream/JRC124621/technical_report_luisa_basemap_2018_v7_final.pdf>`_.

- ``data/bundle/gebco/GEBCO_2014_2D.nc``: A `bathymetric
  <https://en.wikipedia.org/wiki/Bathymetry>`_ data set with a global terrain
  model for ocean and land at 15 arc-second intervals by the `General
  Bathymetric Chart of the Oceans (GEBCO)
  <https://www.gebco.net/data_and_products/gridded_bathymetry_data/>`_.

    .. image:: img/gebco_2019_grid_image.jpg
        :scale: 50 %

    **Source:** `GEBCO
    <https://www.gebco.net/data_and_products/images/gebco_2019_grid_image.jpg>`_

- ``resources/natura.tiff``: confer :ref:`natura`
- ``resources/offshore_shapes.geojson``: confer :ref:`shapes`
- ``resources/regions_onshore_base_s_{clusters}.geojson``: (if not offshore
  wind), confer :ref:`busregions`
- ``resources/regions_offshore_base_s_{clusters}.geojson``: (if offshore wind),
  :ref:`busregions`
- ``"cutouts/" + params["renewable"][{technology}]['cutout']``: :ref:`cutout`
- ``networks/_base_s_{clusters}.nc``: :ref:`base`

Outputs
-------

- ``resources/availability_matrix_{clusters_{technology}.nc``
"""

import functools
import logging
import time

import atlite
import geopandas as gpd
import numpy as np
import xarray as xr
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_renewable_profiles", clusters=100, technology="onwind"
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    noprogress = snakemake.config["run"].get("disable_progressbar", True)
    noprogress = noprogress or not snakemake.config["atlite"]["show_progress"]
    technology = snakemake.wildcards.technology
    params = snakemake.params.renewable[technology]

    cutout = atlite.Cutout(snakemake.input.cutout)
    regions = gpd.read_file(snakemake.input.regions)
    assert not regions.empty, (
        f"List of regions in {snakemake.input.regions} is empty, please "
        "disable the corresponding renewable technology"
    )
    # do not pull up, set_index does not work if geo dataframe is empty
    regions = regions.set_index("name").rename_axis("bus")

    res = params.get("excluder_resolution", 100)
    excluder = atlite.ExclusionContainer(crs=3035, res=res)

    if params["natura"]:
        excluder.add_raster(snakemake.input.natura, nodata=0, allow_no_overlap=True)

    for dataset in ["corine", "luisa"]:
        kwargs = {"nodata": 0} if dataset == "luisa" else {}
        settings = params.get(dataset, {})
        if not settings:
            continue
        if dataset == "luisa" and res > 50:
            logger.info(
                "LUISA data is available at 50m resolution, "
                f"but coarser {res}m resolution is used."
            )
        if isinstance(settings, list):
            settings = {"grid_codes": settings}
        if "grid_codes" in settings:
            codes = settings["grid_codes"]
            excluder.add_raster(
                snakemake.input[dataset], codes=codes, invert=True, crs=3035, **kwargs
            )
        if settings.get("distance", 0.0) > 0.0:
            codes = settings["distance_grid_codes"]
            buffer = settings["distance"]
            excluder.add_raster(
                snakemake.input[dataset], codes=codes, buffer=buffer, crs=3035, **kwargs
            )

    if params.get("ship_threshold"):
        shipping_threshold = (
            params["ship_threshold"] * 8760 * 6
        )  # approximation because 6 years of data which is hourly collected
        func = functools.partial(np.less, shipping_threshold)
        excluder.add_raster(
            snakemake.input.ship_density, codes=func, crs=4326, allow_no_overlap=True
        )

    if params.get("max_depth"):
        # lambda not supported for atlite + multiprocessing
        # use named function np.greater with partially frozen argument instead
        # and exclude areas where: -max_depth > grid cell depth
        func = functools.partial(np.greater, -params["max_depth"])
        excluder.add_raster(snakemake.input.gebco, codes=func, crs=4326, nodata=-1000)

    if params.get("min_depth"):
        func = functools.partial(np.greater, -params["min_depth"])
        excluder.add_raster(
            snakemake.input.gebco, codes=func, crs=4326, nodata=-1000, invert=True
        )

    if "min_shore_distance" in params:
        buffer = params["min_shore_distance"]
        excluder.add_geometry(snakemake.input.country_shapes, buffer=buffer)

    if "max_shore_distance" in params:
        buffer = params["max_shore_distance"]
        excluder.add_geometry(
            snakemake.input.country_shapes, buffer=buffer, invert=True
        )

    logger.info(f"Calculate landuse availability for {technology}...")
    start = time.time()

    kwargs = dict(nprocesses=nprocesses, disable_progressbar=noprogress)
    availability = cutout.availabilitymatrix(regions, excluder, **kwargs)

    duration = time.time() - start
    logger.info(
        f"Completed landuse availability calculation for {technology} ({duration:2.2f}s)"
    )

    # For Moldova and Ukraine: Overwrite parts not covered by Corine with
    # externally determined available areas
    if "availability_matrix_MD_UA" in snakemake.input.keys():
        availability_MDUA = xr.open_dataarray(
            snakemake.input["availability_matrix_MD_UA"]
        )
        availability.loc[availability_MDUA.coords] = availability_MDUA

    availability.to_netcdf(snakemake.output[0])
