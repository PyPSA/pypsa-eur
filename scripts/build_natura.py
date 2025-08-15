# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Rasters the vector data of the `Natura 2000.

<https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas onto all
cutout regions.

Relevant Settings
-----------------

.. code:: yaml

    renewable:
        {technology}:
            cutout:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`renewable_cf`

Inputs
------

- ``data/bundle/natura/Natura2000_end2015.shp``: `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas.

    .. image:: img/natura.png
        :scale: 33 %

Outputs
-------

- ``resources/natura.tiff``: Rasterized version of `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas to reduce computation times.

    .. image:: img/natura.png
        :scale: 33 %

Description
-----------
"""

import logging
import shutil
from pathlib import Path

import atlite
import geopandas as gpd
import rasterio as rio
from _helpers import configure_logging, set_scenario_config
from rasterio.features import geometry_mask
from rasterio.warp import transform_bounds

logger = logging.getLogger(__name__)


def determine_cutout_xXyY(cutout_name):
    """
    Determine the full extent of a cutout.

    Since the coordinates of the cutout data are given as the
    center of the grid cells, the extent of the cutout is
    calculated by adding/subtracting half of the grid cell size.


    Parameters
    ----------
    cutout_name : str
        Path to the cutout.

    Returns
    -------
    A list of extent coordinates in the order [x, X, y, Y].
    """
    cutout = atlite.Cutout(cutout_name)
    assert cutout.crs.to_epsg() == 4326
    x, X, y, Y = cutout.extent
    dx, dy = cutout.dx, cutout.dy
    return [x - dx / 2.0, X + dx / 2.0, y - dy / 2.0, Y + dy / 2.0]


def get_transform_and_shape(bounds, res):
    left, bottom = [(b // res) * res for b in bounds[:2]]
    right, top = [(b // res + 1) * res for b in bounds[2:]]
    shape = int((top - bottom) // res), int((right - left) / res)
    transform = rio.Affine(res, 0, left, 0, -res, top)
    return transform, shape


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_natura_raster")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Move the downloaded archive and unpack it
    shutil.move(snakemake.input["online"], snakemake.output["zip"])
    shutil.unpack_archive(snakemake.output["zip"], snakemake.output["raw"])

    # Find the shapefile in the unpacked nested directory
    shapefile = list(Path(snakemake.output["raw"]).rglob("*.shp"))[0]

    x, X, y, Y = determine_cutout_xXyY(snakemake.input["cutout"])
    bounds = transform_bounds(4326, 3035, x, y, X, Y)
    transform, out_shape = get_transform_and_shape(bounds, res=100)

    # adjusted boundaries
    shapes = gpd.read_file(shapefile).to_crs(3035)
    raster = ~geometry_mask(shapes.geometry, out_shape, transform)
    raster = raster.astype(rio.uint8)

    with rio.open(
        snakemake.output["raster"],
        "w",
        driver="GTiff",
        dtype=rio.uint8,
        count=1,
        transform=transform,
        crs=3035,
        compress="lzw",
        width=raster.shape[1],
        height=raster.shape[0],
    ) as dst:
        dst.write(raster, indexes=1)
