# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Rasters the vector data of the `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas onto all cutout regions.

Relevant Settings
-----------------

.. code:: yaml

    renewable:
        {technology}:
            cutout:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`renewable_cf`

Inputs
------

- ``data/bundle/natura/Natura2000_end2015.shp``: `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas.

    .. image:: ../img/natura.png
        :scale: 33 %

Outputs
-------

- ``resources/natura.tiff``: Rasterized version of `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas to reduce computation times.

    .. image:: ../img/natura.png
        :scale: 33 %

Description
-----------

"""

import logging
from _helpers import configure_logging

import atlite
import geopandas as gpd
import rasterio as rio
from rasterio.features import geometry_mask
from rasterio.warp import transform_bounds

logger = logging.getLogger(__name__)


def determine_cutout_xXyY(cutout_name):
    cutout = atlite.Cutout(cutout_name)
    assert cutout.crs.to_epsg() == 4326
    x, X, y, Y = cutout.extent
    dx, dy = cutout.dx, cutout.dy
    return [x - dx/2., X + dx/2., y - dy/2., Y + dy/2.]


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_continental_raster')
    configure_logging(snakemake)


    cutouts = snakemake.input.cutouts
    xs, Xs, ys, Ys = zip(*(determine_cutout_xXyY(cutout) for cutout in cutouts))
    bounds = transform_bounds(4326, 3035, min(xs), min(ys), max(Xs), max(Ys))
    res = 100 # 100 meter resolution

    # adjusted boundaries
    left, bottom = [(b // res)* res for b in bounds[:2]]
    right, top = [(b // res + 1) * res for b in bounds[2:]]
    out_shape = int((right - left) / res), int((top - bottom) // res)

    transform = rio.Affine(res, 0, left, 0, -res, top)
    shapes = gpd.read_file(snakemake.input.shapes).to_crs(3035)
    raster = ~geometry_mask(shapes.geometry, out_shape[::-1], transform)


    with rio.open(snakemake.output[0], 'w', driver='GTiff', dtype=rio.uint8,
                  count=1, transform=transform, crs=3035,
                  width=raster.shape[1], height=raster.shape[0]) as dst:
        dst.write(raster, indexes=1)

