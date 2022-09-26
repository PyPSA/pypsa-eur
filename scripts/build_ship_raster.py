# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Transforms the global ship density data from
https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-
Traffic-Density to the size of the considered cutout. The global ship density
raster is later used for the exclusion when calculating the offshore
potentials.

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

- ``data/bundle/shipdensity/shipdensity_global.zip``: `Global ship density from <https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density>`.

Outputs
-------

- ``resources/europe_shipdensity_raster.nc``: Reduced version of `Global ship density from <https://datacatalog.worldbank.org/search/dataset/0037580/` to reduce computation time.

Description
-----------
"""

import logging
import os
import zipfile

import xarray as xr
from _helpers import configure_logging
from build_natura_raster import determine_cutout_xXyY

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ship_raster")
    configure_logging(snakemake)

    cutouts = snakemake.input.cutouts
    xs, Xs, ys, Ys = zip(*(determine_cutout_xXyY(cutout) for cutout in cutouts))

    with zipfile.ZipFile(snakemake.input.ship_density) as zip_f:
        zip_f.extract("shipdensity_global.tif")
        with xr.open_rasterio("shipdensity_global.tif") as ship_density:
            ship_density = ship_density.drop(["band"]).sel(
                x=slice(min(xs), max(Xs)), y=slice(max(Ys), min(ys))
            )
            ship_density.to_netcdf(snakemake.output[0])

    os.remove("shipdensity_global.tif")
