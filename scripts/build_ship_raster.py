# SPDX-FileCopyrightText: : 2017-2022 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Transforms the global ship density data from https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density to the size of the considered cutout. The global ship density raster is later used for the exclusion when calculating the offshore potentials.

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

- ``resources/natura.tiff``: Reduced version of `Global ship density from <https://datacatalog.worldbank.org/search/dataset/0037580/` to reduce computation times.

Description
-----------

"""

import logging
from _helpers import configure_logging
from build_natura_raster import determine_cutout_xXyY

import zipfile
import xarray
import os

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_ship_raster')
    configure_logging(snakemake)

    cutouts = snakemake.input.cutouts
    xs, Xs, ys, Ys = zip(*(determine_cutout_xXyY(cutout) for cutout in cutouts))
    
    with zipfile.ZipFile(snakemake.input.ship_density) as zip_f:
        zip_f.extract("shipdensity_global.tif")
        ship_density=xarray.open_dataarray("shipdensity_global.tif", engine="rasterio")
        os.remove("shipdensity_global.tif") 

    ship_density=ship_density.drop("band").sel(x=slice(min(xs),max(Xs)), y=slice(max(Ys),min(ys)))
    
    ship_density.to_netcdf(snakemake.output[0])

