# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Transforms the global ship density data from the `World Bank Data Catalogue.

<https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density>`_
to the size of the considered cutout. The global ship density raster is later
used for the exclusion when calculating the offshore potentials.

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

- ``data/bundle/shipdensity/shipdensity_global.zip``: Global shipping traffic
  density from `World Bank Data Catalogue
  <https://datacatalog.worldbank.org/search/dataset/0037580/>`_.

Outputs
-------

- ``resources/europe_shipdensity_raster.nc``: Reduced version of global shipping
  traffic density from `World Bank Data Catalogue
  <https://datacatalog.worldbank.org/search/dataset/0037580/>`_ to reduce
  computation time.

Description
-----------
"""

import logging
import zipfile
from pathlib import Path

import rioxarray
from _helpers import configure_logging, set_scenario_config
from build_natura_raster import determine_cutout_xXyY

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ship_raster")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    x, X, y, Y = determine_cutout_xXyY(snakemake.input.cutout)

    with zipfile.ZipFile(snakemake.input.ship_density) as zip_f:
        resources = Path(snakemake.output[0]).parent
        fn = "shipdensity_global.tif"
        zip_f.extract(fn, resources)
    with rioxarray.open_rasterio(resources / fn) as ship_density:
        ship_density = ship_density.drop_vars(["band"]).sel(
            x=slice(x, X), y=slice(Y, y)
        )
        ship_density.rio.to_raster(snakemake.output[0])

    (resources / fn).unlink()
