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

- ``data/shipdensity_global.tif``: `Global ship density from <https://datacatalog.worldbank.org/search/dataset/0037580/`.

Description
-----------
"""

import logging
import zipfile
from pathlib import Path

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    output = Path(snakemake.output[0])
    import zipfile

    with zipfile.ZipFile(snakemake.input[0]) as zip_f:
        zip_f.extract(output.name, path=output.parent)
