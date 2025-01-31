# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Transforms the global ship density data from the `World Bank Data Catalogue.

<https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density>`_
to the size of the considered cutout. The global ship density raster is later
used for the exclusion when calculating the offshore potentials.

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

import atlite
import rioxarray
from _helpers import configure_logging, set_scenario_config

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
