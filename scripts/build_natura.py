# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Rasterises the [Natura 2000](https://en.wikipedia.org/wiki/Natura_2000) natural protection areas onto the extent of the weather cutout.

The Natura 2000 shapefile is unpacked, projected to the European LAEA grid
(EPSG:3035) and burned into a binary raster with 100 m resolution covering the
cutout extent. The raster serves as an exclusion layer in the land availability
calculations, which is much faster than intersecting the vector data directly.
"""

import logging
import shutil
from pathlib import Path

import geopandas as gpd
import rasterio as rio
from _helpers import configure_logging, determine_cutout_xXyY, set_scenario_config
from rasterio.features import geometry_mask
from rasterio.warp import transform_bounds

logger = logging.getLogger(__name__)


def get_transform_and_shape(
    bounds: tuple[float, float, float, float], res: float
) -> tuple[rio.Affine, tuple[int, int]]:
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
