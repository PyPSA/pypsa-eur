# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Creates GIS shape files of country boundaries and the Europe bounding shape,
derived from NUTS3 regions and offshore EEZ shapes.
"""

import logging
from operator import attrgetter

import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def country_cover(country_shapes, eez_shapes=None):
    shapes = country_shapes
    if eez_shapes is not None:
        shapes = pd.concat([shapes, eez_shapes])
    europe_shape = shapes.union_all()
    if isinstance(europe_shape, MultiPolygon):
        europe_shape = max(europe_shape.geoms, key=attrgetter("area"))
    return Polygon(shell=europe_shape.exterior)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_shapes")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nuts3_shapes = gpd.read_file(snakemake.input.nuts3_shapes).set_index("index")
    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes).set_index("name")

    country_shapes = nuts3_shapes.groupby("country")["geometry"].apply(
        lambda x: x.union_all()
    )
    country_shapes.crs = nuts3_shapes.crs
    country_shapes.index.name = "name"
    country_shapes.reset_index().to_file(snakemake.output.country_shapes)

    europe_shape = gpd.GeoDataFrame(
        geometry=[country_cover(country_shapes, offshore_shapes.geometry)],
        crs=country_shapes.crs,
    )
    europe_shape.reset_index().to_file(snakemake.output.europe_shape)
