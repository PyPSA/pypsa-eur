# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Creates GIS shape files of offshore exclusive economic zones (EEZ) from the
maritime shapes of the geo_boundaries module.
"""

import logging

import geopandas as gpd

from scripts._helpers import (
    _simplify_polys,
    configure_logging,
    read_geo_boundaries,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def build_offshore_shapes(maritime: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Dissolve the maritime shapes of the geo_boundaries module per country.

    Parameters
    ----------
    maritime : geopandas.GeoDataFrame
        Maritime rows of the module output with an ISO2 ``country`` column.

    Returns
    -------
    geopandas.GeoDataFrame
        Offshore shapes per country, indexed by ISO2 ``name``.
    """
    offshore = maritime.dissolve(by="country")[["geometry"]].rename_axis("name")
    offshore["geometry"] = offshore.geometry.apply(
        _simplify_polys, minarea=0.1, filterremote=False
    )
    return offshore


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_offshore_shapes")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    maritime = read_geo_boundaries(
        snakemake.input.shapes, snakemake.params.countries, "maritime"
    )
    offshore_shapes = build_offshore_shapes(maritime)
    offshore_shapes.reset_index().to_file(snakemake.output.offshore_shapes)
