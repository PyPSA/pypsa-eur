# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Aggregate a per-bus land-eligibility availability matrix (as produced by
``determine_availability_matrix.py`` / ``determine_rock_weathering_availability_matrix.py``)
into eligible land area per network node, in km².

Shared by biochar, afforestation and rock weathering: the per-bus fractional
eligibility (0-1 per cutout grid cell) is multiplied by the cutout grid cell
area and summed per bus, recovering the same eligible-area-in-km² number the
previous per-technology scripts (``build_corine_potentials.py``,
``build_EW_potentials.py``) each computed independently via
``atlite.gis.shape_availability``.

Output: CSV with columns ``area [sqkm]`` (total node area) and
``potential [sqkm]`` (eligible area after exclusions), indexed by node name.
The downstream consumers (``add_biochar``, ``build_afforestation_potentials``,
``add_rock_weathering``) apply their own technology-specific rate (t/km²,
growth rate, etc.) to the ``potential [sqkm]`` column.
"""

import logging

import geopandas as gpd
import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging, load_cutout, set_scenario_config

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_available_land", clusters="39", technology="biochar"
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    technology = snakemake.wildcards.technology

    availability = xr.open_dataarray(snakemake.input.availability_matrix)
    cutout = load_cutout(snakemake.input.cutout)
    regions = (
        gpd.read_file(snakemake.input.regions).set_index("name").rename_axis("bus")
    )

    cell_area = cutout.grid.to_crs(3035).area / 1e6
    cell_area = xr.DataArray(
        cell_area.values.reshape(cutout.shape),
        [cutout.coords["y"], cutout.coords["x"]],
    )

    potential = (availability * cell_area).sum(dim=["x", "y"]).to_pandas()
    area = regions.geometry.to_crs(3035).area / 1e6

    df = pd.DataFrame({"area [sqkm]": area, "potential [sqkm]": potential})
    df.index.name = "node"

    logger.info(
        "Total %s potential: %.0f sqkm (of %.0f sqkm total node area)",
        technology,
        df["potential [sqkm]"].sum(),
        df["area [sqkm]"].sum(),
    )

    df.to_csv(snakemake.output.csv_file)
