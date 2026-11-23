# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build the land area usable for a land-based carbon dioxide removal option per region.

The availability matrix gives the eligible share of each cutout grid cell per
region. It is multiplied with the grid cell area and summed per region.
Optionally, grid cells with an annual mean air temperature below a threshold
are excluded. Only the configured maximum share of the eligible land is used.
"""

import logging

import xarray as xr

from scripts._helpers import configure_logging, load_cutout, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_carbon_dioxide_removal_land", technology="afforestation"
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    availability = xr.open_dataarray(snakemake.input.availability_matrix)
    cutout = load_cutout(snakemake.input.cutout)

    config = snakemake.params.technology
    if "min_mean_temp" in config:
        mean_temp = cutout.data["temperature"].mean("time") - 273.15
        availability = availability.where(mean_temp >= config["min_mean_temp"], 0)

    cell_area = xr.DataArray(
        (cutout.grid.to_crs(3035).area / 1e6).values.reshape(cutout.shape),
        [cutout.coords["y"], cutout.coords["x"]],
    )
    land = (availability * cell_area).sum(["x", "y"]).to_pandas()
    land *= config["max_land_usage"]

    logger.info(
        f"Usable land for {snakemake.wildcards.technology}: {land.sum():.0f} sqkm"
    )
    land.rename_axis("name").rename("area").to_csv(snakemake.output[0])
