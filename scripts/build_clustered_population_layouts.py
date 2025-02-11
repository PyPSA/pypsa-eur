# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build population layouts for all clustered model regions as total as well as
split by urban and rural population.
"""

import logging

import atlite
import geopandas as gpd
import pandas as pd
import xarray as xr
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_clustered_population_layouts", clusters=48)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    cutout = atlite.Cutout(snakemake.input.cutout)

    clustered_regions = (
        gpd.read_file(snakemake.input.regions_onshore).set_index("name").buffer(0)
    )

    I = cutout.indicatormatrix(clustered_regions)  # noqa: E741

    pop = {}
    for item in ["total", "urban", "rural"]:
        pop_layout = xr.open_dataarray(snakemake.input[f"pop_layout_{item}"])
        pop[item] = I.dot(pop_layout.stack(spatial=("y", "x")))

    pop = pd.DataFrame(pop, index=clustered_regions.index)

    pop["ct"] = pop.index.str[:2]
    country_population = pop.total.groupby(pop.ct).sum()
    pop["fraction"] = pop.total / pop.ct.map(country_population)

    pop.to_csv(snakemake.output.clustered_pop_layout)
