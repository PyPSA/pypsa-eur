# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build solar rooftop potentials for all clustered model regions per resource class.
"""

import geopandas as gpd
import pandas as pd
import xarray as xr

from scripts._helpers import load_cutout, set_scenario_config

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_clustered_solar_rooftop_potentials",
            clusters=5,
            configfiles="config/test/config.myopic.yaml",
        )

    set_scenario_config(snakemake)

    cutout = load_cutout(snakemake.input.cutout)

    class_regions = gpd.read_file(snakemake.input.class_regions).set_index(
        ["bus", "bin"]
    )

    I = cutout.indicatormatrix(class_regions)  # noqa: E741

    with xr.open_dataarray(snakemake.input.pop_layout) as pop_layout:
        pop = I.dot(pop_layout.stack(spatial=("y", "x")))

    # add max solar rooftop potential assuming 0.1 kW/m2 and 20 m2/person,
    # i.e. 2 kW/person (population data is in thousands of people) so we get MW
    potentials = 0.1 * 20 * pd.Series(pop, index=class_regions.index)

    potentials.to_csv(snakemake.output.potentials)
