# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Generate placeholder solar rooftop potentials for reduced test runs."""

import geopandas as gpd

from scripts._helpers import configure_logging, set_scenario_config

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_solar_rooftop_potentials")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    class_regions = gpd.read_file(snakemake.input.class_regions)

    potentials = (
        class_regions[["bus", "bin"]]
        .drop_duplicates()
        .sort_values(["bus", "bin"])
        .reset_index(drop=True)
    )
    potentials["potential_mwp"] = 0.0

    potentials.to_csv(snakemake.output.potentials, index=False)
