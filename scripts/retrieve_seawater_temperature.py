# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""

Description
-----------
"""

import logging

import copernicusmarine

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "retrieve_seawater_temperature",
            clusters="39",
            opts="",
            ll="vopt",
            sector_opts="",
            planning_horizons=2050,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    _ = copernicusmarine.subset(
        dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
        start_datetime=snakemake.wildcards.start_snapshot,
        end_datetime=snakemake.wildcards.end_snapshot,
        minimum_longitude=-12,
        maximum_longitude=41,
        minimum_latitude=34,
        maximum_latitude=72,
        variables=["thetao"],
        minimum_depth=5,
        maximum_depth=15,
        output_filename=snakemake.output.seawater_temperature,
    )
