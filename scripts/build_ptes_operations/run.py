# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build operational profiles of pit thermal energy storage (PTES) in district heating networks.

The district heating forward temperature is taken as the PTES top temperature
and clipped at the configured maximum. A binary direct-utilisation profile is
1 where the forward temperature does not exceed this maximum, so the storage
can serve the network directly, and 0 where supplemental heating, for example
by a heat pump, is required. The usable capacity `e_max_pu` scales linearly
with the difference between top and return temperature, following Sorknaes
(2018), normalised by the difference between the configured maximum top and
minimum bottom temperatures.

References
----------
- Sorknaes (2018), [Simulation method for a pit seasonal thermal energy storage system with a heat pump in a district heating system](https://doi.org/10.1016/j.energy.2018.03.152)
"""

import logging

import xarray as xr
from _helpers import set_scenario_config

from scripts.build_ptes_operations.ptes_temperature_approximator import (
    PtesTemperatureApproximator,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_ptes_operations",
            horizon="2030",
        )

    set_scenario_config(snakemake)

    # Load temperature profiles
    logger.info(
        "Loading district heating temperature profiles and constructing PTES temperature approximator"
    )
    # Initialize unified PTES temperature class
    ptes_temperature_approximator = PtesTemperatureApproximator(
        forward_temperature=xr.open_dataarray(
            snakemake.input.central_heating_forward_temperature_profiles
        ),
        return_temperature=xr.open_dataarray(
            snakemake.input.central_heating_return_temperature_profiles
        ),
        max_ptes_top_temperature=snakemake.params.max_ptes_top_temperature,
        min_ptes_bottom_temperature=snakemake.params.min_ptes_bottom_temperature,
    )

    # Get PTES clipped top temperature profiles
    logger.info(
        f"Saving TES top temperature profile to {snakemake.output.ptes_top_temperature_profiles}"
    )
    ptes_temperature_approximator.top_temperature.to_netcdf(
        snakemake.output.ptes_top_temperature_profiles
    )

    # if snakemake.params.enable_supplemental_heating:
    # Get PTES supplemental heating profiles
    logger.info(
        f"Saving PTES direct utilisation profile to {snakemake.output.ptes_direct_utilisation_profiles}"
    )
    ptes_temperature_approximator.direct_utilisation_profile.to_netcdf(
        snakemake.output.ptes_direct_utilisation_profiles
    )

    # if snakemake.params.enable_dynamic_capacity:
    logger.info("Calculating dynamic PTES capacity profiles")

    # Get PTES capacity profiles
    logger.info(
        f"Saving PTES capacity profiles to {snakemake.output.ptes_e_max_pu_profiles}"
    )
    ptes_temperature_approximator.e_max_pu.to_netcdf(
        snakemake.output.ptes_e_max_pu_profiles
    )
