# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Approximate the top temperature of the pit thermal energy storage (PTES), ensuring that the temperature does not
exceed the operational limit.

Determine whether supplemental heating is needed. A binary indicator is generated:
    - 1: The forward temperature is less than or equal to the TES maximum; direct usage is possible.
    - 0: The forward temperature exceeds the TES maximum; supplemental heating (e.g., via a heat pump) is required.

Calculate dynamic PTES capacity profiles based on district heating forward and return flow temperatures.
The linear relation between temperature difference and capacity is taken from Sorknaes (2018).

The capacity of thermal energy storage systems varies with the temperature difference
between the forward and return flows in district heating networks assuming a direct
integration of the storage. This script calculates normalized capacity factors (e_max_pu)
for PTES systems based on these temperature differences.

Relevant Settings
-----------------
.. code:: yaml
    sector
        district_heating:
            ptes:
                dynamic_ptes_capacity: true
                supplemental_heating:
                    enable: true
                max_top_temperature: 90  # Maximum PTES temperature (°C) usable without supplemental heating

Inputs
------
- `resources/<run_name>/forward_temperature.nc`
    Forward temperature profiles for the district heating networks.
- `resources/<run_name>/central_heating_return_temperature_profiles.nc`:
    Return temperature profiles for the district heating networks.

Outputs
-------
- `resources/<run_name>/ptes_top_temperature_profiles.nc`
    Clipped PTES top temperature profile (in °C).
- `resources/<run_name>/ptes_supplemental_heating_profiles.nc`
    Binary indicator for additional heating (1 = direct PTES use, 0 = supplemental heating required).
- `resources/<run_name>/ptes_e_max_pu_profiles.nc`
    Normalized PTES capacity profiles.

Source
------
Sorknæs, P. 2018. "Simulation method for a pit seasonal thermal energy storage system with a heat pump in a district heating system", Energy, Volume 152, https://doi.org/10.1016/j.energy.2018.03.152.
Approximate thermal energy storage (TES) top temperature and identify need for supplemental heating.
"""

import logging

import xarray as xr
from _helpers import set_scenario_config

from scripts.build_ptes_operations.PTESCapacityApproximator import (
    PTESCapacityApproximator,
)
from scripts.build_ptes_operations.PTESSupplementalHeatingApproximator import (
    PTESSupplementalHeatingApproximator,
)
from scripts.build_ptes_operations.PTESTopTemperatureProfileApproximator import (
    PTESTopTemperatureProfileApproximator,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_ptes_operations",
            clusters=5,
            planning_horizons="2030",
        )

    set_scenario_config(snakemake)

    # Load temperature profiles
    logger.info("Loading district heating temperature profiles")

    forward_temp = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    return_temp = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    # Load parameter for max. and bottom PTES temperature.
    max_ptes_top_temperature = snakemake.params.max_ptes_top_temperature
    min_ptes_bottom_temperature = snakemake.params.min_ptes_bottom_temperature
    logger.info(
        f"Using maximum PTES direct usage temperature: {max_ptes_top_temperature}°C"
    )

    # Initialize TES top temperature approximator.
    ptes_top_temperature_profile = PTESTopTemperatureProfileApproximator(
        forward_temperature=forward_temp,
        max_ptes_top_temperature=max_ptes_top_temperature,
    ).clipped_top_temperature

    logger.info(
        f"Saving TES top temperature profile to {snakemake.output.ptes_top_temperature_profiles}"
    )
    ptes_top_temperature_profile.to_netcdf(
        snakemake.output.ptes_top_temperature_profiles
    )

    if snakemake.params.enable_ptes_supplemental_heating_approximatior:
        # Initialize supplemental heating approximator.
        ptes_supplemental_heating_profiles = PTESSupplementalHeatingApproximator(
            forward_temperature=forward_temp,
            max_ptes_top_temperature=max_ptes_top_temperature,
        ).determine_ptes_usage()

        logger.info(
            f"Saving supplemental heating profile to {snakemake.output.ptes_supplemental_heating_profiles}"
        )
        ptes_supplemental_heating_profiles.to_netcdf(
            snakemake.output.ptes_supplemental_heating_profiles
        )

    if snakemake.params.enable_ptes_capacity_approximatior:
        logger.info(
            f"Calculating PTES capacity profiles with max temperature {max_ptes_top_temperature}°C"
        )

        # Create TES capacity approximator
        ptes_e_max_pu = PTESCapacityApproximator(
            top_temperature=ptes_top_temperature_profile,
            bottom_temperature=return_temp,
            max_ptes_top_temperature=max_ptes_top_temperature,
            min_ptes_bottom_temperature=min_ptes_bottom_temperature,
        ).calculate_e_max_pu()

        logger.info(
            f"Saving PTES capacity profiles to {snakemake.output.ptes_e_max_pu_profiles}"
        )
        ptes_e_max_pu.to_netcdf(snakemake.output.ptes_e_max_pu_profiles)
