# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Approximate the top temperature of the TES system, ensuring that the temperature does not exceed the operational limit.

Determine whether supplemental heating is needed. A binary indicator is generated:
    - 1: The forward temperature is less than or equal to the TES maximum; direct usage is possible.
    - 0: The forward temperature exceeds the TES maximum; supplemental heating (e.g., via a heat pump) is required.

Calculate dynamic pit thermal energy storage (PTES) capacity profiles based on
district heating forward and return flow temperatures. The linear relation between
temperature difference and capacity is taken from Sorknaes (2018).

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
- `resources/<run_name>/tes_top_temperature_profile.nc`
    Clipped TES top temperature profile (in °C).
- `resources/<run_name>/additional_heating_indicator.nc`
    Binary indicator for additional heating (1 = direct TES use, 0 = supplemental heating required).
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

from scripts.build_tes_operation.build_tes_top_temperature_profile import (
    BuildTesTopTemperature,
)
from scripts.build_tes_operation.tes_capacity_approximator import (
    TesCapacityApproximator,
)
from scripts.build_tes_operation.tes_supplemental_heating_approximator import (
    TESSupplementalHeatingApproximator,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tes_operation",
            clusters=6,
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
    max_ptes_temperature = snakemake.params.max_ptes_temperature
    min_bottom_temperature = snakemake.params.get("min_bottom_temperature")
    logger.info(
        f"Using maximum PTES direct usage temperature: {max_ptes_temperature}°C"
    )

    # Initialize TES top temperature approximator.
    tes_temp_builder = BuildTesTopTemperature(
        forward_temperature_celsius=forward_temp,
        max_ptes_temperature=max_ptes_temperature,
    )

    # Approximate PTES top temperature profile
    ptes_top_temperature = tes_temp_builder.clipped_top_temperature

    # Initialize supplemental heating approximator.
    supplemental_heating_approximator = TESSupplementalHeatingApproximator(
        forward_temperature_celsius=forward_temp,
        max_ptes_temperature=max_ptes_temperature,
    )

    logger.info(
        f"Calculating PTES capacity profiles with max temperature {max_ptes_temperature}°C"
    )

    # Create TES capacity approximator
    tes_capacity_approximator = TesCapacityApproximator(
        top_temperature=ptes_top_temperature,
        bottom_temperature=return_temp,
        max_top_temperature=max_ptes_temperature,
        min_bottom_temperature=min_bottom_temperature,
    )

    # Calculate e_max_pu
    e_max_pu = tes_capacity_approximator.calculate_e_max_pu()

    # Save output
    logger.info(
        f"Saving PTES capacity profiles to {snakemake.output.ptes_e_max_pu_profiles}"
    )
    e_max_pu.to_netcdf(snakemake.output.ptes_e_max_pu_profiles)

    supplemental_heating_profile = (
        supplemental_heating_approximator.determine_ptes_usage()
    )

    # Save output
    logger.info(
        f"Saving TES top temperature profile to {snakemake.output.tes_top_temperature_profile}"
    )
    ptes_top_temperature.to_netcdf(snakemake.output.tes_top_temperature_profile)

    # Save output
    logger.info(
        f"Saving supplemental heating profile to {snakemake.output.tes_supplemental_heating_profile}"
    )
    supplemental_heating_profile.to_netcdf(
        snakemake.output.tes_supplemental_heating_profile
    )
