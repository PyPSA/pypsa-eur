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
                dynamic_ptes_capacity:
                supplemental_heating:
                    enable:
                max_top_temperature:

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
- `resources/<run_name>/ptes_supplemental_heating_required.nc`
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

from scripts.build_ptes_operations.ptes_temperature_approximator import (
    PtesTemperatureApproximator,
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
