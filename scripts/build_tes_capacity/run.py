# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Calculate dynamic pit thermal energy storage (PTES) capacity profiles based on
district heating forward and return flow temperatures. The linear relation between
temperature difference and capacity is taken from Sorknaes (2018).

The capacity of thermal energy storage systems varies with the temperature difference
between the forward and return flows in district heating networks assuming a direct
integration of the storage. This script calculates normalized capacity factors (e_max_pu)
for PTES systems based on these temperature differences.

Source
------
Sorknæs, P. 2018. "Simulation method for a pit seasonal thermal energy storage system with a heat pump in a district heating system", Energy, Volume 152, https://doi.org/10.1016/j.energy.2018.03.152.

Relevant Settings
-----------------

.. code:: yaml
    sector:
        district_heating:
            dynamic_ptes_capacity: true

Inputs
------
- `resources/<run_name>/central_heating_forward_temperature_profiles.nc`: Forward flow temperature profiles
- `resources/<run_name>/central_heating_return_temperature_profiles.nc`: Return flow temperature profiles

Outputs
-------
- `resources/<run_name>/ptes_e_max_pu_profiles.nc`: Normalized PTES capacity profiles
"""

import logging

import xarray as xr

from scripts._helpers import set_scenario_config
from scripts.build_tes_capacity.tes_capacity_approximator import TesCapacityApproximator

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tes_parameters",
            clusters=48,
            planning_horizons="2050",
        )

    set_scenario_config(snakemake)

    # Load temperature profiles
    logger.info("Loading district heating temperature profiles")
    t_forward = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    t_return = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    # Define operational limits for PTES
    max_top_temperature = snakemake.params.get("max_top_temperature")
    min_bottom_temperature = snakemake.params.get("min_bottom_temperature")

    logger.info(
        f"Calculating PTES capacity profiles with max temperature {max_top_temperature}°C"
    )

    # Create TES capacity approximator
    tes_capacity_approximator = TesCapacityApproximator(
        top_temperature=t_forward,
        bottom_temperature=t_return,
        max_top_temperature=max_top_temperature,
        min_bottom_temperature=min_bottom_temperature,
    )

    # Calculate e_max_pu
    e_max_pu = tes_capacity_approximator.calculate_e_max_pu()

    # Save output
    logger.info(
        f"Saving PTES capacity profiles to {snakemake.output.ptes_e_max_pu_profiles}"
    )
    e_max_pu.to_netcdf(snakemake.output.ptes_e_max_pu_profiles)
