# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Approximate the top temperature of the pit thermal energy storage (PTES), ensuring that the temperature does not
exceed the operational limit.

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
                temperature_dependent_capacity:
                charge_boosting_required:
                discharge_resistive_boosting:
                max_top_temperature:
                min_bottom_temperature:

Inputs
------
- `resources/<run_name>/central_heating_forward_temperature_profiles.nc`
    Forward temperature profiles for the district heating networks.
- `resources/<run_name>/central_heating_return_temperature_profiles.nc`:
    Return temperature profiles for the district heating networks.

Outputs
-------
- `resources/<run_name>/ptes_top_temperature_profile.nc`
    Clipped PTES top temperature profile (in °C).
- `resources/<run_name>/ptes_e_max_pu_profile.nc`
    Normalized PTES capacity profiles.
- `resources/<run_name>/boost_per_discharge_profile.nc` (conditional):
    Discharge temperature boost ratio time series (only if discharge_resistive_boosting is enabled).

Source
------
Sorknæs, P. 2018. "Simulation method for a pit seasonal thermal energy storage system with a heat pump in a district heating system", Energy, Volume 152, https://doi.org/10.1016/j.energy.2018.03.152.
Approximate thermal energy storage (TES) top temperature and identify need for supplemental heating.
"""

import logging

import xarray as xr

from scripts._helpers import set_scenario_config
from scripts.build_ptes_operations.ptes_temperature_approximator import (
    PtesTemperatureApproximator,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_ptes_operations",
            clusters=5,
            planning_horizons="2030",
        )

    set_scenario_config(snakemake)

    # Load temperature profiles
    logger.info(
        "Loading district heating temperature profiles and approximating PTES temperatures"
    )
    logger.info(f"PTES configuration: {snakemake.params}")

    # Discharge boosting is "required" only if boosting via resistive heaters is enabled
    # if you'd like to model boosting via heat pumps, add "ptes" to central heating heat sources and set temperatures accordingly
    discharge_boosting_required: bool = snakemake.params.discharge_resistive_boosting

    # Initialize unified PTES temperature class
    ptes_temperature_approximator = PtesTemperatureApproximator(
        forward_temperature=xr.open_dataarray(
            snakemake.input.central_heating_forward_temperature_profiles
        ),
        return_temperature=xr.open_dataarray(
            snakemake.input.central_heating_return_temperature_profiles
        ),
        top_temperature=snakemake.params.top_temperature,
        bottom_temperature=snakemake.params.bottom_temperature,
        charge_boosting_required=snakemake.params.charge_boosting_required,
        discharge_boosting_required=discharge_boosting_required,
        temperature_dependent_capacity=snakemake.params.temperature_dependent_capacity,
    )

    ptes_temperature_approximator.top_temperature_profile.to_netcdf(
        snakemake.output.ptes_top_temperature_profiles
    )

    ptes_temperature_approximator.e_max_pu.to_netcdf(
        snakemake.output.ptes_e_max_pu_profiles
    )

    ptes_temperature_approximator.boost_per_discharge.to_netcdf(
        snakemake.output.ptes_boost_per_discharge_profiles
    )
