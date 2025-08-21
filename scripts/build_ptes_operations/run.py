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
                dynamic_capacity:
                discharge_boosting_required:
                charge_boosting_required:
                temperature_profile
                max_top_temperature:
                min_bottom_temperature:

Inputs
------
- `resources/<run_name>/forward_temperature.nc`
    Forward temperature profiles for the district heating networks.
- `resources/<run_name>/central_heating_return_temperature_profiles.nc`:
    Return temperature profiles for the district heating networks.
- `resources/<run_name>/ptes_temperature_boost_ratio_profiles.nc`
    Ratio of PTES charge that requires additional heating due to temperature differences.

Outputs
-------
- `resources/<run_name>/ptes_top_temperature_profiles.nc`
    Clipped PTES top temperature profile (in °C).
- `resources/<run_name>/ptes_e_max_pu_profiles.nc`
    Normalized PTES capacity profiles.
- `ptes_temperature_boost_ratio_profiles` (netCDF):
    Charging temperature boost ratio time series.
- `ptes_forward_temperature_boost_ratio_profiles` (netCDF):
    Forward flow temperature boost ratio time series.

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
    TesTemperatureMode,
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

    if (snakemake.params.charge_boosting_required and
            TesTemperatureMode(snakemake.params.ptes_temperature_profile) is TesTemperatureMode.DYNAMIC):
        raise ValueError(
            "Charger boosting cannot be used with 'dynamic' temperature profile"
        )

    # Load temperature profiles
    logger.info(
        f"Loading district heating temperature profiles and approximating PTES temperatures"
    )
    logger.info(
        f"PTES configuration: temperature_profile={snakemake.params.ptes_temperature_profile}, "
        f"charge_boosting_required={snakemake.params.charge_boosting_required}, "
        f"discharge_boosting_required={snakemake.params.discharge_boosting_required}, "
        f"dynamic_capacity={snakemake.params.dynamic_capacity}"
    )

    # Initialize unified PTES temperature class
    ptes_temperature_approximator = PtesTemperatureApproximator(
        forward_temperature=xr.open_dataarray(
            snakemake.input.central_heating_forward_temperature_profiles
        ),
        return_temperature=xr.open_dataarray(
            snakemake.input.central_heating_return_temperature_profiles
        ),
        max_top_temperature=snakemake.params.max_ptes_top_temperature,
        min_bottom_temperature=snakemake.params.min_ptes_bottom_temperature,
        temperature_profile=TesTemperatureMode(snakemake.params.ptes_temperature_profile),
        charge_boosting_required=snakemake.params.charge_boosting_required,
        discharge_boosting_required=snakemake.params.discharge_boosting_required,
        dynamic_capacity=snakemake.params.dynamic_capacity,
    )

    ptes_temperature_approximator.top_temperature.to_netcdf(
        snakemake.output.ptes_top_temperature_profile
    )

    ptes_temperature_approximator.e_max_pu.to_netcdf(
        snakemake.output.ptes_e_max_pu_profile
    )
    
    ptes_temperature_approximator.boost_per_discharge.to_netcdf(
        snakemake.output.boost_per_discharge_profile
    )

    ptes_temperature_approximator.boost_per_charge.to_netcdf(
        snakemake.output.boost_per_charge_profile
    )
