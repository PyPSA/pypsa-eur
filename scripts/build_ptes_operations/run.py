# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build PTES (Pit Thermal Energy Storage) operational profiles.

This script calculates temperature and capacity profiles for pit thermal energy
storage systems integrated with district heating networks. It determines:

1. **Top temperature profiles**: The operational top layer temperature, either
   following the district heating forward temperature (clipped to design limits)
   or a constant value.

2. **Capacity profiles (e_max_pu)**: Normalized storage capacity based on the
   temperature difference between top and bottom layers, relative to the design
   temperature difference. This captures how storage capacity varies when
   operating temperatures deviate from design conditions.

3. **Discharge boosting profiles**: When the PTES top temperature is below the
   required forward temperature, additional heating (boosting) is needed during
   discharge. This profile quantifies the ratio of boost energy to stored energy.

The outputs are used by ``prepare_sector_network.py`` to configure PTES storage
components, charger/discharger links, and optional resistive boosting infrastructure.

Relevant Settings
-----------------
.. code:: yaml

    sector:
        district_heating:
            ptes:
                enable: true
                temperature_dependent_capacity: false # if true, e_max_pu varies with temperature difference (static but scaled if top/bottom are constant)
                charge_boosting_required: false # currently not supported
                discharge_resistive_boosting: false # if true, adds resistive boosting link for discharge boosting and disables heat pump boosting
                top_temperature: 90  # or "forward" for dynamic
                bottom_temperature: 35  # or "return" for dynamic
                design_top_temperature: 90 # used to compute design temperature difference for e_max_pu if temperature_dependent_capacity is true
                design_bottom_temperature: 35 # used to compute design temperature difference for e_max_pu if temperature_dependent_capacity is true

Inputs
------
- ``resources/<run_name>/central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Forward temperature profiles for district heating networks (°C).
- ``resources/<run_name>/central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Return temperature profiles for district heating networks (°C).

Outputs
-------
- ``resources/<run_name>/temp_ptes_top_profiles_base_s_{clusters}_{planning_horizons}.nc``
    PTES top layer temperature profile (°C), clipped to design_top_temperature.
- ``resources/<run_name>/ptes_e_max_pu_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Normalized PTES capacity profiles (0-1 p.u.).
- ``resources/<run_name>/ptes_boost_per_discharge_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Discharge boosting ratio profiles. Values represent the ratio of boost energy
    to discharge energy needed when T_forward > T_top. Only used when resistive
    boosting is enabled.

References
----------
Sorknæs, P. (2018). "Simulation method for a pit seasonal thermal energy storage
system with a heat pump in a district heating system." Energy, Volume 152.
https://doi.org/10.1016/j.energy.2018.03.152
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

    logger.info(
        "Loading district heating temperature profiles and calculating PTES operational profiles"
    )
    logger.info(f"PTES configuration: {snakemake.params}")

    # Discharge boosting profiles are calculated when resistive boosting is enabled.
    # For heat pump-based boosting, add "ptes" to central heating heat sources instead and disable resistive boosting.
    discharge_boosting_required: bool = snakemake.params.discharge_resistive_boosting

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
        design_bottom_temperature=snakemake.params.design_bottom_temperature,
        design_top_temperature=snakemake.params.design_top_temperature,
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
