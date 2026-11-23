# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build operational profiles of pit thermal energy storage (PTES) in district heating networks.

The PTES top temperature is either constant or follows the district heating
forward temperature, clipped at the design top temperature. The bottom
temperature is either constant or follows the return temperature. The usable
capacity `e_max_pu` scales linearly with the difference between top and
bottom temperature, following Sorknaes (2018), normalised by the design
temperature difference. Where the forward temperature exceeds the top
temperature, the boost-per-discharge profile gives the heat that a resistive
heater must add per unit of discharged heat.

References
----------
- Sorknaes (2018), [Simulation method for a pit seasonal thermal energy storage system with a heat pump in a district heating system](https://doi.org/10.1016/j.energy.2018.03.152)
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
            horizon="2030",
        )

    set_scenario_config(snakemake)

    logger.info(
        "Loading district heating temperature profiles and calculating PTES operational profiles"
    )

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
