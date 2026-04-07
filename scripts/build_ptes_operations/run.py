# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build PTES (Pit Thermal Energy Storage) operational profiles.

This script pre-computes temperature-dependent parameters for pit thermal
energy storage systems integrated with district heating networks, using the
unified ``PtesApproximator`` class.

Outputs are consumed by ``prepare_sector_network.py`` and by the COP / heat
source utilisation profile builders.

Relevant Settings
-----------------
.. code:: yaml

    sector:
        district_heating:
            ptes:
                enable: true
                top_temperature: 90
                bottom_temperature: 35
                design_top_temperature: 90
                design_bottom_temperature: 35
                discharge_resistive_boosting: false
                layered:
                    num_layers: 1

Inputs
------
- ``resources/<run_name>/central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``
- ``resources/<run_name>/central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``

Outputs
-------
- ``resources/<run_name>/ptes_layered_params_base_s_{clusters}_{planning_horizons}.nc``
    Pre-computed parameter dataset for the PTES model.

References
----------
Sorknæs, P. (2018). "Simulation method for a pit seasonal thermal energy storage
system with a heat pump in a district heating system." Energy, Volume 152.
https://doi.org/10.1016/j.energy.2018.03.152
"""

import logging

import xarray as xr

from scripts._helpers import set_scenario_config
from scripts.build_ptes_operations.ptes_approximator import PtesApproximator

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

    logger.info("Building PTES operational profiles")
    logger.info(f"PTES configuration: {snakemake.params}")

    forward_temperature = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    return_temperature = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    # Resolve top/bottom temperature to constant floats
    top_temp = float(snakemake.params.top_temperature)
    bottom_temp = float(snakemake.params.bottom_temperature)

    layered_config = snakemake.params.layered

    approximator = PtesApproximator(
        forward_temperature=forward_temperature,
        return_temperature=return_temperature,
        top_temperature=top_temp,
        bottom_temperature=bottom_temp,
        num_layers=layered_config["num_layers"],
        design_top_temperature=snakemake.params.design_top_temperature,
        design_bottom_temperature=snakemake.params.design_bottom_temperature,
        design_standing_losses=0.0,  # placeholder; actual value set from costs data
        interlayer_heat_transfer_coefficient=snakemake.params.interlayer_heat_transfer_coefficient,
    )

    # Write single dataset with all pre-computed PTES parameters
    approximator.to_dataset().to_netcdf(snakemake.output.ptes_operations)
