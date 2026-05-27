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
                    layer_temperatures: [90]

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
from scripts.build_cop_profiles.central_heating_cop_approximator import (
    CentralHeatingCopApproximator,
)
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

    approximator = PtesApproximator(
        forward_temperature=forward_temperature,
        return_temperature=return_temperature,
        layer_temperatures=snakemake.params.layer_temperatures,
        design_top_temperature=snakemake.params.design_top_temperature,
        design_bottom_temperature=snakemake.params.design_bottom_temperature,
        design_standing_losses=0.0,  # placeholder; actual value set from costs data
        interlayer_heat_transfer_coefficient=snakemake.params.interlayer_heat_transfer_coefficient,
        cop=xr.open_dataarray(snakemake.input.cop_profiles),
        conservative_return_layer=snakemake.params.conservative_return_layer,
    )

    # Warn if the chosen layer temperatures place the return-level layer above the
    # actual return temperature, which lets direct discharge create energy. Stay
    # conservative by configuring a layer at/below the return temperature.
    free_lunch_warnings = approximator.discharge_free_lunch_warnings()
    if free_lunch_warnings:
        logger.warning(
            "PTES discretisation may create energy on discharge for %d node(s): "
            "the return-level layer sits ABOVE the district-heating return "
            "temperature. Add a layer at or below the return temperature in "
            "sector.district_heating.ptes.layered.layer_temperatures to stay "
            "conservative.\n%s",
            len(free_lunch_warnings),
            "\n".join(free_lunch_warnings),
        )

    # Write single dataset with all pre-computed PTES parameters, including the
    # booster COP recomputed at the actual evaporator outlet (the reinjection
    # layer depth) so the electricity/source split reflects the real cooling depth.
    dataset = approximator.to_dataset()
    hp_cop_params = snakemake.params.heat_pump_cop_approximation_central_heating
    dataset["booster_cop"] = approximator.booster_cop(
        CentralHeatingCopApproximator,
        refrigerant=hp_cop_params["refrigerant"],
        delta_t_pinch_point=hp_cop_params[
            "heat_exchanger_pinch_point_temperature_difference"
        ],
        isentropic_compressor_efficiency=hp_cop_params[
            "isentropic_compressor_efficiency"
        ],
        heat_loss=hp_cop_params["heat_loss"],
        min_delta_t_lift=hp_cop_params["min_delta_t_lift"],
    )
    # Clear inherited encodings (the booster COP carries the temperature profile's
    # time encoding, which can conflict with the dataset's and overflow on decode)
    # and write the time axis with an explicit, safe reference so it round-trips.
    for name in list(dataset.variables):
        dataset[name].encoding.clear()
    encoding = {}
    if "time" in dataset.coords:
        encoding["time"] = {
            "units": "hours since 1900-01-01 00:00:00",
            "calendar": "proleptic_gregorian",
            "dtype": "int64",
        }
    dataset.to_netcdf(snakemake.output.ptes_operations, encoding=encoding)
