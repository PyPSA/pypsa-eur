# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build availability profiles for direct heat source utilisation (1 in regions and time steps where heat source can be utilised, 0 otherwise).
When direct utilisation is possible, heat pump COPs are set to zero (c.f. `build_cop_profiles`).

Inputs
------
- `resources/<run_name>/central_heating_forward_temperatures_base_s_{clusters}_{planning_horizons}.nc`: Central heating forward temperature profiles

Outputs
-------
- `resources/<run_name>/direct_heat_source_utilisation_profiles_base_s_{clusters}_{planning_horizons}.nc`: Direct heat source utilisation profiles
"""

import logging

import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config
from scripts.definitions.heat_source import HeatSource

logger = logging.getLogger(__name__)


def get_source_temperature(
    snakemake_params: dict, snakemake_input: dict, heat_source_name: str
) -> float | xr.DataArray:
    heat_source = HeatSource(heat_source_name)
    if heat_source.has_constant_temperature:
        try:
            return snakemake_params[f"constant_temperature_{heat_source_name}"]
        except KeyError:
            raise ValueError(
                f"Constant temperature for heat source {heat_source_name} not specified in parameters."
            )

    else:
        if f"temp_{heat_source_name}" not in snakemake_input.keys():
            raise ValueError(
                f"Missing input temperature for heat source {heat_source_name}."
            )
        return xr.open_dataarray(snakemake_input[f"temp_{heat_source_name}"])


def get_direct_utilisation_profile(
    source_temperature: float | xr.DataArray, forward_temperature: xr.DataArray
) -> xr.DataArray | float:
    """
    Get the direct heat source utilisation profile.

    Args:
    ----
    source_temperature: float | xr.DataArray
        The constant temperature of the heat source in degrees Celsius. If `xarray`, indexed by `time` and `region`. If a float, it is broadcasted to the shape of `forward_temperature`.
    forward_temperature: xr.DataArray
        The central heating forward temperature profiles. If `xarray`, indexed by `time` and `region`. If a float, it is broadcasted to the shape of `return_temperature`.

    Returns:
    -------
    xr.DataArray | float
        The direct heat source utilisation profile.

    """
    return xr.where(source_temperature >= forward_temperature, 1.0, 0.0)


def get_preheater_utilisation_profile(
    source_temperature: float | xr.DataArray,
    forward_temperature: xr.DataArray,
    return_temperature: xr.DataArray,
    heat_source_cooling: float,
) -> xr.DataArray | float:
    """
    Get the direct heat source utilisation profile.

    Args:
    ----
    source_temperature: float | xr.DataArray
        The constant temperature of the heat source in degrees Celsius. If `xarray`, indexed by `time` and `region`. If a float, it is broadcasted to the shape of `forward_temperature`.
    forward_temperature: xr.DataArray
        The central heating forward temperature profiles. If `xarray`, indexed by `time` and `region`. If a float, it is broadcasted to the shape of `return_temperature`.
    return_temperature: xr.DataArray
        The central heating return temperature profiles. If `xarray`, indexed by `time` and `region`.

    Returns:
    -------
    xr.DataArray | float
        The direct heat source utilisation profile.

    """
    return xr.where(
        (source_temperature < forward_temperature)
        * (source_temperature > return_temperature),
        heat_source_cooling
        / (source_temperature - return_temperature + heat_source_cooling),
        0.0,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            clusters=48,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    heat_sources: list[str] = snakemake.params.heat_sources

    central_heating_forward_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    central_heating_return_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    xr.concat(
        [
            get_direct_utilisation_profile(
                source_temperature=get_source_temperature(
                    snakemake_params=snakemake.params,
                    snakemake_input=snakemake.input,
                    heat_source_name=heat_source_key,
                ),
                forward_temperature=central_heating_forward_temperature,
            ).assign_coords(heat_source=heat_source_key)
            for heat_source_key in heat_sources
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.heat_source_direct_utilisation_profiles)

    xr.concat(
        [
            get_preheater_utilisation_profile(
                source_temperature=get_source_temperature(
                    heat_source_name=heat_source_key,
                    snakemake_params=snakemake.params,
                    snakemake_input=snakemake.input,
                ),
                forward_temperature=central_heating_forward_temperature,
                return_temperature=central_heating_return_temperature,
                heat_source_cooling=snakemake.params.heat_source_cooling,
            ).assign_coords(heat_source=heat_source_key)
            for heat_source_key in heat_sources
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.heat_source_preheater_utilisation_profiles)
