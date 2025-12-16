# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build heat source utilisation profiles for district heating networks.

This script calculates when and how much heat from various sources (geothermal,
PTES, river water, etc.) can be used, based on the temperature relationship
between the heat source and the district heating network.

Two utilisation modes are calculated:

1. **Direct utilisation**: When the source temperature meets or exceeds the
   forward temperature (T_source ≥ T_forward), the heat source can directly
   supply the district heating network. Profile value is 1.0 (full utilisation)
   or 0.0 (not possible).

2. **Preheater utilisation**: When the source temperature is between the return
   and forward temperatures (T_return < T_source < T_forward), the heat source
   can preheat the return flow before a heat pump lifts it to forward temperature.
   The profile value represents that share of the heat above the return temperature which is utilised to increase the heat pump's sink inflow temperature. The return flow serves as the source inlet.

These profiles are used by ``prepare_sector_network.py`` to configure heat
utilisation links that model cascading temperature use: direct supply when
possible, preheating when beneficial, with heat pumps handling the final lift.

Relevant Settings
-----------------
.. code:: yaml

    sector:
        heat_sources:
            urban central:
                - air
                - geothermal
        district_heating:
            heat_source_cooling: 6  # K
            geothermal:
                constant_temperature_celsius: 65

Inputs
------
- ``resources/<run_name>/central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Forward temperature profiles for district heating networks (°C).
- ``resources/<run_name>/central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Return temperature profiles for district heating networks (°C).
- Heat source temperature profiles (for variable-temperature sources like PTES, air, ground).

Outputs
-------
- ``resources/<run_name>/heat_source_direct_utilisation_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Direct utilisation profiles indexed by (time, name, heat_source).
    Values: 1.0 when T_source ≥ T_forward, 0.0 otherwise.
- ``resources/<run_name>/heat_source_preheater_utilisation_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Preheater utilisation profiles indexed by (time, name, heat_source).
    Values: heat extraction efficiency when T_return < T_source < T_forward, 0.0 otherwise.
"""

import logging

import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config
from scripts.definitions.heat_source import HeatSource

logger = logging.getLogger(__name__)


def get_source_temperature(
    snakemake_params: dict, snakemake_input: dict, heat_source_name: str
) -> float | xr.DataArray:
    """
    Get the temperature profile or constant value for a heat source.

    Parameters
    ----------
    snakemake_params : dict
        Snakemake parameters containing heat_source_temperatures dict for
        constant-temperature sources.
    snakemake_input : dict
        Snakemake input files containing temperature profiles for variable sources.
    heat_source_name : str
        Name of the heat source (e.g., 'geothermal', 'ptes', 'air', 'electrolysis_excess').

    Returns
    -------
    float | xr.DataArray
        Either a constant temperature (float) for sources like geothermal,
        or a DataArray with time-varying temperatures for sources like PTES or air.

    Raises
    ------
    ValueError
        If the required temperature data is not available in params or inputs.
    """
    heat_source = HeatSource(heat_source_name)
    if heat_source.temperature_from_config:
        try:
            return snakemake_params["heat_source_temperatures"][heat_source_name]
        except KeyError:
            raise ValueError(
                f"Constant temperature for heat source '{heat_source_name}' not specified in heat_source_temperatures."
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
    Calculate when a heat source can directly supply district heating.

    Direct utilisation is possible when the source temperature meets or exceeds
    the required forward temperature of the district heating network.

    Parameters
    ----------
    source_temperature : float | xr.DataArray
        Heat source temperature in °C. If float, applies uniformly.
        If DataArray, indexed by (time, name).
    forward_temperature : xr.DataArray
        District heating forward temperature profiles in °C,
        indexed by (time, name).

    Returns
    -------
    xr.DataArray
        Binary profile: 1.0 where T_source ≥ T_forward (direct use possible),
        0.0 otherwise.
    """
    return xr.where(source_temperature >= forward_temperature, 1.0, 0.0)


def get_preheater_utilisation_profile(
    source_temperature: float | xr.DataArray,
    forward_temperature: xr.DataArray,
    return_temperature: xr.DataArray,
    heat_source_cooling: float,
) -> xr.DataArray | float:
    """
    Calculate preheater utilisation efficiency for intermediate-temperature sources.

    When a heat source temperature is between the return and forward temperatures,
    it can preheat the return flow before a heat pump provides the final temperature
    lift. This improves overall efficiency by reducing the heat pump's lift.

    The efficiency represents the fraction of heat extracted from the source that
    goes into preheating (vs. the additional cooling through the heat pump):

        efficiency = (T_source - T_return) / (T_source - T_return + heat_source_cooling)

    Parameters
    ----------
    source_temperature : float | xr.DataArray
        Heat source temperature in °C. If float, applies uniformly.
        If DataArray, indexed by (time, name).
    forward_temperature : xr.DataArray
        District heating forward temperature profiles in °C,
        indexed by (time, name).
    return_temperature : xr.DataArray
        District heating return temperature profiles in °C,
        indexed by (time, name).
    heat_source_cooling : float | xr.DataArray
        Additional temperature drop (K) when extracting heat from the source
        through the heat pump, beyond the preheating contribution.

    Returns
    -------
    xr.DataArray
        Preheater efficiency profile: value in (0, 1) where T_return < T_source < T_forward,
        0.0 otherwise (source too cold or hot enough for direct use).
    """
    return xr.where(
        (source_temperature < forward_temperature)
        * (source_temperature > return_temperature),
        (source_temperature - return_temperature)
        / (source_temperature - return_temperature + heat_source_cooling),
        0.0,
    )


def get_heat_pump_cooling(
    heat_source_name: str,
    default_heat_source_cooling: float,
    snakemake_input: dict,
    return_temperature: xr.DataArray = None,
) -> float | xr.DataArray:
    """
    Get the additional heat source cooling (temperature drop) through the heat pump for a heat source.

    For PTES, this equals the temperature difference between
    return flow and bottom layers (return_temperature - bottom_temperature), which can
    be time-varying. For other sources, uses the default constant value.

    Parameters
    ----------
    heat_source_name : str
        Name of the heat source (e.g., 'ptes', 'geothermal', 'air').
    default_heat_source_cooling : float
        Default heat source cooling in Kelvin, from config.
    snakemake_input : dict
        Snakemake input files, may contain PTES temperature profiles.
    return_temperature : xr.DataArray, optional
        District heating return temperature profiles in °C. Required for PTES.

    Returns
    -------
    float | xr.DataArray
        Heat source cooling in Kelvin. Returns a float for most sources,
        or a DataArray for PTES when temperatures vary with time.

    Raises
    ------
    ValueError
        If heat source is PTES but bottom temperature profile is not provided.
    """
    if heat_source_name == "ptes":
        if "temp_ptes_bottom" not in snakemake_input.keys():
            raise ValueError(
                "PTES heat source requires bottom temperature profile "
                "(temp_ptes_bottom) to calculate heat source cooling."
            )
        if return_temperature is None:
            raise ValueError(
                "PTES heat source requires return_temperature to calculate heat pump cooling."
            )
        ptes_bottom_temperature = xr.open_dataarray(snakemake_input["temp_ptes_bottom"])
        return return_temperature - ptes_bottom_temperature
    return default_heat_source_cooling


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
                heat_source_cooling=get_heat_pump_cooling(
                    heat_source_name=heat_source_key,
                    default_heat_source_cooling=snakemake.params.heat_source_cooling,
                    snakemake_input=snakemake.input,
                    return_temperature=central_heating_return_temperature,
                ),
            ).assign_coords(heat_source=heat_source_key)
            for heat_source_key in heat_sources
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.heat_source_preheater_utilisation_profiles)
