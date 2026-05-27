# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build heat source boosting ratio profiles for district heating networks.

This script calculates the boosting ratio ``b`` for each heat source: the
fraction of HP boost needed per unit of source heat to reach the forward
temperature of the district heating network.

**Boosting ratio profile**: For each heat source and timestep, ``b`` is:

- 0 when T_source ≥ T_forward (direct use, no HP boost needed)
- 1 when T_source < T_return (source cannot preheat return flow)
- (T_forward − T_source) / (T_source − T_return) otherwise, clipped to [0, 1]

The boosting ratio is consumed by ``prepare_sector_network.py`` to set the
efficiencies of the heat source utilisation link
(``bus0=resource → bus1=DH heat, bus2=intermediate``):

- **Preheating sources** (PTES, geothermal):
  ``efficiency = 1 + b/cop``, ``efficiency2 = -b``.
  The source directly contributes ``1`` MW of preheat per source MW and the HP
  supplies the ``b`` MW boost via the intermediate bus, consuming
  ``b/cop`` MW electricity. Energy balance per source MW:
  ``in = 1 (source) + b/cop (elec) = out = 1 + b/cop (DH heat)``.
- **Non-preheating limited sources** (river_water):
  ``efficiency = 0``, ``efficiency2 = 1`` — all source heat is routed to the
  HP cold side via the intermediate bus.

Relevant Settings
-----------------
.. code:: yaml

    sector:
        heat_sources:
            urban central:
                - air
                - geothermal
                - ptes
        district_heating:
            geothermal:
                constant_temperature_celsius: 65
            ptes:
                enable: true
                discharge_resistive_boosting: false

Inputs
------
- ``resources/<run_name>/central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Forward temperature profiles for district heating networks (°C).
- ``resources/<run_name>/central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Return temperature profiles for district heating networks (°C).
- Heat source temperature profiles (for variable-temperature sources like PTES, air, ground).

Outputs
-------
- ``resources/<run_name>/heat_source_boosting_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Boosting ratio profiles indexed by (time, name, heat_source).
    Values in [0, 1]: 0 = direct use, 1 = full HP boosting required.
"""

import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config
from scripts.definitions.heat_source import HeatSource


def get_source_temperature(
    snakemake_params: dict, snakemake_input: dict, heat_source_name: str
) -> float | xr.DataArray:
    """
    Get the temperature profile or constant value for a heat source.

    Parameters
    ----------
    snakemake_params : dict
        Snakemake parameters containing constant temperatures for applicable sources.
    snakemake_input : dict
        Snakemake input files containing temperature profiles for variable sources.
    heat_source_name : str
        Name of the heat source (e.g., 'geothermal', 'ptes', 'air').

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


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_heat_source_utilisation_profiles",
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
