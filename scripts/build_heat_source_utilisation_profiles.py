# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build heat source alpha profiles for district heating networks.

This script calculates the boosting ratio (alpha) for each heat source: how
much heat pump output is needed per unit of source heat to reach the forward
temperature of the district heating network.

**Alpha profile**: For each heat source and timestep, alpha is:

- 0 when T_source ≥ T_forward (direct use possible, no HP boost needed)
- (T_forward − T_source) / delta_T otherwise, where delta_T is:
  - heat_pump_cooling                       if T_source < T_return
  - T_source − T_return + heat_pump_cooling if T_return ≤ T_source < T_forward

heat_pump_cooling is the additional temperature drop achievable by the heat
pump beyond the return temperature (i.e. how far below T_return the HP can
extract from the source). For PTES it equals T_return − T_bottom; for other
sources it is the constant config value ``heat_source_cooling``.

The alpha profile is consumed by ``prepare_sector_network.py`` to set the
efficiencies of the heat source utilisation links:

- bus0 (source) → bus1 (DH heat) at efficiency 1 + alpha
- bus2 (HP output bus) drawn at efficiency −alpha per unit of source

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
- ``resources/<run_name>/heat_source_boosting_profiles_base_s_{clusters}_{planning_horizons}.nc``
    Boosting ratio profiles indexed by (time, name, heat_source).
    Values: 0 when T_source ≥ T_forward; (T_fwd − T_src) / delta_T otherwise.
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


def get_boosting_profile(
    source_temperature: float | xr.DataArray,
    forward_temperature: xr.DataArray,
    return_temperature: xr.DataArray,
) -> xr.DataArray:
    """
    Calculate the boosting ratio: HP heat needed per unit of source heat.

    Alpha represents the ratio of heat pump output required to source heat input
    in order to reach the district heating forward temperature:

        alpha = 0                                      if T_source ≥ T_forward
        alpha = (T_fwd − T_src) / delta_T             otherwise

    where delta_T is:

        delta_T = heat_pump_cooling                       if T_source < T_return
        delta_T = T_source − T_return + heat_pump_cooling if T_return ≤ T_source < T_forward

    For PTES (where heat_pump_cooling = T_return − T_bottom):
        delta_T = T_source − T_bottom (when T_source ≥ T_return, which is the normal case)

    Parameters
    ----------
    source_temperature : float | xr.DataArray
        Heat source temperature in °C.
    forward_temperature : xr.DataArray
        District heating forward temperature in °C, indexed by (time, name).
    return_temperature : xr.DataArray
        District heating return temperature in °C, indexed by (time, name).
    heat_pump_cooling : float | xr.DataArray
        Additional temperature drop achievable by HP beyond return temperature (K).

    Returns
    -------
    xr.DataArray
        Alpha profile: 0 where direct use is possible, positive ratio otherwise.
        Shape matches forward_temperature.
    """
    return xr.where(
        source_temperature >= forward_temperature,
        # no boosting needed if source_temp > forward_temp
        0.0,
        xr.where(
            source_temperature < return_temperature,
            # source does not pre-heat return flow if below return temp
            1,
            # if source is between return and forward temp, it pre-heats return flow (part of heat is utilised directly, part is boosted by HP)
            (forward_temperature - source_temperature)
            / (source_temperature - return_temperature),
        ).clip(min=0, max=1),
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
    ptes_enable: bool = snakemake.params.ptes_enable

    # Validate PTES configuration
    if ptes_enable and "ptes" not in heat_sources:
        raise ValueError(
            "PTES is enabled (district_heating.ptes.enable=true) but 'ptes' "
            "is not in heat_sources.urban_central. PTES requires being listed in heat_sources to create the necessary buses and links for heat discharge to the 'urban central heat' bus."
        )

    central_heating_forward_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    central_heating_return_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    xr.concat(
        [
            get_boosting_profile(
                source_temperature=get_source_temperature(
                    snakemake_params=snakemake.params,
                    snakemake_input=snakemake.input,
                    heat_source_name=heat_source_key,
                ),
                forward_temperature=central_heating_forward_temperature,
                return_temperature=central_heating_return_temperature,
            ).assign_coords(heat_source=heat_source_key)
            for heat_source_key in heat_sources
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.heat_source_boosting_profiles)
