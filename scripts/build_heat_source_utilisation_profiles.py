# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build heat source boosting ratio profiles for district heating networks.

This script calculates the boosting ratio (b) for each heat source: how
much heat pump output is needed per unit of source heat to reach the forward
temperature of the district heating network.

**Boosting ratio profile**: For each heat source and timestep, b is:

- 0 when T_source ≥ T_forward (direct use, no HP boost needed)
- 1 when T_source < T_return (source cannot preheat return flow)
- (T_forward − T_source) / (T_source − T_return) otherwise, clipped to [0, 1]

The boosting ratio is consumed by ``prepare_sector_network.py`` to set the
efficiencies of the heat source utilisation link (forward, p ≥ 0):

- bus0 (source) → bus1 (DH heat) at efficiency (1 − b)
- bus0 (source) → bus2 (HP input bus) at efficiency2 = b

Energy is conserved: (1 − b) + b = 1.

Relevant Settings
-----------------
.. code:: yaml

    sector:
        heat_sources:
            urban central:
                - air
                - geothermal
        district_heating:
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
    Values in [0, 1]: 0 = direct use, 1 = full HP boosting required.
"""

import logging

import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config
from scripts.definitions.heat_source import HeatSource, HeatSourceType

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
        Name of the heat source (e.g., 'geothermal', 'ptes', 'air', 'electrolysis_waste').

    Returns
    -------
    float | xr.DataArray
        Either a constant temperature (float) for sources like geothermal,
        or a DataArray with time-varying temperatures for sources like PTES or air.

    Notes
    -----
    Presence of constant-temperature entries in ``heat_source_temperatures``
    is validated at config load time by ``SectorConfig``.
    """
    heat_source = HeatSource(heat_source_name)
    if heat_source.temperature_from_config:
        return snakemake_params["heat_source_temperatures"][heat_source_name]
    elif heat_source.source_type == HeatSourceType.STORAGE:
        # PTES layer temperatures are constants from the ptes_operations dataset
        if heat_source_name.startswith("ptes layer"):
            layer_idx = int(heat_source_name.split()[-1])
            return float(ptes_ds["layer_temperatures"].sel(layer=layer_idx).item())
        else:
            return float(ptes_ds.attrs["top_temperature"])
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

    The boosting ratio b represents the fraction of source heat that must be
    routed through the heat pump (rather than used directly) to reach the
    district heating forward temperature:

        b = 0                                                if T_source ≥ T_forward
        b = 1                                                if T_source < T_return
        b = (T_forward − T_source) / (T_source − T_return)  otherwise

    The result is clipped to [0, 1].

    Parameters
    ----------
    source_temperature : float | xr.DataArray
        Heat source temperature in °C.
    forward_temperature : xr.DataArray
        District heating forward temperature in °C, indexed by (time, name).
    return_temperature : xr.DataArray
        District heating return temperature in °C, indexed by (time, name).

    Returns
    -------
    xr.DataArray
        Boosting ratio profile in [0, 1]: 0 = direct use, 1 = full HP boosting.
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

    # PTES uses dedicated routing in prepare_sector_network and no longer
    # consumes the generic boosting profile pathway.
    filtered_heat_sources = [
        hs
        for hs in heat_sources
        if not (ptes_enable and HeatSource(hs).source_type == HeatSourceType.STORAGE)
    ]

    # Load PTES operations dataset if enabled
    if ptes_enable:
        ptes_ds = xr.open_dataset(snakemake.input.ptes_operations)
        num_ptes_layers = int(ptes_ds.attrs["num_layers"])
    else:
        ptes_ds = None
        num_ptes_layers = 0

    central_heating_forward_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    central_heating_return_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    if filtered_heat_sources:
        boosting_profiles = xr.concat(
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
                for heat_source_key in filtered_heat_sources
            ],
            dim="heat_source",
        )
    else:
        boosting_profiles = xr.DataArray(
            [],
            coords={"heat_source": []},
            dims=["heat_source"],
        )

    boosting_profiles.to_netcdf(snakemake.output.heat_source_boosting_profiles)
