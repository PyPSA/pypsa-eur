# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build heat pump coefficient-of-performance (COP) profiles for different heat
sources and heating system types. COP values below 1 are set to zero (infeasible
operating conditions).

For central heating (district heating), the approximation is based on
Jensen et al. (2018) using a thermodynamic model (see
:class:`CentralHeatingCopApproximator`). For decentral heating (individual
heating), the approximation uses quadratic regression from Staffell et al.
(2012) (see :class:`DecentralHeatingCopApproximator`).

Relevant Settings
-----------------

.. code:: yaml

    sector:
        heat_pump_sink_T_individual_heating:
        heat_sources:
            urban central:
            urban decentral:
            rural:
        district_heating:
            forward_temperature:
            return_temperature:
            heat_source_cooling:
            heat_pump_cop_approximation:
                refrigerant:
                heat_exchanger_pinch_point_temperature_difference:
                isentropic_compressor_efficiency:
                heat_loss:
                min_delta_t_lift:
            geothermal:
                constant_temperature_celsius:

Inputs
------
- ``resources/<run_name>/central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``:
  District heating forward (supply) temperature profiles.
- ``resources/<run_name>/central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``:
  District heating return temperature profiles.
- ``resources/<run_name>/temp_air_total_base_s_{clusters}.nc``:
  Ambient air temperature (if air source heat pumps are configured).
- ``resources/<run_name>/temp_soil_total_base_s_{clusters}.nc``:
  Ground/soil temperature (if ground source heat pumps are configured).
- ``resources/<run_name>/temp_ptes_top_profiles_base_s_{clusters}_{planning_horizons}.nc``:
  PTES top temperature profiles (if PTES is configured as heat source).

Outputs
-------
- ``resources/<run_name>/cop_profiles_base_s_{clusters}_{planning_horizons}.nc``:
  Heat pump COP profiles with dimensions (time, name, heat_source, heat_system).
"""

import pandas as pd
import xarray as xr

from scripts._helpers import set_scenario_config
from scripts.build_cop_profiles.central_heating_cop_approximator import (
    CentralHeatingCopApproximator,
)
from scripts.build_cop_profiles.decentral_heating_cop_approximator import (
    DecentralHeatingCopApproximator,
)
from scripts.definitions.heat_source import HeatSource
from scripts.definitions.heat_system_type import HeatSystemType


def get_cop(
    heat_system_type: str,
    heat_source: str,
    source_inlet_temperature_celsius: xr.DataArray,
    sink_outlet_temperature_celsius: xr.DataArray = None,
    sink_inlet_temperature_celsius: xr.DataArray = None,
) -> xr.DataArray:
    """
    Calculate the coefficient of performance (COP) for a heating system.

    Uses different approximation methods depending on the heating system type:
    - Central heating: Thermodynamic model from Jensen et al. (2018)
    - Decentral heating: Quadratic regression from Staffell et al. (2012)

    Parameters
    ----------
    heat_system_type : str
        The type of heating system (e.g., "urban central", "urban decentral", "rural").
    heat_source : str
        The heat source used (e.g., "air", "ground", "ptes", "geothermal").
    source_inlet_temperature_celsius : xr.DataArray
        The inlet temperature of the heat source in Celsius.
    sink_outlet_temperature_celsius : xr.DataArray, optional
        The outlet temperature of the heat sink (forward temperature) in Celsius.
        Required for central heating systems.
    sink_inlet_temperature_celsius : xr.DataArray, optional
        The inlet temperature of the heat sink (return temperature) in Celsius.
        Required for central heating systems.

    Returns
    -------
    xr.DataArray
        The calculated COP values. Values below 1 are set to zero.
    """
    if HeatSystemType(heat_system_type).is_central:
        return CentralHeatingCopApproximator(
            sink_outlet_temperature_celsius=sink_outlet_temperature_celsius,
            sink_inlet_temperature_celsius=sink_inlet_temperature_celsius,
            source_inlet_temperature_celsius=source_inlet_temperature_celsius,
            source_outlet_temperature_celsius=source_inlet_temperature_celsius
            - snakemake.params.heat_source_cooling_central_heating,
            refrigerant=snakemake.params.heat_pump_cop_approximation_central_heating[
                "refrigerant"
            ],
            delta_t_pinch_point=snakemake.params.heat_pump_cop_approximation_central_heating[
                "heat_exchanger_pinch_point_temperature_difference"
            ],
            isentropic_compressor_efficiency=snakemake.params.heat_pump_cop_approximation_central_heating[
                "isentropic_compressor_efficiency"
            ],
            heat_loss=snakemake.params.heat_pump_cop_approximation_central_heating[
                "heat_loss"
            ],
            min_delta_t_lift=snakemake.params.heat_pump_cop_approximation_central_heating[
                "min_delta_t_lift"
            ],
        ).cop

    else:
        return DecentralHeatingCopApproximator(
            sink_outlet_temperature_celsius=snakemake.params.heat_pump_sink_T_decentral_heating,
            source_inlet_temperature_celsius=source_inlet_temperature_celsius,
            source_type=heat_source,
        ).cop


def get_source_temperature(
    snakemake_params: dict, snakemake_input: dict, heat_source_name: str
) -> float | xr.DataArray:
    """
    Retrieve the temperature of a heat source.

    Heat sources can have either constant temperatures (specified in config)
    or time-varying temperatures (loaded from input files).

    Parameters
    ----------
    snakemake_params : dict
        Snakemake parameters containing heat_source_temperatures dict.
    snakemake_input : dict
        Snakemake input files containing temperature profile paths.
    heat_source_name : str
        Name of the heat source (e.g., "air", "ground", "geothermal", "ptes", "electrolysis_waste").

    Returns
    -------
    float | xr.DataArray
        Temperature in Celsius. Returns a float for constant-temperature sources
        or an xr.DataArray for time-varying sources.

    Raises
    ------
    ValueError
        If a constant-temperature source lacks its parameter or a time-varying
        source lacks its input file.
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


def get_source_inlet_temperature(
    heat_source_name: str,
    source_temperature: float | xr.DataArray,
    central_heating_return_temperature: xr.DataArray,
) -> float | xr.DataArray:
    """
    Determine the effective source inlet temperature for the heat pump.

    For heat sources with preheating capability (e.g., PTES), when the source
    temperature exceeds the return temperature, a preheater can be used to
    preheat the return flow. In this case, the heat pump sees the return
    temperature as its effective source inlet (since it lifts from there).
    Otherwise, the heat pump draws directly from the source temperature.

    Parameters
    ----------
    heat_source_name : str
        Name of the heat source.
    source_temperature : float | xr.DataArray
        Temperature of the heat source in Celsius.
    central_heating_return_temperature : xr.DataArray
        District heating return temperature in Celsius.

    Returns
    -------
    float | xr.DataArray
        Effective source inlet temperature for the heat pump in Celsius.
    """
    heat_source = HeatSource(heat_source_name)
    if heat_source.requires_preheater:
        # When source temperature > return temperature, preheater is used:
        # heat pump lifts from return temperature (after preheating).
        # When source temperature <= return temperature, no preheating:
        # heat pump draws directly from the source.
        return central_heating_return_temperature.where(
            central_heating_return_temperature < source_temperature, source_temperature
        )
    else:
        return source_temperature


def get_sink_inlet_temperature(
    heat_source_name: str,
    source_temperature: float | xr.DataArray,
    central_heating_return_temperature: xr.DataArray,
    central_heating_forward_temperature: xr.DataArray,
) -> float | xr.DataArray:
    """
    Determine the effective sink inlet temperature for the heat pump.

    For heat sources with preheating capability (e.g., PTES), when the source
    temperature exceeds the return temperature, a preheater raises the return
    flow to forward temperature. The heat pump then lifts from return to forward
    temperature. When preheating is not used (source <= return), the heat pump
    receives water at return temperature and heats it to forward temperature.

    Parameters
    ----------
    heat_source_name : str
        Name of the heat source.
    source_temperature : float | xr.DataArray
        Temperature of the heat source in Celsius.
    central_heating_return_temperature : xr.DataArray
        District heating return temperature in Celsius.
    central_heating_forward_temperature : xr.DataArray
        District heating forward (supply) temperature in Celsius.

    Returns
    -------
    float | xr.DataArray
        Effective sink inlet temperature for the heat pump in Celsius.
    """
    heat_source = HeatSource(heat_source_name)
    if heat_source.requires_preheater:
        # When source temperature > return temperature, preheater is used:
        # preheater raises return flow, heat pump inlet is at forward temp.
        # When source temperature <= return temperature, no preheating:
        # heat pump inlet is at return temperature.
        return central_heating_forward_temperature.where(
            central_heating_return_temperature < source_temperature,
            central_heating_return_temperature,
        )
    else:
        return central_heating_return_temperature


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            clusters=48,
        )

    set_scenario_config(snakemake)

    central_heating_forward_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    central_heating_return_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    cop_all_system_types = []
    for heat_system_type, heat_sources in snakemake.params.heat_sources.items():
        cop_this_system_type = []
        if not heat_sources:
            cop_all_system_types.append(
                central_heating_forward_temperature.expand_dims(
                    heat_source=pd.Index([], name="heat_source")
                ).isel(heat_source=slice(0, 0))
            )
            continue
        for heat_source_name in heat_sources:
            source_temperature_celsius = get_source_temperature(
                snakemake_params=snakemake.params,
                snakemake_input=snakemake.input,
                heat_source_name=heat_source_name,
            )

            source_inlet_temperature_celsius = get_source_inlet_temperature(
                heat_source_name=heat_source_name,
                source_temperature=source_temperature_celsius,
                central_heating_return_temperature=central_heating_return_temperature,
            )

            sink_inlet_temperature_celsius = get_sink_inlet_temperature(
                heat_source_name=heat_source_name,
                source_temperature=source_temperature_celsius,
                central_heating_forward_temperature=central_heating_forward_temperature,
                central_heating_return_temperature=central_heating_return_temperature,
            )

            cop_da = get_cop(
                heat_system_type=heat_system_type,
                heat_source=heat_source_name,
                source_inlet_temperature_celsius=source_inlet_temperature_celsius,
                sink_outlet_temperature_celsius=central_heating_forward_temperature,
                sink_inlet_temperature_celsius=central_heating_return_temperature,
            )
            cop_this_system_type.append(cop_da)
        cop_all_system_types.append(
            xr.concat(
                cop_this_system_type, dim=pd.Index(heat_sources, name="heat_source")
            )
        )

    cop_dataarray = xr.concat(
        cop_all_system_types,
        dim=pd.Index(snakemake.params.heat_sources.keys(), name="heat_system"),
    )

    cop_dataarray.to_netcdf(snakemake.output.cop_profiles)
