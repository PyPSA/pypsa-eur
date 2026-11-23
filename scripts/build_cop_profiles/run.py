# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build heat pump coefficient-of-performance (COP) time series per clustered region, heat source and heat system.

For central heating, the COP follows the approximation of Jensen et al.
(2018) with the district heating forward and return temperatures from
[build_central_heating_temperature_profiles][] as sink temperatures. For
decentral heating, a quadratic regression on the temperature difference
between source and sink from Staffell et al. (2012) is used. Source
temperatures come from the temperature profiles of the respective heat
source or, where configured, from a constant source temperature. A source
that is warmer than the return flow first preheats the return flow in an
ideal heat exchanger, so the heat pump draws heat at the return temperature
and lifts the preheated flow from the source temperature to the forward
temperature. The COP is set to zero where the temperature lift is infeasible
or the approximation falls below one.

References
----------
- Jensen et al. (2018), Heat pump COP, part 2: Generalized COP estimation of heat pump processes, 13th IIR Gustav Lorentzen Conference on Natural Refrigerants
- Staffell et al. (2012), [A review of domestic heat pumps](https://doi.org/10.1039/C2EE22653G)
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
    heat_source_cooling_central_heating: float,
    cop_approximation_central_heating: dict,
    sink_T_decentral_heating: float,
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
    heat_source_cooling_central_heating : float
        Temperature drop of the heat source in central heating.
    cop_approximation_central_heating : dict
        Parameters of the central heating COP approximation.
    sink_T_decentral_heating : float
        Sink outlet temperature of decentral heating.
    sink_outlet_temperature_celsius : xr.DataArray, optional
        The outlet temperature of the heat sink in Celsius.
    sink_inlet_temperature_celsius : xr.DataArray, optional
        The inlet temperature of the heat sink in Celsius.

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
            - heat_source_cooling_central_heating,
            refrigerant=cop_approximation_central_heating["refrigerant"],
            delta_t_pinch_point=cop_approximation_central_heating[
                "heat_exchanger_pinch_point_temperature_difference"
            ],
            isentropic_compressor_efficiency=cop_approximation_central_heating[
                "isentropic_compressor_efficiency"
            ],
            heat_loss=cop_approximation_central_heating["heat_loss"],
            min_delta_t_lift=cop_approximation_central_heating["min_delta_t_lift"],
        ).cop

    else:
        return DecentralHeatingCopApproximator(
            sink_outlet_temperature_celsius=sink_T_decentral_heating,
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
        Snakemake parameters containing constant temperature values.
    snakemake_input : dict
        Snakemake input files containing temperature profile paths.
    heat_source_name : str
        Name of the heat source (e.g., "air", "ground", "geothermal", "ptes").

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
    if heat_source.has_constant_temperature:
        return snakemake_params[f"constant_temperature_{heat_source_name}"]
    else:
        return xr.open_dataarray(snakemake_input[f"temp_{heat_source_name}"])


def get_source_inlet_temperature(
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
    source_temperature : float | xr.DataArray
        Temperature of the heat source in Celsius.
    central_heating_return_temperature : xr.DataArray
        District heating return temperature in Celsius.

    Returns
    -------
    float | xr.DataArray
        Effective source inlet temperature for the heat pump in Celsius.

    Notes
    -----
    We assume ideal heat exchangers with no temperature losses.
    """
    # When source temperature > return temperature, preheater is used:
    # heat pump lifts from return temperature (after preheating).
    # When source temperature <= return temperature, no preheating:
    # heat pump draws directly from the source.
    return xr.where(
        source_temperature > central_heating_return_temperature,
        central_heating_return_temperature,
        source_temperature,
    )


def get_sink_inlet_temperature(
    source_temperature: float | xr.DataArray,
    central_heating_return_temperature: xr.DataArray,
) -> float | xr.DataArray:
    """
    Determine the effective sink inlet temperature for the heat pump.

    When the source temperature exceeds the return temperature, the preheater
    raises the return flow to the source temperature and the heat pump lifts
    from there to forward. Otherwise no preheating is used and the heat pump
    receives return-temperature water directly.

    Parameters
    ----------
    source_temperature : float | xr.DataArray
        Temperature of the heat source in Celsius.
    central_heating_return_temperature : xr.DataArray
        District heating return temperature in Celsius.

    Returns
    -------
    float | xr.DataArray
        Effective sink inlet temperature for the heat pump in Celsius.
    """
    return xr.where(
        source_temperature > central_heating_return_temperature,
        source_temperature,
        central_heating_return_temperature,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            horizon=2030,
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
                source_temperature=source_temperature_celsius,
                central_heating_return_temperature=central_heating_return_temperature,
            )

            sink_inlet_temperature_celsius = get_sink_inlet_temperature(
                source_temperature=source_temperature_celsius,
                central_heating_return_temperature=central_heating_return_temperature,
            )

            cop_da = get_cop(
                heat_system_type=heat_system_type,
                heat_source=heat_source_name,
                source_inlet_temperature_celsius=source_inlet_temperature_celsius,
                heat_source_cooling_central_heating=snakemake.params.heat_source_cooling_central_heating,
                cop_approximation_central_heating=snakemake.params.heat_pump_cop_approximation_central_heating,
                sink_T_decentral_heating=snakemake.params.heat_pump_sink_T_decentral_heating,
                sink_outlet_temperature_celsius=central_heating_forward_temperature,
                sink_inlet_temperature_celsius=sink_inlet_temperature_celsius,
            )
            cop_this_system_type.append(cop_da)
        cop_all_system_types.append(
            xr.concat(
                cop_this_system_type,
                # object dtype since xarray cannot join pandas 3 str-dtype indexes
                dim=pd.Index(heat_sources, name="heat_source", dtype=object),
            )
        )

    cop_dataarray = xr.concat(
        cop_all_system_types,
        dim=pd.Index(snakemake.params.heat_sources.keys(), name="heat_system"),
    )

    cop_dataarray.to_netcdf(snakemake.output.cop_profiles)
