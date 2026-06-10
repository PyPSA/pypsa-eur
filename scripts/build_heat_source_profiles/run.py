# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build heat source preheater profiles and direct-utilisation profiles for urban central heating as well as COP profiles for all heat sources and heat systems.

**Preheater profile**: For each heat source, region and timestep, ``b`` is the
fraction of source heat that pre-heats the return flow:

- 0 when T_source ≥ T_forward (direct use, no HP boost needed)
- 1 when T_source < T_return (source cannot preheat return flow)
- (T_forward - T_source) / (T_source - T_return + heat_pump_cooling) otherwise, clipped to [0, 1]

**Direct-utilisation profile**: For each heat source, region and timestep, 1 if T_source ≥ T_forward (direct use possible), 0 otherwise.

**COP profile**: For each heat source, region and timestep, the COP is the ratio of heat output and electrical input.
COP values below 1 are set to zero (infeasible operating conditions).

For central heating (district heating), the approximation is based on Jensen et al. (2018) using a thermodynamic model (see :class:`CentralHeatingCopApproximator`). For decentral heating (individual heating), the approximation uses quadratic regression from Staffell et al. (2012) (see :class:`DecentralHeatingCopApproximator`).

For pre-heating sources in central heating, the COP is calculated iteratively based on the preheater utilisation and the resulting source cooling, which affects the source temperature and thus the COP. The iteration continues until convergence or a maximum number of iterations is reached.

All profiles are consumed by ``prepare_sector_network.py`` 

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
            heat_source_cooling:
            heat_pump_cooling_iterative:
            log_heat_pump_cooling_iterations:
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
- ``central_heating_forward_temperature_profiles``: DH forward (supply) temperature (°C).
- ``central_heating_return_temperature_profiles``: DH return temperature (°C).
- ``temp_<source>``: temperature profiles for variable-temperature sources (air, ground, PTES).
- ``regions_onshore``: onshore region geometries.

Outputs
-------
- ``cop_profiles``: COP indexed by (time, name, heat_source, heat_system).
- ``heat_source_direct_utilisation_profiles``: 1 where T_source >= T_forward.
- ``heat_source_preheater_utilisation_profiles``: preheater efficiency in [0, 1].
- ``heat_source_cooling_profiles``: source-side cooling applied by the heat pump to preheating sources (K).

With ``log_heat_pump_cooling_iterations`` enabled, the per-iteration cooling
solve is also written to ``heat_pump_cooling_iterations_*.csv`` for inspection.
"""

from pathlib import Path

import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config
from scripts.build_heat_source_profiles.central_heating_cop_approximator import (
    CentralHeatingCopApproximator,
)
from scripts.build_heat_source_profiles.decentral_heating_cop_approximator import (
    DecentralHeatingCopApproximator,
)
from scripts.definitions.heat_source import HeatSource
from scripts.definitions.heat_system_type import HeatSystemType

def get_source_temperature(
    snakemake_params: dict, snakemake_input: dict, heat_source_name: str
) -> float | xr.DataArray:
    """
    Get the temperature profile or constant value for a heat source.

    Sources flagged as constant-temperature (e.g. geothermal) read their value
    from the ``constant_temperature_<source>`` parameter. All other sources are
    time-varying and their profile is loaded from the ``temp_<source>`` input
    file. The distinction is taken from ``HeatSource.has_constant_temperature``.

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


def get_cop(
    heat_system_type: str,
    heat_source: str,
    forward_temperature: xr.DataArray,
    source_temperature: float | xr.DataArray,
    heat_pump_cooling: float | xr.DataArray,
) -> xr.DataArray:
    """
    Coefficient of performance (COP) for a heating system.

    Central heating uses the Jensen et al. (2018) thermodynamic model, decentral
    heating the Staffell et al. (2012) regression. Source and sink inlet
    temperatures follow from the preheating logic; for central systems the source
    outlet is the source inlet cooled by ``heat_pump_cooling``.

    Parameters
    ----------
    heat_system_type : str
        Heating system type, e.g. "urban central", "urban decentral", "rural".
    heat_source : str
        Heat source, e.g. "air", "ground", "ptes", "geothermal".
    forward_temperature : xr.DataArray
        District-heating forward (sink outlet) temperature in °C.
    source_temperature : float | xr.DataArray
        Heat source temperature in °C.
    heat_pump_cooling : float | xr.DataArray
        Source-side cooling depth in K (central systems only).

    Returns
    -------
    xr.DataArray
        COP values, with values below 1 set to zero.
    """
    source_inlet = get_source_inlet_temperature(
        source_temperature, central_heating_return_temperature
    )
    sink_inlet = get_sink_inlet_temperature(
        source_temperature, central_heating_return_temperature
    )

    if HeatSystemType(heat_system_type).is_central:
        cop_params = snakemake.params.heat_pump_cop_approximation_central_heating
        return CentralHeatingCopApproximator(
            sink_outlet_temperature_celsius=forward_temperature,
            sink_inlet_temperature_celsius=sink_inlet,
            source_inlet_temperature_celsius=source_inlet,
            source_outlet_temperature_celsius=source_inlet - heat_pump_cooling,
            refrigerant=cop_params["refrigerant"],
            delta_t_pinch_point=cop_params[
                "heat_exchanger_pinch_point_temperature_difference"
            ],
            isentropic_compressor_efficiency=cop_params[
                "isentropic_compressor_efficiency"
            ],
            heat_loss=cop_params["heat_loss"],
            min_delta_t_lift=cop_params["min_delta_t_lift"],
        ).cop

    return DecentralHeatingCopApproximator(
        sink_outlet_temperature_celsius=snakemake.params.heat_pump_sink_T_decentral_heating,
        source_inlet_temperature_celsius=source_inlet,
        source_type=heat_source,
    ).cop


def get_source_inlet_temperature(
    source_temperature: float | xr.DataArray,
    central_heating_return_temperature: xr.DataArray,
) -> float | xr.DataArray:
    """
    Determine the effective source inlet temperature for the heat pump.

    For pre-heating heat sources (source temp above return, e.g., PTES), the source pre-heats the return flow. In this case, the heat pump sees the return temperature as its effective source inlet (since it lifts from there).
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
    When the source is warmer than the return, the preheater is used and the heat
    pump lifts from the return temperature; otherwise it draws directly from the
    source. We assume ideal heat exchangers with no temperature losses.
    """
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
    heat_pump_cooling: float | xr.DataArray,
) -> xr.DataArray | float:
    """
    Calculate preheater utilisation efficiency for intermediate-temperature sources.

    When a heat source temperature is between the return and forward temperatures,
    it can preheat the return flow before a heat pump provides the final temperature
    lift. This improves overall efficiency by reducing the heat pump's lift.

    The efficiency represents the fraction of heat extracted from the source that
    goes into preheating (vs. the additional cooling through the heat pump):

        efficiency = (T_source - T_return) / (T_source - T_return + heat_pump_cooling)

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
    heat_pump_cooling : float | xr.DataArray
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
        / (source_temperature - return_temperature + heat_pump_cooling),
        0.0,
    )

def compute_heat_pump_cooling(
    heat_system_type: str,
    heat_source: str,
    forward_temperature: xr.DataArray,
    source_temperature: float | xr.DataArray,
    initial_cooling: float = 0.0,
    tolerance: float = 1e-3,
    max_iterations: int = 100,
    iteration_log: list | None = None,
) -> xr.DataArray:
    """
    Source-side cooling applied by the heat pump to a preheating heat source.

    The cooling equals ``(COP - 1) / COP * (T_forward - T_source)``. Since the
    COP itself depends on the cooling, we start from ``initial_cooling`` and
    re-evaluate the COP and cooling in turn until the cooling stops changing
    (within ``tolerance``). Hours that need no boost (T_forward <= T_source) and
    infeasible operating points (COP < 1) get zero cooling.

    If ``iteration_log`` is given, each iteration's full (time, name) state is
    appended to it for inspection; iteration 0 holds the ``initial_cooling``
    starting guess.
    """
    lift = (forward_temperature - source_temperature).clip(min=0)
    cooling = initial_cooling

    for iteration in range(max_iterations):
        cop = get_cop(
            heat_system_type=heat_system_type,
            heat_source=heat_source,
            forward_temperature=forward_temperature,
            source_temperature=source_temperature,
            heat_pump_cooling=cooling,
        )

        if iteration_log is not None:
            # ``* 0 +`` broadcasts scalars onto the (time, name) grid.
            frame = (
                xr.Dataset(
                    {
                        "forward_temperature": forward_temperature,
                        "source_temperature": forward_temperature * 0 + source_temperature,
                        "heat_pump_cooling": forward_temperature * 0 + cooling,
                        "cop": cop,
                    }
                )
                .to_dataframe()
                .reset_index()
            )
            frame.insert(0, "iteration", iteration)
            frame.insert(0, "heat_source", heat_source)
            frame.insert(0, "heat_system", heat_system_type)
            iteration_log.append(frame)

        updated = (1 - 1 / cop.clip(min=1)) * lift
        if float(abs(updated - cooling).max()) < tolerance:
            return updated
        cooling = updated

    raise RuntimeError(
        f"Heat pump cooling did not converge within {max_iterations} "
        f"iterations for '{heat_source}' ({heat_system_type})."
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_heat_source_profiles",
            clusters=48,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Collect the per-iteration cooling solve when logging is enabled.
    cooling_iteration_log = (
        [] if snakemake.params.log_heat_pump_cooling_iterations else None
    )

    central_heating_forward_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    central_heating_return_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    # Cache the cooling assigned to each urban-central source so the COP profile
    # and the preheater utilisation profile stay consistent.
    central_heat_pump_cooling: dict[str, float | xr.DataArray] = {}

    cop_all_system_types = []
    heat_pump_cooling_all_system_types = []
    for heat_system_type, system_heat_sources in snakemake.params.heat_sources.items():
        if not system_heat_sources:
            cop_all_system_types.append(
                central_heating_forward_temperature.expand_dims(
                    heat_source=pd.Index([], name="heat_source")
                ).isel(heat_source=slice(0, 0))
            )
            heat_pump_cooling_all_system_types.append(
                central_heating_forward_temperature.expand_dims(
                    heat_source=pd.Index([], name="heat_source")
                ).isel(heat_source=slice(0, 0))
            )
            continue

        cop_this_system_type = []
        heat_pump_cooling_this_system_type = []
        for heat_source_name in system_heat_sources:
            source_temperature = get_source_temperature(
                snakemake_params=snakemake.params,
                snakemake_input=snakemake.input,
                heat_source_name=heat_source_name,
            )

            # Preheating sources solve the operating-point-consistent cooling
            # (when enabled), seeded from the flat heat_source_cooling value.
            # Otherwise (toggle off, non-preheating, or decentral) the flat
            # heat_source_cooling constant is used, as in the old code.
            if (
                snakemake.params.heat_pump_cooling_iterative
                and HeatSystemType(heat_system_type).is_central
                and HeatSource(heat_source_name).requires_preheater
            ):
                heat_pump_cooling = compute_heat_pump_cooling(
                    heat_system_type=heat_system_type,
                    heat_source=heat_source_name,
                    forward_temperature=central_heating_forward_temperature,
                    source_temperature=source_temperature,
                    initial_cooling=snakemake.params.heat_source_cooling,
                    iteration_log=cooling_iteration_log,
                )
            else:
                heat_pump_cooling = xr.full_like(
                    central_heating_forward_temperature,
                    snakemake.params.heat_source_cooling,
                )

            if heat_system_type == HeatSystemType.URBAN_CENTRAL.value:
                central_heat_pump_cooling[heat_source_name] = heat_pump_cooling

            cop_this_system_type.append(
                get_cop(
                    heat_system_type=heat_system_type,
                    heat_source=heat_source_name,
                    forward_temperature=central_heating_forward_temperature,
                    source_temperature=source_temperature,
                    heat_pump_cooling=heat_pump_cooling,
                )
            )
            heat_pump_cooling_this_system_type.append(heat_pump_cooling)

        cop_all_system_types.append(
            xr.concat(
                cop_this_system_type,
                dim=pd.Index(system_heat_sources, name="heat_source"),
            )
        )
        heat_pump_cooling_all_system_types.append(
            xr.concat(
                heat_pump_cooling_this_system_type,
                dim=pd.Index(system_heat_sources, name="heat_source"),
            )
        )

    cop_profiles = xr.concat(
        cop_all_system_types,
        dim=pd.Index(snakemake.params.heat_sources.keys(), name="heat_system"),
    )
    cop_profiles.to_netcdf(snakemake.output.cop_profiles)

    heat_pump_cooling = xr.concat(
        heat_pump_cooling_all_system_types,
        dim=pd.Index(snakemake.params.heat_sources.keys(), name="heat_system"),
    )
    heat_pump_cooling.to_netcdf(snakemake.output.heat_source_cooling_profiles)

    # Write the full per-node, per-timestep Picard trace only when logging is on.
    if cooling_iteration_log:
        iterations_path = (
            Path(snakemake.output.cop_profiles).parent
            / f"heat_pump_cooling_iterations_base_s_{snakemake.wildcards.clusters}_{snakemake.wildcards.planning_horizons}.csv"
        )
        pd.concat(cooling_iteration_log, ignore_index=True).to_csv(
            iterations_path, index=False
        )

    # Utilisation profiles are built for urban-central heat sources only.
    urban_central_heat_sources = snakemake.params.heat_sources[
        HeatSystemType.URBAN_CENTRAL.value
    ]

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
            for heat_source_key in urban_central_heat_sources
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.heat_source_direct_utilisation_profiles)

    xr.concat(
        [
            get_preheater_utilisation_profile(
                source_temperature=get_source_temperature(
                    snakemake_params=snakemake.params,
                    snakemake_input=snakemake.input,
                    heat_source_name=heat_source_key,
                ),
                forward_temperature=central_heating_forward_temperature,
                return_temperature=central_heating_return_temperature,
                heat_pump_cooling=central_heat_pump_cooling[heat_source_key],
            ).assign_coords(heat_source=heat_source_key)
            for heat_source_key in urban_central_heat_sources
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.heat_source_preheater_utilisation_profiles)
