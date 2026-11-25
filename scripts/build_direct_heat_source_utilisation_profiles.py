# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build availability profiles for direct heat source utilisation (1 in regions and time steps where heat source can be utilised, 0 otherwise).
When direct utilisation is possible, heat pump COPs are set to zero (c.f. `build_cop_profiles`).

A source between the return and the forward temperature cannot be used
directly. It first preheats the return flow in an ideal heat exchanger, which
cools it to the return temperature, and a heat pump then cools it further by
`heat_source_cooling`. The preheating profile is the share of the extracted
heat that the heat exchanger transfers.
"""

import logging

import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def get_profile(
    source_temperature: float | xr.DataArray, forward_temperature: xr.DataArray
) -> xr.DataArray:
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


def get_preheating_profile(
    source_temperature: float | xr.DataArray,
    forward_temperature: xr.DataArray,
    return_temperature: xr.DataArray,
    heat_source_cooling: float,
) -> xr.DataArray:
    """Get the share of source heat that preheats the return flow before the heat pump."""
    preheating = source_temperature - return_temperature
    return xr.where(
        (source_temperature < forward_temperature) & (preheating > 0),
        preheating / (preheating + heat_source_cooling),
        0.0,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_direct_heat_source_utilisation_profiles",
            horizon=2030,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    source_temperatures: dict[str, float | xr.DataArray] = {
        heat_source: settings["constant_temperature_celsius"]
        for heat_source, settings in snakemake.params.limited_heat_sources.items()
        if settings["constant_temperature_celsius"] is not False
    }
    source_temperatures["ptes"] = xr.open_dataarray(snakemake.input.temp_ptes)

    forward_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    return_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    xr.concat(
        [
            get_profile(
                source_temperature=source_temperature,
                forward_temperature=forward_temperature,
            ).assign_coords(heat_source=heat_source)
            for heat_source, source_temperature in source_temperatures.items()
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.direct_heat_source_utilisation_profiles)

    xr.concat(
        [
            get_preheating_profile(
                source_temperature=source_temperature,
                forward_temperature=forward_temperature,
                return_temperature=return_temperature,
                heat_source_cooling=snakemake.params.heat_source_cooling,
            ).assign_coords(heat_source=heat_source)
            for heat_source, source_temperature in source_temperatures.items()
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.heat_source_preheating_profiles)
