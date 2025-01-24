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
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def get_source_temperature(heat_source_key: str):
    """
    Get the constant temperature of a heat source.

    Args:
    ----
    heat_source_key: str
        The key (name) of the heat source.

    Returns:
    -------
    float
        The constant temperature of the heat source in degrees Celsius.

    Raises:
    ------
    ValueError
        If the heat source is unknown (not in `config`).
    """

    if heat_source_key in snakemake.params.heat_utilisation_potentials.keys():
        return snakemake.params.heat_utilisation_potentials[heat_source_key][
            "constant_temperature_celsius"
        ]
    else:
        raise ValueError(
            f"Unknown heat source {heat_source_key}. Must be one of "
            f"{snakemake.params.heat_sources.keys()}."
        )


def get_profile(
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


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            clusters=48,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    direct_utilisation_heat_sources: list[str] = (
        snakemake.params.direct_utilisation_heat_sources
    )

    central_heating_forward_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )

    xr.concat(
        [
            get_profile(
                source_temperature=get_source_temperature(heat_source_key),
                forward_temperature=central_heating_forward_temperature,
            ).assign_coords(heat_source=heat_source_key)
            for heat_source_key in direct_utilisation_heat_sources
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.direct_heat_source_utilisation_profiles)
