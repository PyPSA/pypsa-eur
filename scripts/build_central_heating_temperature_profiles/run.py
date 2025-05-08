# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Approximate district heating forward and return temperature profiles based on
ambient temperature. The method is based on a reference curve from Pieper et
al. 2019, where for ambient temperatures below 0C, the highest possible forward
temperature is assumed and vice versa for temperatures above 10C. Between these
threshold levels, forward temperatures are linearly interpolated.

By default, `max_forward_temperature` from Euroheat DHC Market Outlook 2024 is used; `min_forward_temperature` and `return_temperature` for Germany is used from AGFW-Hauptbericht 2022.
`min_forward_temperature` and `return_temperature` for other countries are extrapolated based on the ratio between `max_forward_temperature` and `min_forward_temperature` and `return_temperature` for those countries not missing (by default only Germany).

Relevant Settings
-----------------

.. code:: yaml
    sector:
        district_heating:
            max_forward_temperature:
            min_forward_temperature:
            return_temperature:
Inputs
------
- `resources/<run_name>/temp_air_total`: Air temperature

Outputs
-------
- `resources/<run_name>/central_heating_temperature_profiles.nc`:

References
----------
- Pieper, et al. (2019): "Assessment of a combination of three heat sources for heat pumps to supply district heating" (https://doi.org/10.1016/j.energy.2019.03.165).
- AGFW (2022): "Hauptbericht 2022" (https://www.agfw.de/zahlen-und-statistiken/agfw-hauptbericht)
"""

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from _helpers import set_scenario_config
from central_heating_temperature_approximator import (
    CentralHeatingTemperatureApproximator,
)


def extrapolate_missing_supply_temperatures_by_country(
    extrapolate_from: dict, extrapolate_to: dict
) -> xr.DataArray:
    """
    Extrapolates missing supply temperatures by country.

    Parameters
    ----------
        extrapolate_from (dict): A dictionary containing supply temperatures to extrapolate from. Should contain all countries.
        extrapolate_to (dict): A dictionary containing supply temperatures to extrapolate to. Where `country` is present, average ratio between `extrapolate_to[country]` and `extrapolate_from[country]` is applied to all countries for which `country` is not present in `extrapolate_from.keys()`  to infer ratio for extrapolation.

    Returns
    -------
        xr.DataArray: A DataArray containing the extrapolated supply temperatures.
    """

    if not all([key in extrapolate_from.keys() for key in extrapolate_to.keys()]):
        raise ValueError(
            "Not all countries in extrapolate_to are present in extrapolate_from."
        )
    # average ratio between extrapolate_from and extrapolate_to for those countries that are in both dictionaries
    extrapolation_ratio = np.mean(
        [extrapolate_to[key] / extrapolate_from[key] for key in extrapolate_to.keys()]
    )

    # apply extrapolation ratio to all keys missing in extrapolate_to
    return {
        key: (
            extrapolate_to[key]
            if key in extrapolate_to.keys()
            else extrapolate_from[key] * extrapolation_ratio
        )
        for key in extrapolate_from.keys()
    }


def get_country_from_node_name(node_name: str) -> str:
    """
    Extracts the country code from a given node name.

    Parameters
    ----------
        node_name (str): The name of the node.

    Returns
    -------
        str: The country code extracted from the node name.
    """
    return node_name[:2]


def map_temperature_dict_to_onshore_regions(
    supply_temperature_by_country: dict,
    regions_onshore: pd.Index,
) -> xr.DataArray:
    """
    Map dictionary of temperatures to onshore regions.

    Missing values are replaced by the mean of all values.

    Parameters
    ----------
    supply_temperature_by_country : dictionary
        Dictionary with temperatures as values and country keys as keys.
    regions_onshore : pd.Index
        Names of onshore regions

    Returns
    -------
    xr.DataArray
        The dictionary values mapped to onshore regions with onshore regions as coordinates.
    """
    return xr.DataArray(
        [
            (
                supply_temperature_by_country[get_country_from_node_name(node_name)]
                if get_country_from_node_name(node_name)
                in supply_temperature_by_country.keys()
                else np.mean(list(supply_temperature_by_country.values()))
            )
            for node_name in regions_onshore.values
        ],
        dims=["name"],
        coords={"name": regions_onshore},
    )


def scale_temperature_to_investment_year(
    temperature_baseyear: dict,
    relative_annual_temperature_reduction: float,
    investment_year: int,
    current_year: int,
) -> dict:
    """
    Scale temperature for each country to investment year by constant
    exponential factor.

    Parameters
    ----------
    temperature_baseyear : dict
        Dictionary with the current temperatures as values and country keys as keys.
    relative_annual_temperature_reduction : float
        The annual reduction rate of the temperature, must be between 0 and 1.
    investment_year : int
        The year in which the investment is planned.
    current_year : int
        The current year.

    Returns
    -------
    dict
        A dictionary with the temperature for each country in the investment year.

    Raises
    ------
    ValueError
        If the investment year is before the current year.
        If the relative annual temperature reduction is not between 0 and 1.
    """

    if investment_year < current_year:
        raise ValueError(
            f"Error: Investment year {investment_year} is before current year {current_year}."
        )
    if (
        relative_annual_temperature_reduction < 0
        or relative_annual_temperature_reduction > 1
    ):
        raise ValueError(
            f"Error: Relative annual temperature reduction must be between 0 and 1, but is {relative_annual_temperature_reduction}."
        )

    return {
        key: temperature_baseyear[key]
        * (
            (1 - relative_annual_temperature_reduction)
            ** (investment_year - current_year)
        )
        for key in temperature_baseyear.keys()
    }


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            clusters=48,
            planning_horizons="2050",
        )

    set_scenario_config(snakemake)

    # reduce temperatures to investment years by constant exponential factor
    max_forward_temperature_investment_year = scale_temperature_to_investment_year(
        temperature_baseyear=snakemake.params.max_forward_temperature_central_heating_baseyear,
        relative_annual_temperature_reduction=snakemake.params.relative_annual_temperature_reduction,
        investment_year=int(snakemake.wildcards.planning_horizons),
        current_year=int(snakemake.params.energy_totals_year),
    )

    min_forward_temperature_investment_year = scale_temperature_to_investment_year(
        temperature_baseyear=snakemake.params.min_forward_temperature_central_heating_baseyear,
        relative_annual_temperature_reduction=snakemake.params.relative_annual_temperature_reduction,
        investment_year=int(snakemake.wildcards.planning_horizons),
        current_year=int(snakemake.params.energy_totals_year),
    )

    return_temperature_investment_year = scale_temperature_to_investment_year(
        temperature_baseyear=snakemake.params.return_temperature_central_heating_baseyear,
        relative_annual_temperature_reduction=snakemake.params.relative_annual_temperature_reduction,
        investment_year=int(snakemake.wildcards.planning_horizons),
        current_year=int(snakemake.params.energy_totals_year),
    )

    # min_forward_temperature and return_temperature contain only values for Germany by default
    # extrapolate missing values based on ratio between max_forward_temperature and min_forward_temperature / return_temperature for Germany (or other countries provided in min_forward_temperature_baseyear and return_temperature_baseyear)
    min_forward_temperature_investment_year = (
        extrapolate_missing_supply_temperatures_by_country(
            extrapolate_from=max_forward_temperature_investment_year,
            extrapolate_to=min_forward_temperature_investment_year,
        )
    )
    return_temperature_investment_year = (
        extrapolate_missing_supply_temperatures_by_country(
            extrapolate_from=max_forward_temperature_investment_year,
            extrapolate_to=return_temperature_investment_year,
        )
    )

    # map forward and return temperatures specified on country-level to onshore regions
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore)["name"]
    snapshots = pd.date_range(freq="h", **snakemake.params.snapshots)
    max_forward_temperature_central_heating_by_node_and_time: xr.DataArray = (
        map_temperature_dict_to_onshore_regions(
            supply_temperature_by_country=max_forward_temperature_investment_year,
            regions_onshore=regions_onshore,
        )
    )
    min_forward_temperature_central_heating_by_node_and_time: xr.DataArray = (
        map_temperature_dict_to_onshore_regions(
            supply_temperature_by_country=min_forward_temperature_investment_year,
            regions_onshore=regions_onshore,
        )
    )
    return_temperature_central_heating_by_node_and_time: xr.DataArray = (
        map_temperature_dict_to_onshore_regions(
            supply_temperature_by_country=return_temperature_investment_year,
            regions_onshore=regions_onshore,
        )
    )

    central_heating_temperature_approximator = CentralHeatingTemperatureApproximator(
        ambient_temperature=xr.open_dataarray(snakemake.input.temp_air_total),
        max_forward_temperature=max_forward_temperature_central_heating_by_node_and_time,
        min_forward_temperature=min_forward_temperature_central_heating_by_node_and_time,
        fixed_return_temperature=return_temperature_central_heating_by_node_and_time,
        lower_threshold_ambient_temperature=snakemake.params.lower_threshold_ambient_temperature,
        upper_threshold_ambient_temperature=snakemake.params.upper_threshold_ambient_temperature,
        rolling_window_ambient_temperature=snakemake.params.rolling_window_ambient_temperature,
    )

    central_heating_temperature_approximator.forward_temperature.to_netcdf(
        snakemake.output.central_heating_forward_temperature_profiles
    )
    central_heating_temperature_approximator.return_temperature.to_netcdf(
        snakemake.output.central_heating_return_temperature_profiles
    )
