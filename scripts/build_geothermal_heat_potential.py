# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build geothermal heat source potentials for district heating networks.

This script maps and aggregates geothermal heat source potentials from LAU-level
data to onshore regions. It converts supply potentials from Manz et al. (2024)
to technical potentials (MW) and scales them based on the actual temperature
differences in the district heating system.

Temperature-Based Scaling
-------------------------
Manz et al. (2024) assume a design temperature difference of 15 K when computing
geothermal heat potentials. This script scales the potentials based on the actual
temperature delta achievable in the district heating system:

a) If source_temperature > forward_temperature (direct utilisation), potentials simply scale linearly in both directions, depending on the actual temperature difference (source_temperature - return_temperature):
   scale_factor = (source_temperature - return_temperature) / 15 K

b) If forward_temperature >= source_temperature > return_temperature (preheating), the actual temperature difference is increased by the additional cooling through the geothermal-sourced heat pump `heat_source_cooling`:
   scale_factor = (source_temperature - return_temperature + heat_source_cooling) / 15 K

c) If source_temperature <= return_temperature (heat pump only), the actual temperature difference is only the exploitation induced by the heat pump:
   scale_factor = heat_source_cooling / 15 K

This results in time-varying heat source power profiles that reflect the actual
extractable heat capacity at each timestep.

Relevant Settings
-----------------
.. code:: yaml

    sector:
        district_heating:
            heat_source_cooling: 6  # K
            heat_source_temperatures:
                geothermal: 65
            dh_areas:
                handle_missing_countries: fill  # or 'ignore' or 'raise'

Inputs
------
- ``data/isi_heat_utilisation_potentials.xlsx``: Heat potentials from Manz et al. (2024)
- ``resources/<run_name>/regions_onshore_base_s_{clusters}.geojson``: Onshore regions
- ``resources/<run_name>/central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``: Forward temperature profiles
- ``resources/<run_name>/central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc``: Return temperature profiles
- ``data/lau_regions.zip``: LAU region geometries

Outputs
-------
- ``resources/<run_name>/heat_source_power_geothermal_base_s_{clusters}_{planning_horizons}.csv``:
  Time-varying geothermal heat source power (MW) with time as index and regions as columns.

Raises
------
ValueError
    If LAU regions in ISI heat potentials are missing from the LAU Regions data.

Source
------
Manz et al. 2024: "Spatial analysis of renewable and excess heat potentials for
climate-neutral district heating in Europe", Renewable Energy, vol. 224, no. 120111,
https://doi.org/10.1016/j.renene.2024.120111
"""

import logging

import geopandas as gpd
import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

ISI_TEMPERATURE_SCENARIOS = {
    65: "low_temperatures",
    85: "high_temperatures",
}
FULL_LOAD_HOURS = 4000
# 3000 for petrothermal

PYPSA_EUR_UNIT = "MWh"

GEOTHERMAL_SOURCE = (
    "Hydrothermal "  # trailing space for hydrothermal necessary to get correct column
)

DESIGN_TEMPERATURE_DIFFERENCE = (
    15  # K - assumed temperature difference for geothermal sources in Manz et al. 2024
)


def scale_heat_source_power(
    heat_source_power: pd.Series,
    forward_temperature: xr.DataArray,
    return_temperature: xr.DataArray,
    source_temperature: float,
    heat_source_cooling: float,
) -> xr.DataArray:
    """
    Scale heat source power based on temperature differences.

    Manz et al. 2024 assume a temperature difference of 15K for geothermal heat
    sources. This function scales the heat source power based on the actual
    temperature difference in the district heating system.

    The scaling logic follows three cases:
    a) If source_temperature > forward_temperature:
       scale_factor = (source_temperature - return_temperature) / 15K
    b) Elif source_temperature > return_temperature:
       scale_factor = (source_temperature - return_temperature + heat_source_cooling) / 15K
    c) Else:
       scale_factor = heat_source_cooling / 15K

    Parameters
    ----------
    heat_source_power : pd.Series
        Base heat source power per region [MW]. Index: region names.
    forward_temperature : xr.DataArray
        Forward temperature profiles [°C]. Dims: (time, name).
    return_temperature : xr.DataArray
        Return temperature profiles [°C]. Dims: (name,).
    source_temperature : float
        Constant geothermal source temperature [°C].
    heat_source_cooling : float
        Temperature drop in heat source when extracting heat via heat pump [K].

    Returns
    -------
    xr.DataArray
        Scaled heat source power [MW]. Dims: (time, name).
    """
    # Ensure alignment of regions
    regions = heat_source_power.index
    forward_temp = forward_temperature.sel(name=regions)
    return_temp = return_temperature.sel(name=regions)

    # Broadcast return_temperature to match forward_temperature dimensions
    return_temp_broadcast = return_temp.broadcast_like(forward_temp)

    # Compute scale factors for each case
    # Case a: source_temp > forward_temp (direct utilisation possible)
    scale_a = (
        source_temperature - return_temp_broadcast
    ) / DESIGN_TEMPERATURE_DIFFERENCE

    # Case b: forward_temp >= source_temp > return_temp (preheating mode)
    scale_b = (
        source_temperature - return_temp_broadcast + heat_source_cooling
    ) / DESIGN_TEMPERATURE_DIFFERENCE

    # Case c: source_temp <= return_temp (HP-only mode)
    scale_c = heat_source_cooling / DESIGN_TEMPERATURE_DIFFERENCE

    # Apply conditional logic
    scale_factor = xr.where(
        source_temperature > forward_temp,
        scale_a,
        xr.where(
            source_temperature > return_temp_broadcast,
            scale_b,
            scale_c,
        ),
    )

    # Convert heat_source_power to DataArray and broadcast
    heat_source_power_da = xr.DataArray(
        heat_source_power.squeeze().values,
        dims=["name"],
        coords={"name": regions},
    )

    # Scale the heat source power
    scaled_power = heat_source_power_da * scale_factor

    # Transpose to (time, name) for consistency with other profiles
    scaled_power = scaled_power.transpose("time", "name")

    return scaled_power


def get_unit_conversion_factor(
    input_unit: str,
    output_unit: str,
    unit_scaling: dict = {"Wh": 1, "kWh": 1e3, "MWh": 1e6, "GWh": 1e9, "TWh": 1e12},
) -> float:
    """
    Get the unit conversion factor between two units.

    Parameters
    ----------
    input_unit : str
        Input unit. Must be one of the keys in `unit_scaling`.
    output_unit : str
        Output unit. Must be one of the keys in `unit_scaling`.
    unit_scaling : dict, optional
        Dictionary of unit scaling factors. Default: {"Wh": 1, "kWh": 1e3, "MWh": 1e6, "GWh": 1e9, "TWh": 1e12}.
    """

    if input_unit not in unit_scaling.keys():
        raise ValueError(
            f"Input unit {input_unit} not allowed. Must be one of {unit_scaling.keys()}"
        )
    elif output_unit not in unit_scaling.keys():
        raise ValueError(
            f"Output unit {output_unit} not allowed. Must be one of {
                unit_scaling.keys()
            }"
        )

    return unit_scaling[input_unit] / unit_scaling[output_unit]


def get_heat_source_power(
    regions_onshore: gpd.GeoDataFrame,
    lau_regions: gpd.GeoDataFrame,
    supply_potentials: pd.Series,
    full_load_hours: float,
    input_unit: str,
    handle_missing_countries: str,
    output_unit: str = "MWh",
) -> pd.Series:
    """
    Map LAU-level supply potentials to onshore regions and convert to power.

    The Manz et al. (2024) supply potentials are provided at LAU level and cover
    EU-27 countries. This function spatially joins them to the onshore regions,
    aggregates, and converts from energy to power using ``full_load_hours``.

    Regions outside the data coverage (e.g. non-EU-27 countries) are handled
    according to ``handle_missing_countries``.

    Parameters
    ----------
    regions_onshore : gpd.GeoDataFrame
        Onshore regions indexed by region name.
    lau_regions : gpd.GeoDataFrame
        LAU region geometries indexed by GISCO_ID.
    supply_potentials : pd.Series
        Supply potentials per LAU region (energy units). Index: GISCO_ID.
    full_load_hours : float
        Full load hours assumed in the supply potential computation.
    input_unit : str
        Unit of the supply potentials.
    handle_missing_countries : str
        Strategy for regions without data:

        - ``"ignore"``: Assume zero potential with a warning.
        - ``"fill"``: Assume zero potential with a warning.
        - ``"raise"``: Raise an error listing the missing regions.
    output_unit : str, optional
        Desired output unit. Default: ``"MWh"``.

    Returns
    -------
    pd.Series
        Heat source power [MW] per onshore region, indexed by region name.

    Raises
    ------
    ValueError
        If ``handle_missing_countries`` is ``"raise"`` and regions are missing.
    """

    unit_conversion_factor = get_unit_conversion_factor(
        input_unit=input_unit,
        output_unit=output_unit,
    )
    scaling_factor = unit_conversion_factor / full_load_hours

    heat_potentials_in_lau = gpd.GeoDataFrame(
        supply_potentials,
        geometry=lau_regions.geometry[supply_potentials.index],
    )

    heat_potentials_in_onshore_regions = gpd.sjoin(
        heat_potentials_in_lau, regions_onshore, how="inner", predicate="intersects"
    )
    heat_potentials_in_onshore_regions_aggregated = (
        heat_potentials_in_onshore_regions.groupby("name").sum(numeric_only=True)
    )

    heat_source_power = heat_potentials_in_onshore_regions_aggregated * scaling_factor

    # Handle regions without data (e.g. non-EU-27 countries)
    non_covered = regions_onshore.index.difference(heat_source_power.index)
    if not non_covered.empty:
        if handle_missing_countries == "raise":
            raise ValueError(
                f"No geothermal heat potentials for {len(non_covered)} regions: "
                f"{non_covered.to_list()}. "
                f"The Manz et al. data covers EU-27 only. "
                f"Set sector.district_heating.dh_areas.handle_missing_countries "
                f"to 'fill' or 'ignore' to fill these with zero."
            )
        elif handle_missing_countries == "fill":
            logger.warning(
                f"No geothermal potentials for {len(non_covered)} regions "
                f"(outside Manz et al. coverage). Filling with zero: "
                f"{non_covered.to_list()}"
            )
        elif handle_missing_countries == "ignore":
            logger.warning(
                f"No geothermal potentials for {len(non_covered)} regions "
                f"(outside Manz et al. coverage). Assuming no potential: "
                f"{non_covered.to_list()}"
            )

    heat_source_power = heat_source_power.reindex(regions_onshore.index, fill_value=0)

    return heat_source_power


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_geothermal_heat_potential",
            clusters=48,
            planning_horizons=2040,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # get onshore regions and index them by region name
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).to_crs("EPSG:4326")
    regions_onshore.index = regions_onshore.name
    regions_onshore.drop(columns=["name"], inplace=True)

    # get LAU regions and index them by LAU-ID
    lau_regions = gpd.read_file(
        f"{snakemake.input.lau_regions}!LAU_RG_01M_2019_3035.geojson",
        crs="EPSG:3035",
    ).to_crs("EPSG:4326")
    lau_regions.index = lau_regions.GISCO_ID

    # temperature scenario that was assumed by Manz et al. when computing potentials is 65C (default) or 85C
    this_temperature_scenario = ISI_TEMPERATURE_SCENARIOS[
        snakemake.params.source_temperature
    ]

    # get heat potentials, index them by LAU-ID and get the geothermal potentials
    isi_heat_potentials = pd.read_excel(
        snakemake.input.isi_heat_potentials,
        sheet_name="Matching_results",
        index_col=0,
        header=[0, 1],
    )
    input_unit = isi_heat_potentials[
        (
            "Unnamed: 2_level_0",
            "Unit",
        )
    ][0]
    geothermal_supply_potentials = isi_heat_potentials[
        (
            f"Supply_potentials_{this_temperature_scenario}",
            GEOTHERMAL_SOURCE,
        )
    ].drop(index="Total")

    # check if all LAU regions in ISI heat potentials are present in LAU Regions data
    if not geothermal_supply_potentials.index.isin(lau_regions.index).all():
        raise ValueError(
            "Some LAU regions in ISI heat potentials are missing from the LAU Regions data."
        )

    # get heat source power by mapping heat potentials to onshore regions and scaling from supply potentials to technical potentials
    heat_source_power = get_heat_source_power(
        regions_onshore=regions_onshore,
        lau_regions=lau_regions,
        supply_potentials=geothermal_supply_potentials,
        full_load_hours=FULL_LOAD_HOURS,
        input_unit=input_unit,
        handle_missing_countries=snakemake.params.handle_missing_countries,
    )

    # Load temperature profiles for scaling
    forward_temperature = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    return_temperature = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    # Scale heat source power based on temperature differences
    scaled_heat_source_power = scale_heat_source_power(
        heat_source_power=heat_source_power,
        forward_temperature=forward_temperature,
        return_temperature=return_temperature,
        source_temperature=snakemake.params.source_temperature,
        heat_source_cooling=snakemake.params.heat_source_cooling,
    )

    # Convert to DataFrame (time x regions) and save as CSV
    scaled_heat_source_power.to_pandas().to_csv(snakemake.output.heat_source_power)
