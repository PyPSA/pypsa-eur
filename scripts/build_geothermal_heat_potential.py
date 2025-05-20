# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build heat source potentials for a given heat source.

This script maps and aggregates geothermal heat source potentials `onshore_regions`. Input data is provided on LAU-level and is aggregated to the onshore regions.
It scales the heat source utilisation potentials to technical potentials by dividing the utilisation potentials by the full load hours of the heat source, also taking into account the energy unit set for the respective source in the config.


Relevant Settings
-----------------
.. code:: yaml
    sector:
        district_heating:
            limited_heat_sources:
                geothermal:
                    constant_temperature_celsius

Inputs
------
- `resources/<run_name>/regions_onshore.geojson`
- `resources/<run_name>/lau_regions.geojson`
- `resources/<run_name>/isi_heat_potentials.xlsx`

Outputs
-------
- `resources/<run_name>/heat_source_technical_potential_{heat_source}_base_s_{clusters}.csv`

Raises
------
- ValueError if some LAU regions in ISI heat potentials are missing from the LAU Regions data.

Source
----------
- Manz et al. 2024: "Spatial analysis of renewable and excess heat potentials for climate-neutral district heating in Europe", Renewable Energy, vol. 224, no. 120111, https://doi.org/10.1016/j.renene.2024.120111
"""

import logging

import geopandas as gpd
import pandas as pd

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


def identify_non_covered_regions(
    regions_onshore: gpd.GeoDataFrame, heat_source_power: pd.DataFrame
) -> pd.Index:
    """
    Identify regions without heat source power data.

    Parameters
    ----------
    regions_onshore : gpd.GeoDataFrame
        GeoDataFrame of the onshore regions, indexed by region name.
    heat_source_power : pd.DataFrame
        Heat source power data, indexed by region name.

    Returns
    -------
    pd.Index
        Index of regions that have no heat source power data.
    """
    return regions_onshore.index.difference(heat_source_power.index)


def get_heat_source_power(
    regions_onshore: gpd.GeoDataFrame,
    supply_potentials: gpd.GeoDataFrame,
    full_load_hours: float,
    input_unit: str,
    output_unit: str = "MWh",
    ignore_missing_regions: bool = False,
) -> pd.DataFrame:
    """
    Get the heat source power from supply potentials.

    Note
    ----
    Broadcasts to repeat constant heat source power across snapshots.

    Parameters
    ----------
    regions_onshore : gpd.GeoDataFrame
        GeoDataFrame of the onshore regions.
    supply_potentials : gpd.GeoDataFrame
        GeoDataFrame of the heat source supply potentials.
    full_load_hours : float
        Full load hours assumed in the supply potential computation. Used to scale the supply potentials to technical potentials.
    input_unit : str
        Unit of the supply potentials. Used to convert to the output unit.
    output_unit : str, optional
        Unit of the technical potentials. Default: "MWh".

    Returns
    -------
    pd.DataFrame
        Heat source power in the onshore regions. Indexed by name (onshore region).
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

    non_covered_regions = identify_non_covered_regions(
        regions_onshore, heat_source_power
    )

    not_eu_27 = [
        "GB",
        "UA",
        "MD",
        "AL",
        "RS",
        "BA",
        "ME",
        "MK",
        "XK",
    ]

    if not non_covered_regions.empty:
        if all(non_covered_regions.str.contains("|".join(not_eu_27))):
            if ignore_missing_regions:
                logger.warning(
                    f"The onshore regions outside EU 27 ({non_covered_regions.to_list()}) have no heat source power. Filling with zeros."
                )
                heat_source_power = heat_source_power.reindex(
                    regions_onshore.index, fill_value=0
                )
            else:
                raise ValueError(
                    f"The onshore regions outside EU 27 {non_covered_regions.to_list()} have no heat source power. Set the ignore_missing_regions parameter in the config to true if you want to include these countries in your analysis despite missing geothermal data."
                )
        else:
            raise ValueError(
                f"The onshore regions {non_covered_regions.to_list()} have no heat source power. The pre-processing of the potential data might be faulty."
            )

    return heat_source_power


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_geothermal_heat_potential",
            clusters=48,
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
        snakemake.params.constant_temperature_celsius
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

    # get heat source power by mapping heat potentials to onshore regions and scaling to from supply potentials to technical potentials
    heat_source_power = get_heat_source_power(
        regions_onshore=regions_onshore,
        supply_potentials=geothermal_supply_potentials,
        full_load_hours=FULL_LOAD_HOURS,
        input_unit=input_unit,
        ignore_missing_regions=snakemake.params.ignore_missing_regions,
    )

    heat_source_power.to_csv(snakemake.output[0])
