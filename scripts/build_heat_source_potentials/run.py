# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build heat source potentials for a given heat source.

This script maps and aggregates heat source potentials per heat source to `onshore_regions` using `OnshoreRegionData`.
It scales the heat source utilisation potentials to technical potentials by dividing the utilisation potentials by the full load hours of the heat source, also taking into account the energy unit set for the respective source in the config.


Relevant Settings
-----------------
.. code:: yaml
    sector:
        district_heating:
            heat_utilisation_potentials:
    {heat_source}


Inputs
------
- `resources/<run_name>/regions_onshore.geojson`
- `resources/<run_name>/heat_source_utilisation_potentials/<heat_source>.gpkg`

Outputs
-------
- `resources/<run_name>/heat_source_technical_potential_{heat_source}_base_s_{clusters}.csv`
"""

import geopandas as gpd
from _helpers import set_scenario_config

from scripts.build_heat_source_potentials.onshore_region_data import OnshoreRegionData


def get_unit_conversion_factor(
    input_unit: str,
    output_unit: str,
    unit_scaling: dict = {"Wh": 1, "kWh": 1e3, "MWh": 1e6, "GWh": 1e9, "TWh": 1e12},
) -> float:
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


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_heat_source_potentials",
            clusters=48,
        )

    set_scenario_config(snakemake)

    regions_onshore = gpd.read_file(snakemake.input.regions_onshore)
    heat_source_utilisation_potential = gpd.read_file(
        snakemake.input.utilisation_potential
    )

    unit_conversion_factor = get_unit_conversion_factor(
        input_unit=snakemake.params.heat_utilisation_potentials[
            snakemake.wildcards.heat_source
        ]["unit"],
        output_unit="MWh",
    )
    scaling_factor = (
        unit_conversion_factor
        / snakemake.params.heat_utilisation_potentials[snakemake.wildcards.heat_source][
            "full_load_hours"
        ]
    )

    heat_source_technical_potential = OnshoreRegionData(
        onshore_regions=regions_onshore,
        data=heat_source_utilisation_potential,
        column_name=snakemake.params.heat_utilisation_potentials[
            snakemake.wildcards.heat_source
        ]["column_name"],
        scaling_factor=scaling_factor,
    ).data_in_regions_scaled

    heat_source_technical_potential.to_csv(snakemake.output[0])
