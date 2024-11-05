# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
"""

import geopandas as gpd
from _helpers import set_scenario_config
from OnshoreRegionData import OnshoreRegionData


def get_unit_conversion_factor(
    input_unit: str,
    output_unit: str,
    unit_scaling: dict = {"Wh": 1, "kWh": 1e3, "MWh": 1e6, "GWh": 1e9, "TWh": 1e12},
) -> float:

    if input_unit not in unit_scaling.keys():
        raise ValueError(
            f"Input unit {input_unit} not allowed. Must be one of {
                unit_scaling.keys()}"
        )
    elif output_unit not in unit_scaling.keys():
        raise ValueError(
            f"Output unit {output_unit} not allowed. Must be one of {
                unit_scaling.keys()}"
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

    heat_source_technical_potential = {}
    for (
        heat_source,
        heat_source_features,
    ) in snakemake.params.fraunhofer_heat_sources.items():

        heat_source_utilisation_potential = gpd.read_file(snakemake.input[heat_source])

        heat_source_technical_potential[heat_source] = OnshoreRegionData(
            onshore_regions=regions_onshore,
            data=heat_source_utilisation_potential,
            column_name=heat_source_features["column_name"],
            scaling_factor=get_unit_conversion_factor(
                input_unit=heat_source_features["unit"], output_unit="MWh"
            )
            / heat_source_features["full_load_hours"],
        ).data_in_regions_scaled

        heat_source_technical_potential[heat_source].to_csv(
            snakemake.output[heat_source]
        )
