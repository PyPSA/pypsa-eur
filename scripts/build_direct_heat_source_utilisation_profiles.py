# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""

"""

import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from typing import List


from _helpers import set_scenario_config


def get_source_temperature(heat_source_key: str):

    if (
        heat_source_key
        in snakemake.params.fraunhofer_heat_utilisation_potentials.keys()
    ):
        return snakemake.params.fraunhofer_heat_utilisation_potentials[heat_source_key][
            "constant_temperature_celsius"
        ]
    else:
        raise ValueError(
            f"Unknown heat source {heat_source_key}. Must be one of {
                snakemake.params.heat_sources.keys()}."
        )


def get_profile(
    source_temperature: float | xr.DataArray, forward_temperature: xr.DataArray
):

    return xr.where(source_temperature >= forward_temperature, 1.0, 0.0)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            clusters=48,
        )

    set_scenario_config(snakemake)

    direct_utilisation_heat_sources: List[str] = (
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
            ).assign_coords(
                heat_source=heat_source_key
            )
            for heat_source_key in direct_utilisation_heat_sources
        ],
        dim="heat_source",
    ).to_netcdf(snakemake.output.direct_heat_source_utilisation_profiles)
