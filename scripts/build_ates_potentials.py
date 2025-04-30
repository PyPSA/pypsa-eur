# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import geopandas as gpd
import pandas as pd
import xarray as xr
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def mwh_ates_per_m2(
    aquifer_volumetric_heat_capacity: float,
    fraction_of_aquifer_area_available: float,
    effective_screen_length: float,
    hot_well_temperature: float,
    cold_well_temperature: float,
    kwh_per_kj: float = 1 / 3600,
    mwh_per_kwh: float = 1 / 1e3,
):
    """
    Calculate the TWh per m2 for ATES systems.
    """

    return (
        aquifer_volumetric_heat_capacity
        * fraction_of_aquifer_area_available
        * effective_screen_length
        * (hot_well_temperature - cold_well_temperature)
        * kwh_per_kj
        * mwh_per_kwh
    )


def suitable_aquifers(
    aquifer_shapes: gpd.GeoDataFrame,
    suitable_aquifer_types: list,
):
    """
    Filter the aquifer shapes by the suitable aquifer types.
    """

    return aquifer_shapes[aquifer_shapes["AQUIF_NAME"].isin(suitable_aquifer_types)]


def ates_potential_per_onshore_region(
    suitable_aquifers: gpd.GeoDataFrame,
    regions_onshore: gpd.GeoDataFrame,
    dh_areas: gpd.GeoDataFrame,
    dh_area_buffer: float,
    mwh_per_m2: float,
):
    """
    Calculate the area of the aquifer shapes in square kilometers.
    """
    ret_val = regions_onshore.copy(deep=True)
    ret_val.index = ret_val["name"]
    ret_val.drop(columns=["name"], inplace=True)

    suitable_aquifers_in_onshore_regions = gpd.overlay(
        suitable_aquifers, regions_onshore, how="intersection"
    )

    dh_areas_buffered = dh_areas.copy(deep=True)
    dh_areas_buffered["geometry"] = dh_areas_buffered.geometry.buffer(dh_area_buffer)

    aquifers_in_dh_areas = (
        gpd.overlay(
            dh_areas_buffered, suitable_aquifers_in_onshore_regions, how="intersection"
        )
        .groupby("name")["geometry"]
        .apply(lambda x: x.area.sum())
    )

    ret_val["ates_potential"] = aquifers_in_dh_areas * mwh_per_m2
    return ret_val


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ates_potential", clusters="48")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # get onshore regions and index them by region name
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore)
    regions_onshore = regions_onshore.to_crs(epsg=3035)

    aquifer_shapes = gpd.read_file(snakemake.input.aquifer_shapes_shp)
    aquifer_shapes = aquifer_shapes.to_crs(regions_onshore.crs)

    dh_areas = gpd.read_file(snakemake.input.dh_areas)
    dh_areas = dh_areas.to_crs(regions_onshore.crs)

    forward_temperature_profiles = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    return_temperature_profiles = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    mwh_per_m2 = mwh_ates_per_m2(
        aquifer_volumetric_heat_capacity=snakemake.params.aquifer_volumetric_heat_capacity,
        fraction_of_aquifer_area_available=snakemake.params.fraction_of_aquifer_area_available,
        effective_screen_length=snakemake.params.effective_screen_length,
        hot_well_temperature=forward_temperature_profiles.mean(dim="time").to_pandas(),
        cold_well_temperature=return_temperature_profiles.to_pandas(),
    )

    suitable_aquifers = suitable_aquifers(
        aquifer_shapes=aquifer_shapes,
        suitable_aquifer_types=snakemake.params.suitable_aquifer_types,
    )

    ates_potentials = ates_potential_per_onshore_region(
        suitable_aquifers=suitable_aquifers,
        regions_onshore=regions_onshore,
        dh_areas=dh_areas,
        dh_area_buffer=snakemake.params.dh_area_buffer,
        mwh_per_m2=mwh_per_m2,
    )

    ates_potentials.to_csv(snakemake.output.ates_potentials, index_label="name")
