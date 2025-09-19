# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Build yearly average ambient air temperature raster data.

This script processes the ambient air temperature from the cutout data and calculates
the yearly average for each grid cell. It only keeps temperatures within onshore regions
and stores them as a unary region with coordinates longitude, latitude, and name.
"""

import logging

import geopandas as gpd
import numpy as np
import shapely
import shapely.vectorized as sv
import xarray as xr
from _helpers import configure_logging, set_scenario_config

LATITUDE = "latitude"
LONGITUDE = "longitude"

logger = logging.getLogger(__name__)


def get_data_in_geometry(
    data: xr.DataArray,
    geometry: shapely.geometry.polygon.Polygon,
) -> xr.DataArray:
    """
    Get the mask for the geometry border.

    Args:
        data (xr.DataArray): The data array.
        geometry: The geometry.

    Returns:
        xr.DataArray: The mask.
    """

    # Extract coordinate values from ds (note: coordinate names match those in ds)
    lon2d, lat2d = np.meshgrid(data[LONGITUDE], data[LATITUDE])
    # Create a boolean mask for grid points within the border buffer.
    mask = sv.contains(geometry, lon2d, lat2d)

    # Convert to an xarray DataArray with matching dims and coords.
    mask_da = xr.DataArray(
        mask,
        dims=[LATITUDE, LONGITUDE],
        coords={
            LATITUDE: data[LATITUDE].values,
            LONGITUDE: data[LONGITUDE].values,
        },
    )

    return data.where(mask_da)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_ambient_air_temperature_yearly_average",
            clusters="48",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    ambient_temperature = xr.open_mfdataset(snakemake.input.cutout).temperature - 273.15

    # Load onshore regions
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore)
    regions_onshore.set_index("name", inplace=True)

    # Calculate yearly average temperature for each grid cell
    # and rename the coordinates to match the expected format
    average_temperature_in_cutout = ambient_temperature.mean(dim="time").rename(
        {"y": LATITUDE, "x": LONGITUDE}
    )

    # Get data in unary region
    average_temperature_in_all_onshore_regions = get_data_in_geometry(
        average_temperature_in_cutout, regions_onshore.geometry.union_all()
    )

    # add onshore_region as additional coordinate
    average_temperature_by_region = xr.concat(
        [
            get_data_in_geometry(
                data=average_temperature_in_all_onshore_regions,
                geometry=regions_onshore.loc[region_name].geometry,
            )
            for region_name in regions_onshore.index
        ],
        dim=regions_onshore.index,
    )

    # Save the result
    average_temperature_by_region.to_netcdf(
        snakemake.output.average_ambient_air_temperature
    )
