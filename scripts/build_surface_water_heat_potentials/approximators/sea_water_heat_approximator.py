# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
import logging

import geopandas as gpd
import shapely
import xarray as xr

from scripts.build_surface_water_heat_potentials.approximators.surface_water_heat_approximator import (
    SurfaceWaterHeatApproximator,
)

logger = logging.getLogger(__name__)

# Constants
INF = 1e9  # Marker for unusable temperature values


class SeaWaterHeatApproximator(SurfaceWaterHeatApproximator):
    """
    Approximator for sea water heat potential calculations.

    This class extends SurfaceWaterHeatApproximator to handle sea water temperature
    data for heat pump applications. It processes water temperature data to determine
    usable heat potential from sea water sources.
    """

    def __init__(
        self,
        water_temperature: xr.DataArray,
        region: shapely.geometry.polygon.Polygon | gpd.GeoSeries,
        min_inlet_temperature: float = 1,
    ) -> None:
        # buffer the region geometry by half the data resolution
        # This way, offshore data points just outside the region are included
        self.water_temperature = water_temperature
        self.region_geometry = region.geometry.boundary.buffer(
            self._data_resolution / 1.5
        )
        self.min_outlet_temperature = min_inlet_temperature

        # Validate inputs and potentially reproject data
        self._validate_and_reproject_input()

        # Create masked data for processing
        self._clip_data_to_region()

    def _validate_and_reproject_input(self) -> None:
        """
        Validate input data and ensure proper CRS alignment.

        Updates self.water_temperature with properly projected data if needed.

        Raises
        ------
        ValueError
            If inputs are invalid or incompatible
        """
        # Check if data has rio attribute and CRS information
        # Ensure data has rioxarray capabilities
        if not hasattr(self.water_temperature, "rio"):
            raise ValueError("water temperature must have rioxarray capabilities")

        # Ensure data has CRS information
        if not self.water_temperature.rio.crs:
            raise ValueError(
                "water temperature must have CRS information (use rio.write_crs)"
            )

        # Project data to target CRS if needed
        if self.water_temperature.rio.crs.to_epsg() != self.EPSG:
            try:
                self.water_temperature = self.water_temperature.rio.reproject(
                    f"EPSG:{self.EPSG}"
                )
                logger.info(f"Reprojected water_temperature to EPSG:{self.EPSG}")
            except Exception as e:
                raise ValueError(f"Failed to reproject water_temperature: {str(e)}")

    def _clip_data_to_region(self) -> None:
        """
        Mask water temperature to the geometry.
        """

        self._water_temperature_in_region = self.water_temperature.rio.clip(
            self.region_geometry, drop=False
        )

    def get_spatial_aggregate(self) -> xr.Dataset:
        """
        Get the spatial aggregate of water temperature.

        Returns
        -------
        xr.Dataset
            Dataset containing average_temperature
        """
        average_water_temperature = self._water_temperature_in_region.mean(
            dim=["x", "y"], skipna=True
        )

        # Combine into a single dataset and apply cut-off temperature
        return xr.Dataset(
            data_vars={
                "average_temperature": self._get_usable_water_temperature(
                    water_temperature=average_water_temperature
                ),
            }
        )

    def get_temporal_aggregate(self) -> xr.Dataset:
        """
        Get the temporal aggregate of water temperature.

        Returns
        -------
        xr.Dataset
            Dataset containing average_temperature
        """
        average_water_temperature = self._water_temperature_in_region.mean(
            dim=[self.TIME], skipna=True
        )

        # Combine into a single dataset
        # Don't apply cut-off temperature here because this is only used for plotting
        # and analysis
        return xr.Dataset(data_vars={"average_temperature": average_water_temperature})

    def _get_usable_water_temperature(
        self, water_temperature: xr.DataArray
    ) -> xr.DataArray:
        """
        Get the usable water temperature.

        Parameters
        ----------
        water_temperature : xr.DataArray
            Water temperature data

        Returns
        -------
        xr.DataArray
            Usable water temperature with minimum threshold applied
        """
        return xr.where(
            water_temperature > self.min_outlet_temperature,
            water_temperature,
            -INF,
        )
