# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
import logging
from abc import ABC
from functools import cached_property
from typing import Union

import geopandas as gpd
import numpy as np
import shapely
import xarray as xr

logger = logging.getLogger(__name__)


class SurfaceWaterHeatApproximator(ABC):
    """
    A class for calculating heat source potential for seawater-sourced heat pumps.

    This class encapsulates the full workflow for loading oceanographic data,
    calculating heat potential, resampling, masking to regions of interest,
    and visualizing results.

    Attributes:
        results (xr.Dataset): The results of the calculation. Contains `total_power` and `average_temperature`. Coordinates are `time`.
    """

    TIME = "time"
    EPSG = 3035

    def __init__(
        self,
        volume_flow: xr.DataArray,
        water_temperature: xr.DataArray,
        region: Union[shapely.geometry.polygon.Polygon, gpd.GeoSeries],
        max_relative_volume_flow: float = 1.0,
        delta_t_max: float = 4,
        min_outlet_temperature: float = 1,
        min_distance_meters: int = 2000,
    ) -> None:
        """
        Initialize the SurfaceWaterHeatApproximator. This is an abstract class and should not be instantiated directly.

        Parameters
        ----------
        volume_flow : xr.DataArray
            Volume flow data
        water_temperature : xr.DataArray
            Water temperature data
        region : Union[shapely.geometry.polygon.Polygon, gpd.GeoSeries]
            Region of interest geometry
        max_relative_volume_flow : float, optional
            Maximum relative volume flow, by default 1.0
        delta_t_max : float, optional
            Maximum temperature difference, by default 4
        min_outlet_temperature : float, optional
            Minimum outlet temperature, by default 1
        min_distance_meters : int, optional
            Minimum distance between projects in meters, by default 2000
        """
        # Set instance variables
        self.volume_flow = volume_flow
        self.water_temperature = water_temperature
        self.region = region
        self.max_relative_volume_flow = max_relative_volume_flow
        self.delta_t_max = delta_t_max
        self.min_outlet_temperature = min_outlet_temperature
        self.min_distance_meters = min_distance_meters

        # Validate inputs and potentially reproject data
        # self._validate_and_reproject_input()

        # All expensive computations are now lazy via cached_property

    def get_spatial_aggregate(self) -> xr.Dataset:
        """
        Get the spatial aggregate of water temperature and power.

        Returns
        -------
        xr.Dataset
            Dataset containing total_power and average_temperature
        """
        total_power = self._power_sum_spatial * self._scaling_factor

        # Calculate power-weighted average temperature using cached sum
        average_water_temperature = (
            self._water_temperature_in_region * self._power_in_region
        ).sum(dim=["x", "y"]) / (self._power_sum_spatial + 0.001)

        # Combine into a single dataset
        return xr.Dataset(
            data_vars={
                "total_power": total_power,
                "average_temperature": average_water_temperature,
            }
        )

    def get_temporal_aggregate(self) -> xr.Dataset:
        """
        Get the temporal aggregate of water temperature and power.

        Returns
        -------
        xr.Dataset
            Dataset containing total_energy and average_temperature
        """
        total_energy = self._power_sum_temporal * self._scaling_factor

        # Calculate power-weighted average temperature using cached sum
        average_water_temperature = (
            self._water_temperature_in_region * self._power_in_region
        ).sum(dim=[self.TIME]) / (self._power_sum_temporal + 0.001)

        # Combine into a single dataset
        return xr.Dataset(
            data_vars={
                "total_energy": total_energy,
                "average_temperature": average_water_temperature,
            }
        )

    def _validate_and_reproject_input(self) -> None:
        """
        Validate input data and ensure proper CRS alignment.

        Updates self.volume_flow and self.water_temperature with properly
        projected data if needed.

        Raises:
            ValueError: If inputs are invalid or incompatible
        """
        # Check if data has rio attribute and CRS information
        for name, data in [
            ("water_temperature", self.water_temperature),
            ("volume_flow", self.volume_flow),
        ]:
            # Ensure data has rioxarray capabilities
            if not hasattr(data, "rio"):
                raise ValueError(f"{name} must have rioxarray capabilities")

            # Ensure data has CRS information
            if not data.rio.crs:
                raise ValueError(
                    f"{name} must have CRS information (use rio.write_crs)"
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

        if self.volume_flow.rio.crs.to_epsg() != self.EPSG:
            try:
                self.volume_flow = self.volume_flow.rio.reproject(f"EPSG:{self.EPSG}")
                logger.info(f"Reprojected volume_flow to EPSG:{self.EPSG}")
            except Exception as e:
                raise ValueError(f"Failed to reproject volume_flow: {str(e)}")

        # Check that datasets have the same dimensions
        if not set(self.water_temperature.dims) == set(self.volume_flow.dims):
            raise ValueError(
                f"Dimensions mismatch: water_temperature has {self.water_temperature.dims}, "
                f"volume_flow has {self.volume_flow.dims}"
            )

        # Check that x and y coordinates match
        for coord in ["x", "y", self.TIME]:
            if (
                coord in self.water_temperature.coords
                and coord in self.volume_flow.coords
            ):
                if not np.array_equal(
                    self.water_temperature[coord].values, self.volume_flow[coord].values
                ):
                    raise ValueError(
                        f"{coord} coordinates don't match between datasets"
                    )
            else:
                raise ValueError(
                    f"{coord} coordinate '{coord}' not found in both datasets"
                )

        # For region geometry, we just check the type
        # if not isinstance(self.region, shapely.geometry.multipolygon.MultiPolygon):
        #     raise ValueError(f"region_geometry must be a shapely MultiPolygon, got {type(self.region)}")

    @cached_property
    def _volume_flow_in_region(self) -> xr.DataArray:
        """
        Cache clipped volume flow data.

        Returns
        -------
        xr.DataArray
            Volume flow data clipped to region
        """
        return self.volume_flow.rio.clip(self.region.geometry, drop=False)

    @cached_property
    def _water_temperature_in_region(self) -> xr.DataArray:
        """
        Cache clipped water temperature data.

        Returns
        -------
        xr.DataArray
            Water temperature data clipped to region
        """
        return self.water_temperature.rio.clip(self.region.geometry, drop=False)

    @cached_property
    def _data_resolution(self) -> float:
        """
        Cache resolution calculation based on dataset resolution.
        Assumes data is in EPSG:3035 (meters).
        """
        # Get resolution directly from rio
        x_res, y_res = self.water_temperature.rio.resolution()

        # Average resolution in meters (EPSG:3035 uses meters)
        return (abs(x_res) + abs(y_res)) / 2

    @cached_property
    def _scaling_factor(self) -> float:
        """Cache scaling factor calculation."""
        return self._data_resolution / self.min_distance_meters

    @cached_property
    def _power_in_region(self) -> xr.DataArray:
        """
        Cache power calculation from flow and temperature.

        Returns:
            xr.DataArray: Power in MW.
        """
        # Constants for power calculation
        density_water = 1000  # kg/m^3
        heat_capacity_water = 4.18  # kJ/kg/K
        mw_per_kw = 1 / 1000

        # Pre-calculate conversion factor
        conversion_factor = density_water * heat_capacity_water * mw_per_kw

        # Mean Volume flow for the area of interest
        usable_volume_flow = self.max_relative_volume_flow * self._volume_flow_in_region

        # Calculate temperature difference for approximation of the heat flow
        delta_t = (
            self._water_temperature_in_region - self.min_outlet_temperature
        ).clip(max=self.delta_t_max, min=0)

        # Calculate heat flow with single combined operation
        return usable_volume_flow * delta_t * conversion_factor

    @cached_property
    def _power_sum_spatial(self) -> xr.DataArray:
        """
        Cache the expensive spatial sum of power.

        Returns
        -------
        xr.DataArray
            Spatial sum of power over x and y dimensions
        """
        return self._power_in_region.sum(dim=["x", "y"])

    @cached_property
    def _power_sum_temporal(self) -> xr.DataArray:
        """
        Cache the expensive temporal sum of power.

        Returns
        -------
        xr.DataArray
            Temporal sum of power over time dimension
        """
        return self._power_in_region.sum(dim=[self.TIME])
