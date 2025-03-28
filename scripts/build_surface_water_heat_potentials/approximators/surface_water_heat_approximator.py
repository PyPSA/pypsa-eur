from abc import ABC

import numpy as np
import shapely
import shapely.vectorized as sv
import xarray as xr


class SurfaceWaterHeatApproximator(ABC):
    """
    A class for calculating heat source potential for seawater-sourced heat pumps.

    This class encapsulates the full workflow for loading oceanographic data,
    calculating heat potential, resampling, masking to regions of interest,
    and visualizing results.

    Attributes:
        results (xr.Dataset): The results of the calculation. Contains `total_power` and `average_temperature`. Coordinates are `time`.
    """

    METERS_PER_DEGREE: int = 111111
    LONGITUDE = "longitude"
    LATITUDE = "latitude"

    def __init__(
        self,
        volume_flow: xr.DataArray,
        water_temperature: xr.DataArray,
        region_geometry: shapely.geometry.polygon.Polygon,
        max_relative_volume_flow: float = 1.0,
        delta_t_max: float = 4,
        min_outlet_temperature: float = 1,
        min_distance_meters: int = 2000,
    ):
        """
        Initialize the SeawaterThermalApproximator. This is an abstract class and should not be instantiated directly.

        Args:
            volume_flow (xr.DataArray): The volume flow of the water.
            water_temperature (xr.DataArray): The water temperature data.
            region_geometry (shapely.geometry.polygon.Polygon): The geometry of the region of interest.
            max_relative_volume_flow (float): The maximum relative volume flow.
            density_water (float): The density of water.
            heat_capacity_water (float): The heat capacity of water.
            delta_t_max (float): The maximum temperature difference.
            min_outlet_temperature (float): The minimum outlet temperature.
            min_distance_meters (int): The minimum distance in meters between two projects.
        """

        self.volume_flow = volume_flow
        self.water_temperature = water_temperature
        self.region_geometry = region_geometry
        self.max_relative_volume_flow = max_relative_volume_flow
        self.delta_t_max = delta_t_max
        self.min_outlet_temperature = min_outlet_temperature
        self.min_distance_meters = min_distance_meters

        self._validate_input(
            volume_flow=volume_flow, water_temperature=water_temperature
        )

        self.results = self.get_results()

    def get_results(self):
        """Do some computations and return a result."""

        boxed_volume_flow = self._get_boxed_data(data=self.volume_flow)
        boxed_water_temperature = self._get_boxed_data(data=self.water_temperature)

        mask = self.get_geometry_mask(data=boxed_volume_flow)

        masked_volume_flow = boxed_volume_flow.where(mask)
        masked_water_temperature = boxed_water_temperature.where(mask)

        masked_power = self.get_power(
            volume_flow=masked_volume_flow,
            temperature=masked_water_temperature,
        )

        total_power = masked_power.sum(
            dim=[self.LATITUDE, self.LONGITUDE]
        ) * self._get_scaling_factor(data=self.volume_flow)

        # Calculate power-weighted average temperature
        average_water_temperature = (masked_water_temperature * masked_power).sum(
            dim=[self.LATITUDE, self.LONGITUDE]
        ) / (masked_power.sum(dim=[self.LATITUDE, self.LONGITUDE]) + 0.001)

        # Combine into a single dataset
        return xr.Dataset(
            data_vars={
                "total_power": total_power,
                "average_temperature": average_water_temperature,
            }
        )

    def _validate_input(
        self, volume_flow: xr.DataArray, water_temperature: xr.DataArray
    ):
        if not volume_flow.coords.equals(water_temperature.coords):
            raise ValueError(
                "Coordinates of `volume_flow` and `temperature` do not match."
            )

    def _get_data_resolution(self, data: xr.DataArray) -> float:
        return (
            data[self.LONGITUDE].diff(self.LONGITUDE).mean().values
            * self.METERS_PER_DEGREE
            + data[self.LATITUDE].diff(self.LATITUDE).mean().values
            * self.METERS_PER_DEGREE
        ) / 2

    def _get_scaling_factor(self, data: xr.DataArray) -> float:
        return self._get_data_resolution(data) / self.min_distance_meters

    def _get_boxed_data(self, data: xr.DataArray) -> xr.Dataset:
        """
        Get the dataset boxed to the geometry.

        Args:
            data (xr.DataArray): Data to box.
            geometry (gpd.GeoDataFrame): The geometry.

        Returns:
            xr.Dataset: The boxed dataset.
        """

        # bound data to onshore region
        return data.sel(
            **{
                self.LONGITUDE: slice(
                    self.region_geometry.bounds[0], self.region_geometry.bounds[2]
                ),
                self.LATITUDE: slice(
                    self.region_geometry.bounds[1], self.region_geometry.bounds[3]
                ),
            }
        )

    def get_power(
        self,
        volume_flow: xr.DataArray,
        temperature: xr.DataArray,
        density_water: float = 1000,
        heat_capacity_water: float = 4.18,
    ) -> xr.DataArray:
        """
        Get the power from flow and temperature.

        Args:
            volume_flow (xr.DataArray): The volume flow.
            temperature (xr.DataArray): The temperature.
            max_relative_volume_flow (float): The maximum relative volume flow.
            density_water (float): The density of water.
            heat_capacity_water (float): The heat capacity of water.
            delta_t_max (float): The maximum temperature difference.
            min_outlet_temperature (float): The minimum outlet temperature.

        Returns:
            xr.DataArray: Power.
        """
        # Mean Volume flow for the area of interest
        usable_volume_flow = self.max_relative_volume_flow * volume_flow
        # Calculate temperature difference for approximation of the heat flow
        delta_t = (temperature - self.min_outlet_temperature).clip(max=self.delta_t_max)
        # Calculate heat flow
        return usable_volume_flow * density_water * heat_capacity_water * delta_t

    def get_geometry_mask(
        self,
        data: xr.DataArray,
    ) -> xr.DataArray:
        """
        Get the mask for the geometry border.

        Args:
            ds (xr.Dataset): The dataset.
            geometry: The geometry.
            scale_buffer (float): The scale buffer.

        Returns:
            xr.DataArray: The mask.
        """

        # Extract coordinate values from ds (note: coordinate names match those in ds)
        lon2d, lat2d = np.meshgrid(data[self.LONGITUDE], data[self.LATITUDE])
        # Create a boolean mask for grid points within the border buffer.
        mask = sv.contains(self.region_geometry, lon2d, lat2d)

        # Convert to an xarray DataArray with matching dims and coords.
        return xr.DataArray(
            mask,
            dims=[self.LATITUDE, self.LONGITUDE],
            coords={
                self.LATITUDE: data[self.LATITUDE],
                self.LONGITUDE: data[self.LONGITUDE],
            },
        )
