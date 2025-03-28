import xarray as xr
import geopandas as gpd
import numpy as np
import shapely.vectorized as sv
import shapely


METERS_PER_DEGREE: int = 111111
LATITUDE = "latitude"
LONGITUDE = "longitude"
TIME = "time"
TEMPERATURE = "thetao"
EASTWARD_VELOCITY = "uo"
NORTHWARD_VELOCITY = "vo"
NORMED_VELOCITY = "normed_velocity"
POWER = "power"


class SeawaterThermalApproximator:
    """
    A class for calculating heat source potential for seawater-sourced heat pumps.

    This class encapsulates the full workflow for loading oceanographic data,
    calculating heat potential, resampling, masking to regions of interest,
    and visualizing results.

    Attributes:
        results (xr.Dataset): The results of the calculation. Contains `total_power` and `average_temperature`. Coordinates are `time`.
    """

    def __init__(
        self,
        data: xr.Dataset,
        geometry: shapely.geometry.polygon.Polygon,
        max_relative_volume_flow: float = 1.0,
        density_water: float = 1000,
        heat_capacity_water: float = 4.18,
        delta_t_max: float = 4,
        min_outlet_temperature: float = 1,
        min_distance_meters: int = 2000,
    ):
        """
        Initialize the SeawaterThermalApproximator.

        Args:
            data (xr.Dataset): The raw oceanographic dataset.
            geometry (gpd.GeoDataFrame): The geometry of the region of interest.
            max_relative_volume_flow (float): The maximum relative volume flow.
            density_water (float): The density of water.
            heat_capacity_water (float): The heat capacity of water.
            delta_t_max (float): The maximum temperature difference.
            min_outlet_temperature (float): The minimum outlet temperature.
            min_distance_meters (int): The minimum distance in meters between two projects.
        """

        data = self._get_boxed_dataset(data, geometry)
        data[NORMED_VELOCITY] = self._get_normed_velocity(
            data[EASTWARD_VELOCITY], data[NORTHWARD_VELOCITY]
        )
        data[POWER] = self.get_power(
            data,
            max_relative_volume_flow,
            density_water,
            heat_capacity_water,
            delta_t_max,
            min_outlet_temperature,
        )
        data = data[[POWER, TEMPERATURE]]

        mask = self.get_mask_geometry_border(
            ds=data, geometry=geometry
        )

        masked_power = data[POWER].where(mask)
        masked_temp = data[TEMPERATURE].where(mask)

        total_power = masked_power.sum(dim=["latitude", "longitude", "depth"])

        # Scale by minimum distance between projects
        resolution_data = (data.longitude.diff("longitude").mean().values * METERS_PER_DEGREE + data.latitude.diff("latitude").mean().values * METERS_PER_DEGREE) / 2
        scaling_factor = resolution_data / min_distance_meters
        total_power = total_power * scaling_factor


        # Calculate power-weighted average temperature
        weighted_temp = (masked_temp * masked_power).sum(
            dim=["latitude", "longitude", "depth"]
        ) / (masked_power.sum(dim=["latitude", "longitude", "depth"]) + 0.001)

        # Combine into a single dataset
        self.results = xr.Dataset(
            data_vars={"total_power": total_power, "average_temperature": weighted_temp}
        )

    @staticmethod
    def _get_normed_velocity(
        eastward_velocity: xr.DataArray, northward_velocity: xr.DataArray
    ) -> xr.DataArray:
        """
        Calculate the normed velocity.

        Args:
            eastward_velocity (xr.DataArray): The eastward velocity.
            northward_velocity (xr.DataArray): The northward velocity.

        Returns:
            xr.DataArray: The normed velocity.
        """
        return np.sqrt(eastward_velocity**2 + northward_velocity**2)

    @staticmethod
    def _get_boxed_dataset(
        ds: xr.Dataset, geometry: shapely.geometry.polygon.Polygon
    ) -> xr.Dataset:
        """
        Get the dataset boxed to the geometry.

        Args:
            ds (xr.Dataset): The dataset.
            geometry (gpd.GeoDataFrame): The geometry.

        Returns:
            xr.Dataset: The boxed dataset."""
        # bound ds to onshore regions
        return ds.sel(
            longitude=slice(geometry.bounds[0], geometry.bounds[2]),
            latitude=slice(geometry.bounds[1], geometry.bounds[3]),
        )

    @staticmethod
    def get_power(
        ds: xr.Dataset,
        max_relative_volume_flow: float,
        density_water: float,
        heat_capacity_water: float,
        delta_t_max: float,
        min_outlet_temperature: float,
    ) -> xr.DataArray:
        """
        Get the power from flow and temperature.

        Args:
            ds (xr.Dataset): The dataset.
            max_relative_volume_flow (float): The maximum relative volume flow.
            density_water (float): The density of water.
            heat_capacity_water (float): The heat capacity of water.
            delta_t_max (float): The maximum temperature difference.
            min_outlet_temperature (float): The minimum outlet temperature.

        Returns:
            xr.DataArray: Power.
        """
        # normalise flow velocity
        ds[NORMED_VELOCITY] = np.sqrt(ds["uo"] ** 2 + ds["vo"] ** 2)
        # Mean Volume flow for the area of interest
        volume_flow = max_relative_volume_flow * ds[NORMED_VELOCITY]
        # Calculate temperature difference for approximation of the heat flow
        delta_t = (ds[TEMPERATURE] - min_outlet_temperature).clip(max=delta_t_max)
        # Calculate heat flow
        return volume_flow * density_water * heat_capacity_water * delta_t

    @staticmethod
    def get_mask_geometry_border(
        ds: xr.Dataset,
        geometry: shapely.geometry.polygon.Polygon,
        scale_buffer: float = 0.5,
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

        # get resolution
        resolution_ds = ds["latitude"].values[1] - ds["latitude"].values[0]
        # Extract coordinate values from ds (note: coordinate names match those in ds)
        lon2d, lat2d = np.meshgrid(ds["longitude"], ds["latitude"])
        # Create a buffer equal to resolution of ds around the region's border
        border_buffer = geometry.boundary.buffer(resolution_ds * scale_buffer)

        # Create a boolean mask for grid points within the border buffer.
        mask = sv.contains(border_buffer, lon2d, lat2d)

        # Convert to an xarray DataArray with matching dims and coords.
        return (
            xr.DataArray(
                mask,
                dims=["latitude", "longitude"],
                coords={"latitude": ds["latitude"], "longitude": ds["longitude"]},
            )
            * 1.0
        )

    @staticmethod
    def get_total_power(da_power_in_region: xr.DataArray) -> xr.DataArray:
        """
        Get the total power.

        Args:
            da_power_in_region (xr.DataArray): The power in the region.

        Returns:
            xr.DataArray: The total power.

        """
        return da_power_in_region.sum(dim=["latitude", "longitude"])

    @staticmethod
    def get_average_temperature(
        da_temperature_in_region: xr.DataArray, da_power_in_region: xr.DataArray
    ) -> xr.DataArray:
        """
        Get the average temperature.

        Args:
            da_temperature_in_region (xr.DataArray): The temperature in the region.
            da_power_in_region (xr.DataArray): The power in the region.

        Returns:
            xr.DataArray: The average temperature.
        """
        return (
            da_temperature_in_region
            * da_power_in_region.sum(dim=["time", "latitude", "longitude"])
        ).sum(dim=["latitude", "longitude"]) / da_power_in_region.sum(
            dim=["time", "latitude", "longitude"]
        )
