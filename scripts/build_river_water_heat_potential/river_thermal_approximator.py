import xarray as xr
import numpy as np
import shapely as sh
import geopandas as gpd
import pandas as pd

from functools import cached_property

METERS_PER_DEGREE: int = 111111


class RiverThermalApproximator:
    """
    Approximates river thermal properties based on geometry and HERA datasets.


    Parameters
    ----------
    geometry : geopandas.GeoDataFrame
        GeoDataFrame containing the geometry of the region.
        Must have exactly one row.
    hera_river_discharge : xarray.Dataset
        HERA dataset for river discharge.
    hera_ambient_temperature : xarray.Dataset
        HERA dataset for ambient temperature.
    target_snapshots : pandas.DatetimeIndex
        Time snapshots for resampling or indexing the data.
    delta_t_max : float, optional
        Maximum allowed cooling of the river at inlet point [K]. Default is 4.
    minimum_inlet_temperature : float, optional
        Minimum inlet temperature into the river [°C]. Default is 1.
    max_relative_river_volume_flow : float, optional
        Maximum relative allowed volume flow to be taken from a river (0-1). Default is 0.1.
    density_water : float, optional
        Density of water [kg/m³]. Default is 1000.
    heat_capacity_water : float, optional
        Heat capacity of water [kJ/(kg·K)]. Default is 4.18.
    moving_average_num_days : int, optional
        Number of days for the moving average of ambient temperature in river areas.
        Default is 13.
    dim_name_latitude : str, optional
        Name of the latitude dimension in HERA data. Default is "lat".
    dim_name_longitude : str, optional
        Name of the longitude dimension in HERA data. Default is "lon".
    dim_name_time : str, optional
        Name of the time dimension in HERA data. Default is "time".
    decimal_precision_coordinates : int, optional
        Precision for coordinate matching between HERA datasets. Default is 10.
    min_meters_between_heat_sources : float, optional
        Minimum spatial separation in meters between identified heat sources. Default is 6000.

    Attributes
    ----------
    total_river_power : xarray.DataArray
        Aggregated river power across the best locations, resampled to the target snapshots.
    mean_river_temperature : xarray.DataArray
        Weighted average river temperature at the best locations, resampled to the target snapshots.

    Methods
    -------
    approximate_river_power(river_discharge, river_temperature, max_relative_river_volume_flow, minimum_outlet_temperature, max_river_cooling, density_water=1000, heat_capacity_water=4.18) : xarray.DataArray
        Static method to estimate the heat flow (kW) that can be extracted from the river.
    get_hera_data_in_geometry(hera_data_array: xarray.DataArray) : xarray.DataArray
        Mask and return the input data that lies within the specified geometry.

    Raises
    ------
    ValueError
        If `geometry` is not a GeoDataFrame with exactly one row.
        If coordinates of `hera_river_discharge` and `hera_ambient_temperature` do not match.
    """

    def __init__(
        self,
        geometry: gpd.GeoDataFrame,
        hera_river_discharge: xr.Dataset,
        hera_ambient_temperature: xr.Dataset,
        target_snapshots: pd.DatetimeIndex,
        delta_t_max: float = 4,
        minimum_inlet_temperature: float = 1,
        max_relative_river_volume_flow: float = 0.1,
        density_water: float = 1000,
        heat_capacity_water: float = 4.18,
        moving_average_num_days: int = 13,
        dim_name_latitude: str = "lat",
        dim_name_longitude: str = "lon",
        dim_name_time: str = "time",
        decimal_precision_coordinates=10,
        min_meters_between_heat_sources: float = 6000,
    ) -> None:

        self._target_snapshots = target_snapshots
        self._geometry = geometry
        self._decimal_precision_coordinates = decimal_precision_coordinates
        self._dim_name_latitude = dim_name_latitude
        self._dim_name_longitude = dim_name_longitude
        self._dim_name_time = dim_name_time
        self._delta_t_max = delta_t_max
        self._minimum_inlet_temperature = minimum_inlet_temperature
        self._max_relative_river_volume_flow = max_relative_river_volume_flow
        self._density_water = density_water
        self._heat_capacity_water = heat_capacity_water
        self._moving_average_num_days = moving_average_num_days
        self._min_meters_between_heat_sources = min_meters_between_heat_sources

        self._hera_river_discharge = self._round_hera_coordinates(
            self._get_hera_data_in_box(hera_dataset=hera_river_discharge)
        )
        self._hera_ambient_temperature = self._round_hera_coordinates(
            self._get_hera_data_in_box(hera_dataset=hera_ambient_temperature)
        )
        self._validate_input()

    @property
    def input_resolution(self):
        return (
            self._hera_river_discharge[self._dim_name_longitude]
            .diff(self._dim_name_longitude)
            .mean()
            .values
            * METERS_PER_DEGREE
            + self._hera_river_discharge[self._dim_name_latitude]
            .diff(self._dim_name_latitude)
            .mean()
            .values
            * METERS_PER_DEGREE
        ) / 2

    def _validate_input(self):
        """
        Validates the input datasets and geometry.

        Raises
        ------
        ValueError
            If `geometry` is not a GeoDataFrame with exactly one row.
            If coordinates of `hera_river_discharge` and `hera_ambient_temperature` do not match.
            If dimensions of `hera_river_discharge` do not match self._dim_name_longitude, self._dim_name_latitude, self._dim_name_time.
        """
        if (
            not isinstance(self._geometry, gpd.GeoDataFrame)
            or len(self._geometry.index) != 1
        ):
            raise ValueError("Geometry has to be a GeoDataFrame with exactly one row")
        if not self._hera_ambient_temperature.coords.equals(
            self._hera_river_discharge.coords
        ):
            raise ValueError(
                "Coordinates of hera_river_discharge and hera_ambient_temperature do not match. Consider decreasing `decimal_precision_coordinates`."
            )

        if not {
            self._dim_name_time,
            self._dim_name_latitude,
            self._dim_name_longitude,
        }.issubset(self._hera_river_discharge.dims):
            raise ValueError(
                f"HERA data has unexpected dimension names: {self._hera_river_discharge.dims}. Expected dimensions: {self._dim_name_time}, {self._dim_name_latitude}, {self._dim_name_longitude}"
            )

    @cached_property
    def total_river_power(self):
        """
        Aggregate the river power across best locations.

        Best locations are coarsed HERA coordinates, whereas the coarsening is done to a resolution of `min_meters_between_heat_sources` x `min_meters_between_heat_sources` and for each coarse cell, the river power is the max. river power of those HERA locations in the respective coarse cell (c.f. `_best_locations()`).

        Returns
        -------
        xr.DataArray
            DataArray containing the river power at the specified location.
        """

        return (
            self._resample_to_target_index(
                self._river_power_in_geometry.sum(
                    dim=[self._dim_name_latitude, self._dim_name_longitude]
                )
            )
            * self.input_resolution
            / self._min_meters_between_heat_sources
        )

    @cached_property
    def mean_river_temperature(self) -> xr.DataArray:
        """Calculate the mean river temperature in `geometry` at `self._best_locations`, weighted by river power in best locations.

        Returns
        -------
        xr.DataArray
            Series containing the mean river temperature.
        """

        return self._resample_to_target_index(
            (self._river_temperature_in_geometry * self._river_power_in_geometry).sum(
                dim=[self._dim_name_latitude, self._dim_name_longitude]
            )
            / (
                self._river_power_in_geometry.sum(
                    dim=[self._dim_name_latitude, self._dim_name_longitude]
                )
                + 0.001
            )
        )

    @cached_property
    def _river_power_in_geometry(self) -> xr.DataArray:
        """Calculate the river power in `geometry`.

        Approximation based on the mean river discharge and temperature using :meth:`_approximamate_mean_river_power.

        Returns
        -------
        xr.DataArray
           xr.DataArray containing the mean river power.
        """
        return self._approximate_river_power(
            river_discharge=self._river_discharge_in_geometry,
            river_temperature=self._river_temperature_in_geometry,
            max_relative_river_volume_flow=self._max_relative_river_volume_flow,
            minimum_outlet_temperature=self._minimum_inlet_temperature,
            max_river_cooling=self._delta_t_max,
        )

    @cached_property
    def _river_temperature_in_geometry(self):
        """
        Load and process the HERA-dataset with ambient temperature to get the river temperature for a given time

        Returns
        -------
        xr.DataArray
            DataArray containing the river temperature.
        """

        return self._approximate_river_temperature(
            ambient_temperature_moving_average=self._ambient_temperature_in_rivers_moving_average
        )

    @cached_property
    def _ambient_temperature_in_rivers_moving_average(self) -> xr.DataArray:
        """Calculate the moving average of the ambient temperature in river areas in `geometry` using :meth:`_get_moving_average`.

        Returns
        -------
        xr.DataArray
            DataArray containing the moving average of the ambient temperature in river areas.
        """
        return xr.DataArray(
            self._get_moving_average(da=self._ambient_temperature_in_rivers),
            name="river_temperature_moving_average",
        )

    @cached_property
    def _ambient_temperature_in_rivers(self) -> xr.DataArray:
        """Return the ambient temperature in river areas in `geometry`.

        Returns
        -------
        xr.DataArray
            DataArray containing the ambient temperature in river areas."""
        return xr.DataArray(
            self._ambient_temperature_in_geometry.where(
                self._river_discharge_in_geometry.notnull()
            ),
            name="ambient_temperature_in_river",
        )

    @cached_property
    def _ambient_temperature_in_geometry(self) -> xr.DataArray:
        """Return the ambient temperature in `geometry`.

        Returns
        -------
        xr.DataArray
            DataArray containing the ambient temperature.
        """
        return xr.DataArray(
            self._get_hera_data_in_geometry(
                hera_data_array=self._hera_ambient_temperature
            ).ta6,
            name="ambient_temperature_in_river",
        )

    @cached_property
    def _river_discharge_in_geometry(self) -> xr.DataArray:
        """Return the ambient temperature in `geometry`.

        Returns
        -------
        xr.DataArray
            DataArray containing the ambient temperature.
        """
        return xr.DataArray(
            self._get_hera_data_in_geometry(
                hera_data_array=self._hera_river_discharge
            ).dis,
            name="river_discharge",
        )

    @cached_property
    def geometry(self) -> sh.geometry:
        """
        Extract the geometry of the region.

        Returns
        -------

        sh.geometry
            Geometry of the region.
        """

        return self._geometry.geometry.iloc[0]

    @cached_property
    def _bounds_x(self):
        """
        Extract bounding x-coordinates from the geographical boundaries of `geometry`, that create a bounding box.

        Returnt
        -------
        slice
            Bounding x-coordinates.
        """

        return slice(self._geometry.bounds.minx.min(), self._geometry.bounds.maxx.max())

    @cached_property
    def _bounds_y(self):
        """
        Extract bounding y-coordinates from the geographical boundaries of `geometry`, that create a bounding box.

        Returns
        -------
        slice
            Bounding y-coordinates.
        """

        return slice(self._geometry.bounds.miny.min(), self._geometry.bounds.maxy.max())

    @staticmethod
    def _approximate_river_temperature(
        ambient_temperature_moving_average: xr.DataArray,
        k1=-0.957,
        k2=28.212,
        k3=12.434,
        k4=0.137,
    ) -> xr.DataArray:
        """
        Apply the formula for derivation of the river temperature from the ambient temperature (Triebs & Tsatsaronis 2022: Estimating the local renewable potentials for the transformation of district heating systems, ECOS 2022, pp. 479-490: https://orbit.dtu.dk/en/publications/proceedings-of-ecos-2022-the-35th-international-conference-on-eff)

        Parameters:
        -----------
        ambient_temperature_moving_average: xr.DataArray
            DataArray containing the moving average of the ambient temperature in river areas.
        k1, k2, k3, k4: float
            Regression coefficients for the approximation of the river temperature.
        """
        return k1 + (k2 / (1 + np.exp(k4 * (k3 - ambient_temperature_moving_average))))

    @staticmethod
    def _approximate_river_power(
        river_discharge: xr.DataArray,
        river_temperature: xr.DataArray,
        max_relative_river_volume_flow: float,
        minimum_outlet_temperature: float,
        max_river_cooling: float,
        density_water=1000,
        heat_capacity_water=4.18,
    ) -> xr.DataArray:
        """
        Calculate the mean heat flow in kW, that can be extracted from the river. Several boundary conditions have to be specifided.
        The heat flow will be averaged over the whole area of interest.

        Parameters
        ----------
        river_discharge: xr.DataArray
            DataArray containing the river discharge in m³/s.
        river_temperature: xr.DataArray
            DataArray containing the river temperature in °C.
        max_relative_river_volume_flow: float
            Maximum relative volume flow to be taken from the river.
        minimum_outlet_temperature: float
            Minimum inlet temperature into the river in °C.
        max_river_cooling: float
            Maximum allowed cooling of the river at the inlet point in K.
        density_water : float, optional
            Density of water [kg/m³], by default 1000.
        heat_capacity_water : float, optional
            Specific heat capacity of water [kJ/(kg·K)], by default 4.18.

        Returns
        -------
        xr.DataArray
            DataArray containing the usable mean heat flow in kW.
        """

        # Mean Volume flow for the area of interest
        v_dot = max_relative_river_volume_flow * river_discharge

        # Mean river temperature for the area of interest
        t_source = river_temperature

        # Calculate temperature difference for approximation of the heat flow
        delta_t = (t_source - minimum_outlet_temperature).clip(max=max_river_cooling)

        # Calculate heat flow
        return v_dot * density_water * heat_capacity_water * delta_t

    def _get_hera_data_in_geometry(self, hera_data_array: xr.DataArray) -> xr.DataArray:
        """
        Mask the input dataset for data, that lies within a given boundary.

        The boundary is given in the shape of a polygon and is more complex than the bounding box, that has been applied before.
        Therefore a mask with boolean values is used to find all the data, that lies within the boundary-polygon.

        Parameters:
        -----------
        hera_dataset: xarray.Dataset
            HERA dataset for river discharge, ambient or river temperature.

        Returns:
        --------
        xarray.Dataset
            Dataset containing the masked HERA data.
        """

        # Extract the coordinates from the xarray dataset
        lats = hera_data_array[self._dim_name_latitude].values
        lons = hera_data_array[self._dim_name_longitude].values

        # Create a 2D-Grid with all combinations of lat and lon
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        # Convert the Grid into a list of points
        points = np.vstack([lon_grid.ravel(), lat_grid.ravel()]).T

        # Create a GeoSeries out of the points
        # TODO: use geopandas points_from_xy or similar
        gdf_points = gpd.GeoSeries([sh.geometry.Point(xy) for xy in points])

        # Check which points lie within the boundary-polygon
        mask = gdf_points.within(self.geometry)

        # reshape the mask into the shape of the Grid
        mask_2d = mask.values.reshape(lat_grid.shape)

        # Convert teh mask into an xarray-DataArray
        mask_array = xr.DataArray(
            mask_2d,
            dims=[self._dim_name_latitude, self._dim_name_longitude],
            coords={self._dim_name_latitude: lats, self._dim_name_longitude: lons},
        )

        # Filter the dataframe. The filtered Dataframe will only contain data for the Area Of interest (AOI), defined by polygonial boundaries
        masked_dataset = hera_data_array.where(mask_array)

        # round coordinates for spatial matching
        masked_dataset.coords[self._dim_name_latitude] = masked_dataset.coords[
            self._dim_name_latitude
        ].round(self._decimal_precision_coordinates)
        masked_dataset.coords[self._dim_name_longitude] = masked_dataset.coords[
            self._dim_name_longitude
        ].round(self._decimal_precision_coordinates)

        return masked_dataset  # .compute()

    def _get_hera_data_in_box(self, hera_dataset: xr.Dataset) -> xr.Dataset:
        """
        Mask HERA data for given boundaries in time and space in a square bounding box.

        Parameters:
        -----------
        hera_dataset: xarray.Dataset
            HERA dataset for river discharge or ambient temperature.

        Returns:
        --------
        xarray.Dataset
            Dataset containing the HERA data for the given boundaries.
        """
        # Select the data for the given boundaries
        hera_dataset = hera_dataset.sortby(
            [self._dim_name_time, self._dim_name_latitude, self._dim_name_longitude]
        )
        return hera_dataset.sel(lon=self._bounds_x, lat=self._bounds_y)

    def _get_moving_average(self, da: xr.DataArray) -> xr.DataArray:
        """
        Calculate the moving average for the air temperature over a given number of days.

        Parameters:
        -----------
        da: xarray.DataArray
        num_days: int
            Number of days to be averaged about. Default: num_days = 13.

        Returns:
        --------
        xr.DataArray
            DataArray containing the moving average `da`.

        """

        # Time resolution (hours per time step)
        time_res = int(da.time.dt.hour.frequency)
        # Time steps per day
        t_step_day = 24 / time_res
        # Window for moving average
        window = int(t_step_day * self._moving_average_num_days)
        # Calculate the mean ambient temperature
        return da.rolling(time=window, min_periods=1).mean()

    def _resample_to_target_index(self, data: xr.DataArray) -> xr.DataArray:
        """
        Resample and interpolate the input data to match the target time index.

        Parameters
        ----------
        data : xr.DataArray
            The input data array with a time dimension.

        Returns
        -------
        xr.DataArray
            The data resampled and interpolated to match the target time index.

        Raises
        ------
        ValueError
            If more than one dimension is present in `data`.
        """
        if len(data.dims) > 1:
            raise ValueError("Data array must have only one dimension")
        # Ensure data is sorted by its index (assuming 'time' is the index)
        data = data.sortby(self._dim_name_time)

        # Reindex to the target time index and interpolate
        resampled_series = (
            data.to_pandas()
            .reindex(self._target_snapshots, method="nearest")
            .interpolate()
        )

        # Convert back to xarray
        return xr.DataArray(
            resampled_series.values,
            dims=[self._dim_name_time],
            coords={self._dim_name_time: self._target_snapshots},
        )

    def _round_hera_coordinates(self, hera_data_array: xr.DataArray) -> xr.DataArray:
        """
        Round the coordinates of the HERA dataset to the defined precision.

        Parameters
        ----------
        hera_data_array : xr.DataArray
            HERA dataset.

        Returns
        -------
        xr.DataArray
            HERA dataset with rounded coordinates.
        """
        hera_data_array.coords[self._dim_name_latitude] = hera_data_array.coords[
            self._dim_name_latitude
        ].round(self._decimal_precision_coordinates)
        hera_data_array.coords[self._dim_name_longitude] = hera_data_array.coords[
            self._dim_name_longitude
        ].round(self._decimal_precision_coordinates)
        return hera_data_array
