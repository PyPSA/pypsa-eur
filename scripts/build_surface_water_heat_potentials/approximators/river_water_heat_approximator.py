# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
import warnings
from typing import Union

import geopandas as gpd
import numpy as np
import shapely
import xarray as xr

from scripts.build_surface_water_heat_potentials.approximators.surface_water_heat_approximator import (
    SurfaceWaterHeatApproximator,
)


class RiverWaterHeatApproximator(SurfaceWaterHeatApproximator):
    """
    River water heat approximator for district heating systems.

    Parameters are mostly based on expert input and Triebs 2023: "Untersuchung der zukünftigen Fernwärmeversorgung unter Unsicherheit bei Berücksichtigung technischer, ökonomischer und ökologischer Randbedingungen". # codespell:ignore unter

    Min_distance of 25km is based on Jung et al.: "Estimation of
    Temperature Recovery Distance and the Influence of Heat Pump Discharge on
    Fluvial Ecosystems".
    """

    def __init__(
        self,
        volume_flow: xr.DataArray,
        ambient_temperature: xr.DataArray,
        region: Union[shapely.geometry.polygon.Polygon, gpd.GeoSeries],
        max_relative_volume_flow: float = 1.0,
        delta_t_max: float = 1,
        min_outlet_temperature: float = 1,
        min_distance_meters: int = 25000,
    ) -> None:
        water_temperature = self._approximate_river_temperature(
            ambient_temperature=ambient_temperature
        )
        water_temperature = water_temperature.rio.write_crs(
            f"EPSG:{ambient_temperature.rio.crs.to_epsg()}"
        )

        super().__init__(
            volume_flow=self._round_coordinates(volume_flow),
            water_temperature=self._round_coordinates(water_temperature),
            region=region,
            max_relative_volume_flow=max_relative_volume_flow,
            delta_t_max=delta_t_max,
            min_outlet_temperature=min_outlet_temperature,
            min_distance_meters=min_distance_meters,
        )

    def _round_coordinates(
        self, da: xr.DataArray, decimal_precision: int = 4
    ) -> xr.DataArray:
        """
        Round the coordinates of the HERA dataset to the defined precision.

        Parameters
        ----------
        da : xr.DataArray
            HERA dataset.

        Returns
        -------
        xr.DataArray
            HERA dataset with rounded coordinates.
        """
        if "x" in da.coords and "y" in da.coords:
            return da.assign_coords(
                x=da.x.round(decimal_precision), y=da.y.round(decimal_precision)
            )
        else:
            raise ValueError(
                "The DataArray does not contain the expected coordinates 'x' and 'y'."
            )

    @staticmethod
    def _approximate_river_temperature(
        ambient_temperature: xr.DataArray,
        moving_average_num_days: int = 13,
        k1: float = -0.957,
        k2: float = 28.212,
        k3: float = 12.434,
        k4: float = 0.137,
    ) -> xr.DataArray:
        """
        Apply the formula for derivation of the river temperature from the ambient temperature.

        Based on Triebs & Tsatsaronis 2022: Estimating the local renewable potentials
        for the transformation of district heating systems, ECOS 2022, pp. 479-490.

        Parameters
        ----------
        ambient_temperature : xr.DataArray
            DataArray containing ambient temperature in river areas
        moving_average_num_days : int, optional
            Number of days for moving average, by default 13
        k1, k2, k3, k4 : float, optional
            Regression coefficients for the approximation of the river temperature

        Returns
        -------
        xr.DataArray
            Approximated river temperature
        """
        # Time steps per day
        time_steps_per_day = int(24 / float(ambient_temperature.time.dt.hour.frequency))
        # Window for moving average
        window = time_steps_per_day * moving_average_num_days
        if window > len(ambient_temperature.time):
            window = len(ambient_temperature.time)
            warnings.warn(
                f"Moving average window of {moving_average_num_days} days in river water temperature approximation exceeds the available time steps ({len(ambient_temperature.time)}). Falling back to the maximum available time steps ({window} hours)."
            )
        # Calculate the mean ambient temperature
        ambient_temperature_moving_average = ambient_temperature.rolling(
            time=window, min_periods=1
        ).mean()

        return k1 + (k2 / (1 + np.exp(k4 * (k3 - ambient_temperature_moving_average))))
