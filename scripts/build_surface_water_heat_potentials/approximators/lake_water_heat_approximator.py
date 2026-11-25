# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Lake water heat approximator for district heating systems."""

from functools import cached_property

import geopandas as gpd
import xarray as xr

from scripts.build_surface_water_heat_potentials.approximators.river_water_heat_approximator import (
    RiverWaterHeatApproximator,
)
from scripts.build_surface_water_heat_potentials.approximators.surface_water_heat_approximator import (
    SurfaceWaterHeatApproximator,
)

SECONDS_PER_YEAR = 8760 * 3600


class LakeWaterHeatApproximator(SurfaceWaterHeatApproximator):
    """
    Approximate the heat that lakes can supply to district heating.

    The lake volume inside the region is cooled by at most `delta_t_max` per
    year, following Triebs (2023). A lake that lies partly inside the region
    contributes the share of its volume that equals its share of area inside
    the region. The water temperature follows the approximation for rivers
    by Triebs & Tsatsaronis (2022).

    Parameters
    ----------
    ambient_temperature : xr.DataArray
        Ambient air temperature in EPSG:3035.
    region : gpd.GeoSeries
        Region of interest in EPSG:3035.
    lake_shapes : gpd.GeoDataFrame
        Lake polygons from HydroLAKES in EPSG:3035 with `Vol_total` in million m³.
    delta_t_max : float, optional
        Maximum cooling of the lake volume per year in K.
    min_outlet_temperature : float, optional
        Minimum outlet water temperature in °C.
    """

    def __init__(
        self,
        ambient_temperature: xr.DataArray,
        region: gpd.GeoSeries,
        lake_shapes: gpd.GeoDataFrame,
        delta_t_max: float = 1,
        min_outlet_temperature: float = 1,
    ) -> None:
        self.lake_shapes = lake_shapes
        super().__init__(
            volume_flow=xr.zeros_like(ambient_temperature),  # not used
            water_temperature=ambient_temperature,
            region=region,
            delta_t_max=delta_t_max,
            min_outlet_temperature=min_outlet_temperature,
        )

    @cached_property
    def _volume_flow_in_region(self) -> float:
        """Volume flow in m³/s that cools the lake volume in the region once per year."""
        lakes = self.lake_shapes.assign(lake_area=self.lake_shapes.area)
        region = gpd.GeoDataFrame(geometry=self.region.geometry)
        parts = gpd.overlay(lakes, region, how="intersection")
        volume = (parts["Vol_total"] * parts.area / parts["lake_area"]).sum() * 1e6
        return volume / SECONDS_PER_YEAR

    @cached_property
    def _water_temperature_in_region(self) -> xr.DataArray:
        """Mean water temperature over time, approximated from the air temperature in the region."""
        air_temperature = self.water_temperature.rio.clip(
            self.region.geometry.buffer(self._data_resolution), drop=True
        )
        return RiverWaterHeatApproximator._approximate_river_temperature(
            ambient_temperature=air_temperature
        ).mean(dim=("x", "y"))

    def get_spatial_aggregate(self) -> xr.Dataset:
        """Get total power in MW and water temperature in °C over time."""
        return xr.Dataset(
            {
                "total_power": self._power_in_region,
                "average_temperature": self._water_temperature_in_region,
            }
        )
