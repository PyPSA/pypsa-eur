# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
from functools import cached_property
from typing import Union

import geopandas as gpd
import shapely
import xarray as xr

from scripts.build_surface_water_heat_potentials.approximators.surface_water_heat_approximator import (
    RiverWaterHeatApproximator,
    SurfaceWaterHeatApproximator,
)


class LakeWaterHeatApproximator(SurfaceWaterHeatApproximator):
    """
    River water heat approximator for district heating systems.

    Parameters are mostly based on expert input and Triebs 2023: "Untersuchung der zukünftigen Fernwärmeversorgung unter Unsicherheit bei Berücksichtigung technischer, ökonomischer und ökologischer Randbedingungen". # codespell:ignore unter

    Min_distance of 25km is based on Jung et al.: "Estimation of
    Temperature Recovery Distance and the Influence of Heat Pump Discharge on
    Fluvial Ecosystems".
    """

    def __init__(
        self,
        ambient_temperature: xr.DataArray,
        region: Union[shapely.geometry.polygon.Polygon, gpd.GeoSeries],
        lake_shapes: gpd.GeoDataFrame,
        delta_t_max: float = 1,
        min_outlet_temperature: float = 1,
    ) -> None:
        self.ambient_temperature = ambient_temperature
        self.lake_shapes = lake_shapes

        water_temperature = self._approximate_lake_temperature(
            ambient_temperature=ambient_temperature
        )
        water_temperature = water_temperature.rio.write_crs(
            f"EPSG:{ambient_temperature.rio.crs.to_epsg()}"
        )

        super().__init__(
            volume_flow=self._volume_flow_in_region,
            water_temperature=self._round_coordinates(water_temperature),
            region=region,
            max_relative_volume_flow=1.0,
            delta_t_max=delta_t_max,
            min_outlet_temperature=min_outlet_temperature,
        )

    @staticmethod
    def _round_coordinates(**kwargs) -> xr.DataArray:
        return RiverWaterHeatApproximator._round_coordinates(**kwargs)

    @staticmethod
    def _approximate_lake_temperature(**kwargs) -> xr.DataArray:
        return RiverWaterHeatApproximator._approximate_river_temperature(**kwargs)

    @cached_property
    def _volume_flow(self):
        raise NotImplementedError()

    @cached_property
    def _volume_flow_in_region(self):
        lake_volume_km3 = self._lake_parts_in_regions.groupby("name")["Vol_total"].sum()
        lake_volume_m3 = lake_volume_km3 * 1e9
        hourly_volume_m3 = lake_volume_m3 / 8760
        return hourly_volume_m3

    @cached_property
    def _lake_parts_in_region(self) -> gpd.GeoDataFrame:
        """Get lake parts within the defined region."""
        lake_parts = gpd.overlay(
            self.lake_shapes,
            gpd.GeoDataFrame(geometry=[self.region], crs=self.lake_shapes.crs),
            how="intersection",
        )
        return lake_parts

    def _air_temperature_in_lakes(
        self, eligible_lake_parts: gpd.GeoDataFrame
    ) -> xr.DataArray:
        """Get air temperature averaged over lake parts."""

        return self.ambient_temperature.rio.clip(
            eligible_lake_parts.geometry, eligible_lake_parts.crs, drop=True
        )

    @cached_property
    def _water_temperature_in_region_raster(self) -> gpd.GeoDataFrame:
        air_temperature_in_lakes = self._air_temperature_in_lakes(
            self._lake_parts_in_region
        )
        return self._approximate_lake_temperature(air_temperature_in_lakes)

    @cached_property
    def _water_temperature_in_region(self) -> xr.DataArray:
        return self._water_temperature_in_region_raster.mean(dim=("lat", "lon"))

    @cached_property
    def _power_sum_spatial(self) -> None:
        raise NotImplementedError("Lake power is not spatially resolved.")

    @cached_property
    def _power_sum_temporal(self) -> None:
        raise NotImplementedError("Lake power is not spatially resolved.")

    @cached_property
    def _power_in_region(self):
        raise NotImplementedError("Lake power is not spatially resolved.")
