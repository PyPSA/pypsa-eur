# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
from functools import cached_property
from typing import Union

import geopandas as gpd
import shapely
import xarray as xr

from scripts.build_surface_water_heat_potentials.approximators.river_water_heat_approximator import (
    RiverWaterHeatApproximator,
)
from scripts.build_surface_water_heat_potentials.approximators.surface_water_heat_approximator import (
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
        min_outlet_temperature: float = -100,
    ) -> None:
        self.ambient_temperature = ambient_temperature

        self.lake_shapes = lake_shapes
        lake_shapes.geometry = lake_shapes.geometry.map(
            lambda g: type(g)([p.exterior for p in g.geoms])
        ).simplify(  # drop holes
            0.001, preserve_topology=True
        )  # simplify (~100 m)

        self.region = region  # Set early so _lake_parts_in_region can access it

        water_temperature = self._approximate_lake_temperature(
            ambient_temperature=ambient_temperature
        )
        water_temperature = water_temperature.rio.write_crs(
            f"EPSG:{ambient_temperature.rio.crs.to_epsg()}"
        )

        super().__init__(
            volume_flow=self._volume_flow_in_region,
            water_temperature=water_temperature,
            region=region,
            max_relative_volume_flow=1.0,
            delta_t_max=delta_t_max,
            min_outlet_temperature=min_outlet_temperature,
        )

    @property
    def _scaling_factor(self) -> float:
        return 1.0  # No scaling needed for lakes

    @staticmethod
    def _approximate_lake_temperature(**kwargs) -> xr.DataArray:
        return RiverWaterHeatApproximator._approximate_river_temperature(**kwargs)

    @cached_property
    def _volume_flow_in_region(self) -> float:
        """
        Calculate total volume flow in m³/s from lake volume.

        Vol_total in HydroLAKES is in MCM (million cubic meters = 1e6 m³).
        """
        lake_volume_km3 = self._lake_parts_in_region["Vol_total"].sum()
        lake_volume_m3 = lake_volume_km3 * 1e6  # MCM to m³
        # Convert yearly turnover to per-hour flow (m³/h)
        volume_flow_m3_per_h = lake_volume_m3 / 8760
        return volume_flow_m3_per_h

    @cached_property
    def _lake_parts_in_region(self) -> gpd.GeoDataFrame:
        """Get lake parts within the defined region."""
        # Handle both GeoSeries and single geometry
        if isinstance(self.region, gpd.GeoSeries):
            region_gdf = gpd.GeoDataFrame(geometry=self.region, crs=self.region.crs)
        else:
            region_gdf = gpd.GeoDataFrame(
                geometry=[self.region], crs=self.lake_shapes.crs
            )
        lake_parts = gpd.overlay(
            self.lake_shapes,
            region_gdf,
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
    def _water_temperature_in_region_raster(self) -> xr.DataArray:
        air_temperature_in_lakes = self._air_temperature_in_lakes(
            self._lake_parts_in_region
        )
        return self._approximate_lake_temperature(
            ambient_temperature=air_temperature_in_lakes
        )

    @cached_property
    def _water_temperature_in_region(self) -> xr.DataArray:
        return self._water_temperature_in_region_raster.mean(dim=("x", "y"))

    @cached_property
    def _power_sum_spatial(self) -> xr.DataArray:
        """
        Cache the expensive spatial sum of power.

        Returns
        -------
        xr.DataArray
            Spatial sum of power over x and y dimensions
        """
        return self._power_in_region

    def get_spatial_aggregate(self) -> xr.Dataset:
        """
        Get the spatial aggregate of water temperature and power.

        Returns
        -------
        xr.Dataset
            Dataset containing total_power and average_temperature
        """
        # Data is already spatially aggregated
        # Combine into a single dataset
        return xr.Dataset(
            data_vars={
                "total_power": self._power_sum_spatial,
                "average_temperature": self._water_temperature_in_region,
            }
        )
