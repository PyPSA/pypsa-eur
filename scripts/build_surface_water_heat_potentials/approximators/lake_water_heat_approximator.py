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
        min_outlet_temperature: float = 1,
    ) -> None:
        self.ambient_temperature = ambient_temperature
        self.lake_shapes = lake_shapes
        self.region = region  # Set early so _lake_parts_in_region can access it

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
    def _round_coordinates(
        da: xr.DataArray, decimal_precision: int = 4
    ) -> xr.DataArray:
        return RiverWaterHeatApproximator._round_coordinates(da, decimal_precision)

    @staticmethod
    def _approximate_lake_temperature(**kwargs) -> xr.DataArray:
        return RiverWaterHeatApproximator._approximate_river_temperature(**kwargs)

    @cached_property
    def _total_volume_flow_m3_per_s(self) -> float:
        """
        Calculate total volume flow in m³/s from lake volume.

        Vol_total in HydroLAKES is in MCM (million cubic meters = 1e6 m³).
        """
        lake_volume_mcm = self._lake_parts_in_region["Vol_total"].sum()
        lake_volume_m3 = lake_volume_mcm * 1e6  # MCM to m³
        # Convert yearly turnover to per-second flow (m³/s)
        volume_flow_m3_per_s = lake_volume_m3 / (8760 * 3600)
        return volume_flow_m3_per_s

    @cached_property
    def _volume_flow_in_region(self) -> xr.DataArray:
        """
        Create a volume flow raster distributed uniformly over lake pixels.

        The total flow is divided by the number of valid pixels so that
        sum(flow_raster) = total_flow.
        """
        total_flow = self._total_volume_flow_m3_per_s

        # Use water temperature raster as template for the shape (before averaging)
        template = self._water_temperature_in_region_raster

        # Count valid (non-NaN) pixels in the raster
        # Use a 2D slice to count pixels (same count for all time steps)
        template_2d = template.isel(time=0) if "time" in template.dims else template
        n_valid_pixels = template_2d.notnull().sum().values

        if n_valid_pixels == 0:
            return xr.full_like(template, 0.0)

        # Distribute flow uniformly so sum = total_flow
        flow_per_pixel = total_flow / n_valid_pixels

        # Create raster with flow_per_pixel where template is valid, 0 elsewhere
        flow_raster = xr.where(template.notnull(), flow_per_pixel, 0.0)
        return flow_raster

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
        # Use x, y dims (EPSG:3035 projected) instead of lat, lon
        return self._water_temperature_in_region_raster.mean(dim=("x", "y"))
