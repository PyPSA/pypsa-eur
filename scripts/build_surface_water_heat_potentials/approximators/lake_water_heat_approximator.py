# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Lake water heat approximator for district heating systems."""

from functools import cached_property

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
    Lake water heat approximator for district heating systems.

    Estimates the thermal potential of lakes as a heat source for district
    heating applications. Uses lake volume data from HydroLAKES and ambient
    temperature from HERA to compute available heating power.

    The water temperature approximation follows the methodology from
    Triebs & Tsatsaronis 2022 (originally developed for rivers).

    Parameters
    ----------
    ambient_temperature : xr.DataArray
        Ambient air temperature data (HERA), used to approximate water temperature.
    region : Union[shapely.geometry.polygon.Polygon, gpd.GeoSeries]
        Geographic region of interest.
    lake_shapes : gpd.GeoDataFrame
        Lake polygons from HydroLAKES with Vol_total attribute.
    delta_t_max : float, optional
        Maximum temperature difference for heat extraction, by default 1 K.
    min_outlet_temperature : float, optional
        Minimum outlet water temperature, by default 1 °C.

    References
    ----------
    .. [1] Messager et al. (2016). Estimating the volume and age of water
           stored in global lakes. Nature Communications, 13603.
    .. [2] Triebs & Tsatsaronis (2022). Estimating the local renewable
           potentials for the transformation of district heating systems.
    """

    def __init__(
        self,
        ambient_temperature: xr.DataArray,
        region: shapely.geometry.polygon.Polygon | gpd.GeoSeries,
        lake_shapes: gpd.GeoDataFrame,
        delta_t_max: float = 1,
        min_outlet_temperature: float = 1,
    ) -> None:
        self.ambient_temperature = ambient_temperature
        self.lake_shapes = self._simplify_lake_geometries(lake_shapes)
        self.region = region

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
        """Scaling factor for power calculation (disabled for lakes)."""
        return 1.0

    @staticmethod
    def _simplify_lake_geometries(
        lake_shapes: gpd.GeoDataFrame, tolerance: float = 0.001
    ) -> gpd.GeoDataFrame:
        """
        Simplify lake polygon geometries for faster spatial operations.

        Removes interior holes and simplifies exterior boundaries.

        Parameters
        ----------
        lake_shapes : gpd.GeoDataFrame
            Lake polygons from HydroLAKES.
        tolerance : float, optional
            Simplification tolerance in CRS units (~100m for EPSG:3035),
            by default 0.001.

        Returns
        -------
        gpd.GeoDataFrame
            Lake shapes with simplified geometries.
        """
        from shapely.geometry import MultiPolygon, Polygon

        def remove_holes(geom):
            """Remove interior holes from polygon geometries."""
            if isinstance(geom, Polygon):
                return Polygon(geom.exterior)
            elif isinstance(geom, MultiPolygon):
                return MultiPolygon([Polygon(p.exterior) for p in geom.geoms])
            return geom

        lake_shapes = lake_shapes.copy()
        lake_shapes.geometry = lake_shapes.geometry.map(remove_holes).simplify(
            tolerance, preserve_topology=True
        )
        return lake_shapes

    @staticmethod
    def _approximate_lake_temperature(**kwargs) -> xr.DataArray:
        """
        Approximate lake water temperature from ambient air temperature.

        Uses the river temperature approximation from Triebs & Tsatsaronis 2022,
        which applies a moving average and empirical formula to derive water
        temperature from ambient air temperature.

        Parameters
        ----------
        **kwargs
            Keyword arguments passed to
            ``RiverWaterHeatApproximator._approximate_river_temperature``.

        Returns
        -------
        xr.DataArray
            Approximated lake water temperature.
        """
        return RiverWaterHeatApproximator._approximate_river_temperature(**kwargs)

    @cached_property
    def _volume_flow_in_region(self) -> float:
        """
        Calculate volume flow rate from total lake volume.

        Assumes annual turnover of lake volume, converted to hourly flow rate.
        Vol_total in HydroLAKES is in MCM (million cubic meters = 1e6 m³).

        Returns
        -------
        float
            Volume flow rate in m³/h.
        """
        lake_volume_mcm = self._lake_parts_in_region["Vol_total"].sum()
        lake_volume_m3 = lake_volume_mcm * 1e6
        volume_flow_m3_per_h = lake_volume_m3 / 8760
        return volume_flow_m3_per_h

    @cached_property
    def _lake_parts_in_region(self) -> gpd.GeoDataFrame:
        """
        Get lake polygons intersected with the region.

        Returns
        -------
        gpd.GeoDataFrame
            Lake parts that fall within the defined region.
        """
        if isinstance(self.region, gpd.GeoSeries):
            region_gdf = gpd.GeoDataFrame(geometry=self.region, crs=self.region.crs)
        else:
            region_gdf = gpd.GeoDataFrame(
                geometry=[self.region], crs=self.lake_shapes.crs
            )
        return gpd.overlay(self.lake_shapes, region_gdf, how="intersection")

    def _air_temperature_in_lakes(
        self, eligible_lake_parts: gpd.GeoDataFrame
    ) -> xr.DataArray:
        """
        Clip ambient temperature raster to lake areas.

        Parameters
        ----------
        eligible_lake_parts : gpd.GeoDataFrame
            Lake polygons to clip temperature data to.

        Returns
        -------
        xr.DataArray
            Ambient temperature clipped to lake areas.
        """
        return self.ambient_temperature.rio.clip(
            eligible_lake_parts.geometry, eligible_lake_parts.crs, drop=True
        )

    @cached_property
    def _water_temperature_in_region_raster(self) -> xr.DataArray:
        """
        Compute lake water temperature raster for the region.

        Returns
        -------
        xr.DataArray
            Water temperature raster with spatial (x, y) and time dimensions.
        """
        air_temperature_in_lakes = self._air_temperature_in_lakes(
            self._lake_parts_in_region
        )
        return self._approximate_lake_temperature(
            ambient_temperature=air_temperature_in_lakes
        )

    @cached_property
    def _water_temperature_in_region(self) -> xr.DataArray:
        """
        Compute spatially-averaged lake water temperature.

        Returns
        -------
        xr.DataArray
            Mean water temperature over time (spatial dimensions averaged).
        """
        return self._water_temperature_in_region_raster.mean(dim=("x", "y"))

    @cached_property
    def _power_sum_spatial(self) -> xr.DataArray:
        """
        Get total power (already spatially aggregated for lakes).

        Returns
        -------
        xr.DataArray
            Total thermal power over time.
        """
        return self._power_in_region

    def get_spatial_aggregate(self) -> xr.Dataset:
        """
        Get spatially aggregated heat potential results.

        Returns
        -------
        xr.Dataset
            Dataset with variables:
            - total_power : Total thermal power [MW] over time.
            - average_temperature : Mean water temperature [°C] over time.
        """
        return xr.Dataset(
            data_vars={
                "total_power": self._power_sum_spatial,
                "average_temperature": self._water_temperature_in_region,
            }
        )
