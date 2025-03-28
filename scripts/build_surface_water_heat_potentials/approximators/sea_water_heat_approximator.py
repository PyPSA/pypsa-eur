import shapely
import xarray as xr

from scripts.build_surface_water_heat_potentials.approximators.surface_water_heat_approximator import (
    SurfaceWaterHeatApproximator,
)


class SeaWaterHeatApproximator(SurfaceWaterHeatApproximator):
    def __init__(
        self,
        water_temperature: xr.DataArray,
        region_geometry: shapely.geometry.polygon.Polygon,
        max_relative_volume_flow: float = 1.0,
        delta_t_max: float = 4,
        min_outlet_temperature: float = 1,
        min_distance_meters: int = 2000,
    ):
        # buffer the region geometry by half the data resolution
        # This way, offshore data points just outside the region are included
        self.region_geometry = region_geometry.buffer(
            self._get_data_resolution(data=water_temperature) / 2
        )

        self.water_temperature = water_temperature
        self.delta_t_max = delta_t_max
        self.min_outlet_temperature = min_outlet_temperature
        self.min_distance_meters = min_distance_meters

        self.results = self.get_results()

    def get_results(self):
        """Do some computations and return a result."""

        boxed_water_temperature = self._get_boxed_data(data=self.water_temperature)

        mask = self.get_geometry_mask(data=boxed_water_temperature)
        masked_water_temperature = boxed_water_temperature.where(mask)

        # Combine into a single dataset
        return xr.Dataset(
            data_vars={
                "average_temperature": masked_water_temperature.mean(
                    dim=["latitude", "longitude"], skipna=True
                ),
            }
        )
