# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MITimport shapely
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
        min_inlet_temperature: float = 1,
    ):
        # buffer the region geometry by half the data resolution
        # This way, offshore data points just outside the region are included
        self.region_geometry = region_geometry.boundary.buffer(
            self._get_data_resolution(data=water_temperature)
            / self.METERS_PER_DEGREE
            / 1.5
        )

        self.water_temperature = water_temperature
        self.min_outlet_temperature = min_inlet_temperature
        self._mask_to_geometry()

    def _mask_to_geometry(self):
        """Mask water temperature to the geometry."""

        boxed_water_temperature = self._get_boxed_data(data=self.water_temperature)

        mask = self.get_geometry_mask(data=boxed_water_temperature)
        self.masked_water_temperature = boxed_water_temperature.where(mask)

    def get_spatial_aggregate(self):
        """Get the spatial aggregate of water temperature."""
        average_water_temperature = self.masked_water_temperature.mean(
            dim=[self.LATITUDE, self.LONGITUDE], skipna=True
        )

        # Combine into a single dataset and apply cut-off temperature
        return xr.Dataset(
            data_vars={
                "average_temperature": self._get_usable_water_temperature(
                    water_temperature=average_water_temperature
                ),
            }
        )

    def get_temporal_aggregate(self):
        """Get the temporal aggregate of water temperature."""
        average_water_temperature = self.masked_water_temperature.mean(
            dim=[self.TIME], skipna=True
        )

        # Combine into a single dataset
        # Don't apply cut-off temperature here because this is only used for plotting
        # and analysis
        return xr.Dataset(data_vars={"average_temperature": average_water_temperature})

    def _get_usable_water_temperature(
        self, water_temperature: xr.DataArray
    ) -> xr.DataArray:
        """Get the usable water temperature."""
        return xr.where(
            water_temperature > self.min_outlet_temperature,
            water_temperature,
            -1e9,  # absolute zero
        )
