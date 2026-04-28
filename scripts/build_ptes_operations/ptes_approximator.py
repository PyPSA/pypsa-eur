# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import numpy as np
import pandas as pd
import xarray as xr


class PtesApproximator:
    def __init__(
        self,
        forward_temperature: xr.DataArray,
        return_temperature: xr.DataArray,
        layer_temperatures: list[float] | np.ndarray,
        design_top_temperature: float,
        design_bottom_temperature: float,
        design_standing_losses: float,
        interlayer_heat_transfer_coefficient: float,
        cop: xr.DataArray,
        rho=1000,
        # Specific heat in MWh/kg/K to keep rho*c_p*DeltaT in MWh/m3.
        c_p=4184 / 3.6e9,
    ):
        self.forward_temperature = forward_temperature
        self.return_temperature = return_temperature

        layer_temperatures = np.asarray(layer_temperatures, dtype=float)

        self._layer_temperatures = layer_temperatures
        self.top_temperature = float(layer_temperatures[0])
        self.bottom_temperature = float(layer_temperatures[-1])
        self.num_layers = int(layer_temperatures.size)
        self.design_top_temperature = design_top_temperature
        self.design_bottom_temperature = design_bottom_temperature
        self.design_standing_losses = design_standing_losses
        self.interlayer_heat_transfer_coefficient = interlayer_heat_transfer_coefficient
        self.rho = rho
        self.c_p = c_p
        self._cop = cop

        self._time_dim = [d for d in forward_temperature.dims if d != "name"][0]

    @property
    def cop(self):
        return xr.concat(
            [
                self._cop.sel(
                    heat_system="urban central", heat_source=f"ptes layer {i}"
                ).drop("heat_source")
                for i in range(self.num_layers)
            ],
            dim=pd.Index(np.arange(self.num_layers), name="layer"),
        ).drop("heat_system")

    @property
    def layer_temperatures(self) -> np.ndarray:
        """Fixed temperature of each layer [°C], hottest first."""
        return self._layer_temperatures

    @property
    def layer_temperatures_da(self) -> xr.DataArray:
        """Fixed temperature of each layer [°C], hottest first."""
        return xr.DataArray(
            self._layer_temperatures,
            dims=["layer"],
            coords={"layer": np.arange(len(self._layer_temperatures))},
        )

    @property
    def charging_availability(self) -> xr.DataArray:
        """
        Binary charging availability Φ⁺_{l,t}.

        For each timestep, only the layer whose temperature is closest to the forward temperature can receive charge. Returns a DataArray with dimensions (snapshot, name, layer).
        """

        closest_layer_to_forward = np.abs(
            self.forward_temperature - self.layer_temperatures_da
        ).argmin(dim="layer")  # (snapshot, name, layer)

        return xr.concat(
            [closest_layer_to_forward == l for l in range(self.num_layers)],
            dim=pd.Index(np.arange(self.num_layers), name="layer"),
        )

    @property
    def preheater_efficiency(self) -> xr.DataArray:
        """
        Conversion factor from volume to energy flow for heat exchange with return temperature. For each layer, this is (T_layer - T_return) * rho * c_p, but only if T_layer > T_return (otherwise 0). Returns a DataArray with dimensions (snapshot, name, layer).
        """

        return xr.where(
            self.return_temperature < self.layer_temperatures_da,
            self.rho
            * self.c_p
            * (self.layer_temperatures_da - self.return_temperature),
            0,
        )

    @property
    def heat_pump_efficiency(self) -> xr.DataArray:
        """
        Conversion factor from volume to energy flow for heat exchange with return temperature. For each layer, this is (T_layer - T_return) * rho * c_p, but only if T_layer > T_return (otherwise 0). Returns a DataArray with dimensions (snapshot, name, layer).
        """

        return xr.where(
            self.layer_temperatures_da > self.forward_temperature,
            0.0,
            xr.where(
                self.layer_temperatures_da > self.return_temperature,
                self.rho
                * self.c_p
                * (self.forward_temperature - self.layer_temperatures_da),
                self.rho
                * self.c_p
                * (self.forward_temperature - self.return_temperature),
            ),
        )

    @property
    def layer_needs_boosting(self) -> xr.DataArray:
        """
        Binary indicator whether a layer requires boosting to reach forward temperature.

        Returns dimensions (snapshot, name, layer).
        """
        return xr.where(self.forward_temperature > self.layer_temperatures_da, 1, 0)

    @property
    def charger_efficiency(self) -> xr.DataArray:
        """Heat-to-volume conversion efficiency for PTES charging [MWh to m3]."""
        return xr.where(
            self.charging_availability, self.charging_availability / self.m3_to_mwh, 0
        )

    @property
    def m3_to_mwh(self) -> np.ndarray:
        """Per-layer conversion from volume state to energy equivalent [MWh/m3]."""
        return (
            self.rho * self.c_p * (self.layer_temperatures_da - self.bottom_temperature)
        )

    @property
    def e_max_pu(self) -> float:
        """Scalar PTES capacity scaling against design temperature spread."""
        design_delta_t = self.design_top_temperature - self.design_bottom_temperature
        if design_delta_t <= 0:
            raise ValueError(
                "design_top_temperature must be larger than design_bottom_temperature"
            )

        operational_delta_t = self.top_temperature - self.bottom_temperature
        return max(operational_delta_t / design_delta_t, 0.0)

    @property
    def return_temp_layer(self) -> xr.DataArray:
        """Closest layer index to node-wise mean return temperature."""

        return (
            np.abs(self.layer_temperatures_da - self.return_temperature)
            .argmin(dim="layer")
            .astype(int)
        )

    @property
    def preheater_return_layer(self) -> xr.DataArray:
        """
        Reinjection layer index for the preheater's no-boost return volume.

        For each source layer l and node:
          a) T_l > T_ret and l is not the closest-to-T_ret layer: target = closest-to-T_ret
          b) T_l > T_ret and l is the closest-to-T_ret layer: target = min(l + 1, num_layers - 1)
             (next colder modelled layer; prevents the self-loop that produces phantom heat)
          c) T_l <= T_ret: no physical feedback; placeholder = num_layers - 1.
             The preheater's efficiency3 (= 1 - layer_needs_boosting) masks this flow
             to 0 because any layer colder than T_ret always needs boosting
             (T_fwd > T_ret >= T_l).

        Returns dimensions (layer, name).
        """
        ret_layer = self.return_temp_layer
        next_colder = xr.where(
            ret_layer < self.num_layers - 1, ret_layer + 1, self.num_layers - 1
        )
        layer_idx = self.layer_temperatures_da.layer

        return xr.where(
            layer_idx < ret_layer,
            ret_layer,
            xr.where(
                layer_idx == ret_layer,
                next_colder,
                self.num_layers - 1,
            ),
        ).astype(int)

    @property
    def hp_return_layer(self) -> xr.DataArray:
        """
        Reinjection layer index for PTES layer heat pumps.

        For each source layer l and node, select an "other" layer according to:
        a) if T_l < T_ret: target = T_l - dT
        b) else: target = T_ret - dT

        Then choose the candidate layer with minimum |T_other - target|,
        where candidates are strictly colder than the source layer.
        Returns dimensions (layer, name).
        """

        delta_t = xr.where(
            self.return_temperature < self.layer_temperatures_da,
            (self.forward_temperature - self.layer_temperatures_da)
            * (1 - 1 / self.cop.clip(min=1)),
            (self.forward_temperature - self.return_temperature)
            * (1 - 1 / self.cop.clip(min=1)),
        )

        target_temp = xr.where(
            self.return_temperature < self.layer_temperatures_da,
            self.return_temperature - delta_t,
            self.layer_temperatures_da - delta_t,
        )

        # take min over time (target layer must be constant)
        target_temp = target_temp.min(dim="time")

        ret_val = []
        for l in range(self.num_layers):
            closest_layer = abs(
                self.layer_temperatures_da - target_temp.sel(layer=l)
            ).argmin(dim="layer")
            closest_layer = xr.where(
                (closest_layer <= l) & (closest_layer < self.num_layers - 1),
                closest_layer + 1,
                closest_layer,
            )
            ret_val.append(closest_layer)

        return xr.concat(
            ret_val, dim=pd.Index(np.arange(self.num_layers), name="layer")
        )

    @property
    def interlayer_heat_transfer_coefficients(self) -> xr.DataArray:
        """."""
        return xr.full_like(
            self.layer_temperatures_da, self._interlayer_heat_transfer_coefficient()
        )

    def _interlayer_heat_transfer_coefficient(
        self,
        conductivity=0.6,
        density=1000,
        heat_capacity=4184,
        storage_height=15,
        thermocline=1,
        seconds_per_hour=3600,
    ) -> np.ndarray:
        """
        Interlayer heat transfer coefficient relative to layer state-of-charge. Requires "nominal" layer height, which is assumed to be storage_height distributed equally across layers. Returns array of length num_layers - 1, where each entry corresponds to the coefficient between layer l and l+1.

        Args:
        ----
        conductivity: Thermal conductivity of the storage medium [W/m/K]
        density: Density of the storage medium [kg/m3]
        heat_capacity: Specific heat capacity of the storage medium [J/kg/K]
        storage_height: Total height of the storage medium [m]
        seconds_per_hour: Number of seconds in an hour (for unit conversion)
        """
        return 0.001
        layer_height = storage_height / (
            self.num_layers + 1
        )  # assuming equal layer heights

        return (
            conductivity
            / (density * heat_capacity * layer_height)
            * seconds_per_hour
            * (self.layer_temperatures[:-1] - self.layer_temperatures[1:])
            / (self.layer_temperatures[:-1] - self.bottom_temperature)
            / thermocline
        )

    @property
    def standing_losses(self) -> float:
        """Tbd."""
        return np.full(self.num_layers, self.design_standing_losses)

    def to_dataset(self) -> xr.Dataset:
        """
        Export all pre-computed parameters as a single xr.Dataset.

        Variables:
        - layer_temperatures: (layer,)
        - return_temperature: (snapshot, name)
        - return_layer_index: (name,)
        - m3_to_mwh: (layer,)
        - hx_volume_to_return_temp: (snapshot, name, layer)
        - hx_volume_to_forward_temp: (snapshot, name, layer)
        - charger_efficiency: (snapshot, name, layer)
        - charging_availability: (snapshot, name, layer)
        - interlayer_transfer_coefficients: (layer_pair,)
        - standing_losses: (layer,)
        - e_max_pu: scalar
        """
        layer_coord = np.arange(self.num_layers)

        ds = xr.Dataset(
            {
                "layer_temperatures": xr.DataArray(
                    self.layer_temperatures,
                    dims=["layer"],
                    coords={"layer": layer_coord},
                ),
                "return_temperature_layer": self.return_temp_layer,
                "preheater_return_layer": self.preheater_return_layer,
                "hp_return_layer": self.hp_return_layer,
                "m3_to_mwh": xr.DataArray(
                    self.m3_to_mwh,
                    dims=["layer"],
                    coords={"layer": layer_coord},
                ),
                "preheater_efficiency": self.preheater_efficiency,
                "heat_pump_efficiency": self.heat_pump_efficiency,
                "layer_needs_boosting": self.layer_needs_boosting,
                "charger_efficiency": self.charger_efficiency,
                "charging_availability": self.charging_availability,
                "interlayer_heat_transfer_coefficients": self.interlayer_heat_transfer_coefficients,
                "standing_losses": xr.DataArray(
                    self.standing_losses,
                    dims=["layer"],
                    coords={"layer": layer_coord},
                ),
                "e_max_pu": self.e_max_pu,
            },
            attrs={
                "num_layers": self.num_layers,
                "top_temperature": self.top_temperature,
                "bottom_temperature": self.bottom_temperature,
                "storage_height": 15,  # m, assumed for interlayer coefficient calculation
                "design_top_temperature": self.design_top_temperature,
                "design_bottom_temperature": self.design_bottom_temperature,
                "design_standing_losses": self.design_standing_losses,
                "rho": self.rho,
                "c_p": self.c_p,
            },
        )
        return ds
