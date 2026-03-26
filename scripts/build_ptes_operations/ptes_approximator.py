# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import numpy as np
import xarray as xr


class PtesApproximator:
    def __init__(
        self,
        forward_temperature: xr.DataArray,
        return_temperature: xr.DataArray,
        top_temperature: float,
        bottom_temperature: float,
        num_layers: int,
        design_top_temperature: float,
        design_bottom_temperature: float,
        design_standing_losses: float,
        interlayer_heat_transfer_coefficient: float,
    ):
        self.forward_temperature = forward_temperature
        self.return_temperature = return_temperature
        self.top_temperature = top_temperature
        self.bottom_temperature = bottom_temperature
        self.num_layers = num_layers
        self.design_top_temperature = design_top_temperature
        self.design_bottom_temperature = design_bottom_temperature
        self.design_standing_losses = design_standing_losses
        self.interlayer_heat_transfer_coefficient = interlayer_heat_transfer_coefficient

        self._time_dim = [d for d in forward_temperature.dims if d != "name"][0]

    @property
    def layer_temperatures(self) -> np.ndarray:
        """Calculate the fixed temperature of each layer [°C], hottest first."""
        return np.linspace(
            self.top_temperature, self.bottom_temperature, self.num_layers + 1
        )[:-1]

    @property
    def volume_weights(self) -> np.ndarray:
        """
        Volume weight W_l per layer.

        W_l = (T_top - T_bottom) / (T_l - T_bottom).
        A unit of energy in a colder layer occupies more volume.
        W_1 = 1 for the hottest layer; W_l ≥ 1 for colder layers.
        """
        weights = (self.design_top_temperature - self.design_bottom_temperature) / (
            self.layer_temperatures - self.bottom_temperature
        )
        return weights

    @property
    def charging_availability(self) -> xr.DataArray:
        """
        Binary charging availability Φ⁺_{l,t}.

        For each timestep, only the layer whose temperature is closest to the forward temperature can receive charge. Returns a DataArray with dimensions (snapshot, name, layer).
        """
        # |T_l - T_fwd_t| for each layer
        # Broadcast: T_fwd is (snapshot, name), T_layers is (layer,)
        layer_dim = xr.DataArray(
            self.layer_temperatures,
            dims=["layer"],
            coords={"layer": np.arange(len(self.layer_temperatures))},
        )
        abs_diff = np.abs(
            layer_dim - self.forward_temperature
        )  # (snapshot, name, layer)

        # Find the layer index with minimum distance at each (snapshot, name)
        closest_layer = abs_diff.argmin(dim="layer")  # (snapshot, name)

        # Build binary availability: 1 where layer == closest_layer
        layer_indices = xr.DataArray(
            np.arange(len(self.layer_temperatures)),
            dims=["layer"],
            coords={"layer": np.arange(len(self.layer_temperatures))},
        )

        binary_availability = (layer_indices == closest_layer).astype(float)
        binary_availability = binary_availability.transpose(
            self._time_dim, "name", "layer"
        )

        return binary_availability

    @property
    def interlayer_heat_transfer_coefficients(self) -> np.ndarray:
        """."""
        return np.full(self.num_layers - 1, self.interlayer_heat_transfer_coefficient)

    def _interlayer_heat_transfer_coefficients(
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
    def boost_ratios(self) -> xr.DataArray:
        """
        Resistive boost ratio α^RH_{l,t} per layer per timestep.

        α^RH_{l,t} = max(0, (T_fwd_t - T_l) / (T_l - T_bottom))

        Returns DataArray with dimensions (snapshot, name, layer).
        """
        layer_dim = xr.DataArray(
            self.layer_temperatures,
            dims=["layer"],
            coords={"layer": np.arange(len(self.layer_temperatures))},
        )

        numerator = self.forward_temperature - layer_dim  # (snapshot, name, layer)
        denominator = layer_dim - self.bottom_temperature  # (layer,)

        alpha = (numerator / denominator).clip(min=0)
        return alpha

    @property
    def e_max_pu(self) -> float:
        """
        Normalized storage capacity as fraction of design capacity.

        e_max_pu = (T_top - T_bottom) / (T_design_top - T_design_bottom),
        clipped to non-negative.
        """
        delta_t = self.top_temperature - self.bottom_temperature
        max_delta_t = self.design_top_temperature - self.design_bottom_temperature
        return max(0, delta_t / max_delta_t)

    @property
    def standing_losses(self) -> float:
        """Tbd."""
        return np.full(self.num_layers, self.design_standing_losses)

    def to_dataset(self) -> xr.Dataset:
        """
        Export all pre-computed parameters as a single xr.Dataset.

        Variables:
        - layer_temperatures: (layer,)
        - volume_weights: (layer,)
        - charging_availability: (snapshot, name, layer)
        - interlayer_transfer_coefficients: (layer_pair,)
        - boost_ratios: (snapshot, name, layer)
        - standing_losses: (layer,)
        """
        layer_coord = np.arange(self.num_layers)
        pair_coord = np.arange(self.num_layers - 1)

        ds = xr.Dataset(
            {
                "layer_temperatures": xr.DataArray(
                    self.layer_temperatures,
                    dims=["layer"],
                    coords={"layer": layer_coord},
                ),
                "volume_weights": xr.DataArray(
                    self.volume_weights,
                    dims=["layer"],
                    coords={"layer": layer_coord},
                ),
                "charging_availability": self.charging_availability,
                "interlayer_heat_transfer_coefficients": xr.DataArray(
                    self.interlayer_heat_transfer_coefficients,
                    dims=["layer_pair"],
                    coords={"layer_pair": pair_coord},
                ),
                "boost_ratios": self.boost_ratios,
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
            },
        )
        return ds
