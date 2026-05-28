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
        conservative_return_layer: bool = False,
        rho=1000,
        # Specific heat in MWh/kg/K to keep rho*c_p*DeltaT in MWh/m3.
        c_p=4184 / 3.6e9,
    ):
        self.forward_temperature = forward_temperature
        self.return_temperature = return_temperature
        # When True, the return-level and heat-pump reinjection layers are chosen
        # as the warmest layer at/below the target temperature (floor), so direct
        # and boosted discharge never deposit spent volume above the return level
        # -> no energy generation on discharge. When False, the closest layer is
        # used (and discharge_free_lunch_warnings flags any free-lunch risk).
        self.conservative_return_layer = bool(conservative_return_layer)

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
        # build_cop_profiles only expands to per-layer "ptes layer {i}" labels when
        # num_layers > 1; the single-layer build keeps the plain "ptes" source. Fall
        # back to it here instead of selecting a non-existent "ptes layer 0".
        def _heat_source(i: int) -> str:
            return f"ptes layer {i}" if self.num_layers > 1 else "ptes"

        return xr.concat(
            [
                self._cop.sel(
                    heat_system="urban central", heat_source=_heat_source(i)
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
        """
        Layer index used as the return level.

        With ``conservative_return_layer`` the warmest layer at or below the
        return temperature is chosen (floor), so direct discharge never deposits
        spent volume above the return level (no energy generation). If no layer is
        at/below the return temperature the coldest (bottom) layer is used.
        Otherwise the closest layer is chosen (nearest).
        """
        layer_t = self.layer_temperatures_da
        return_t = self.return_temperature
        if self._time_dim in getattr(return_t, "dims", ()):
            # Floor against the coldest return temperature to stay conservative
            # at every timestep.
            return_t = return_t.min(dim=self._time_dim)

        if self.conservative_return_layer:
            at_or_below = layer_t <= return_t
            floored = at_or_below.argmax(dim="layer")  # warmest layer with T <= T_ret
            return xr.where(
                at_or_below.any(dim="layer"), floored, self.num_layers - 1
            ).astype(int)

        return np.abs(layer_t - return_t).argmin(dim="layer").astype(int)

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

        layer_t = self.layer_temperatures_da
        layer_idx = layer_t.layer
        ret_val = []
        for l in range(self.num_layers):
            target = target_temp.sel(layer=l)
            colder = layer_idx > l  # candidates must be strictly colder than source
            if self.conservative_return_layer:
                # Warmest strictly-colder layer at/below the target depth (floor),
                # so the boosted volume is never deposited above its reference
                # temperature -> the heat pump can never create energy. Fall back
                # to the bottom layer if none is cold enough.
                at_or_below = (layer_t <= target) & colder
                chosen = xr.where(
                    at_or_below.any(dim="layer"),
                    at_or_below.argmax(dim="layer"),
                    self.num_layers - 1,
                )
            else:
                closest_layer = abs(layer_t - target).argmin(dim="layer")
                chosen = xr.where(
                    (closest_layer <= l) & (closest_layer < self.num_layers - 1),
                    closest_layer + 1,
                    closest_layer,
                )
            ret_val.append(chosen)

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

    def discharge_free_lunch_warnings(self, tol: float = 0.5) -> list[str]:
        """
        Flag nodes whose discretisation risks energy generation on discharge.

        Direct (non-boosted) discharge cools a layer from ``T_layer`` down to the
        district-heating return temperature and deposits the spent volume at the
        return-level layer, delivering ``rho*c_p*(T_layer - T_return)`` of heat
        while debiting the store by ``rho*c_p*(T_layer - T_return_layer)``. If the
        chosen (nearest) return-level layer sits *above* the actual return
        temperature (``T_return_layer > T_return``) the discharge creates
        ``rho*c_p*(T_return_layer - T_return)`` of heat per m3 — a free lunch.

        To stay conservative the layer temperatures should be chosen (in config)
        so the return-level layer is at or below the return temperature. This
        returns a human-readable warning per node where that does not hold (empty
        list if the discretisation is safe everywhere).
        """
        layer_t = self.layer_temperatures_da
        return_t = self.return_temperature
        # Per-node return temperature reduced to a scalar (mean over time if the
        # return profile is time-resolved).
        if self._time_dim in getattr(return_t, "dims", ()):
            return_t = return_t.mean(dim=self._time_dim)
        return_layer_t = layer_t.sel(layer=self.return_temp_layer)
        gap = (return_layer_t - return_t).compute()
        warnings = []
        for name in np.atleast_1d(gap["name"].values):
            g = float(gap.sel(name=name))
            if g > tol:
                warnings.append(
                    f"  {name}: return-level layer {float(return_layer_t.sel(name=name)):.1f} C "
                    f"> return T {float(return_t.sel(name=name)):.1f} C by {g:.1f} K"
                )
        return warnings

    def booster_cop(
        self,
        cop_approximator_cls,
        refrigerant,
        delta_t_pinch_point,
        isentropic_compressor_efficiency,
        heat_loss,
        min_delta_t_lift,
    ) -> xr.DataArray:
        """
        Per-layer booster COP recomputed at the actual evaporator outlet.

        The build_cop_profiles COP fixes the source outlet a shallow
        ``heat_source_cooling`` below the inlet. For the PTES booster the spent
        volume is in fact cooled all the way down to the reinjection
        (``hp_return``) layer, so the evaporator glide is ``T_source_inlet ->
        T_deposit``. Recomputing the COP with that outlet makes it reflect the
        real cooling depth (deeper cooling -> lower COP -> more electricity).

        Energy conservation is unaffected (heat = evaporator + electricity holds
        for any COP); this only sharpens the electricity/source split. Returns
        dimensions (layer, time, name).
        """
        layer_t = self.layer_temperatures_da
        ret = self.return_temperature
        hp_return = self.hp_return_layer  # (layer, name)
        cops = []
        for layer in range(self.num_layers):
            t_layer = float(self.layer_temperatures[layer])
            # preheater raises the return flow to T_layer when T_layer > T_ret, so
            # the heat pump source inlet is min(T_layer, T_ret) and its sink inlet
            # max(T_layer, T_ret).
            source_inlet = xr.where(ret < t_layer, ret, t_layer)
            sink_inlet = xr.where(ret < t_layer, t_layer, ret)
            source_outlet = layer_t.sel(layer=hp_return.sel(layer=layer)).drop_vars(
                "layer", errors="ignore"
            )
            # the evaporator can only cool (deposit at/below the inlet); clamp for
            # the non-boosted layers whose hp_return value is unused anyway.
            source_outlet = np.minimum(source_outlet, source_inlet)
            cop_layer = cop_approximator_cls(
                sink_outlet_temperature_celsius=self.forward_temperature,
                sink_inlet_temperature_celsius=sink_inlet,
                source_inlet_temperature_celsius=source_inlet,
                source_outlet_temperature_celsius=source_outlet,
                refrigerant=refrigerant,
                delta_t_pinch_point=delta_t_pinch_point,
                isentropic_compressor_efficiency=isentropic_compressor_efficiency,
                heat_loss=heat_loss,
                min_delta_t_lift=min_delta_t_lift,
            ).cop
            cops.append(cop_layer.drop_vars("layer", errors="ignore"))
        return xr.concat(cops, dim=pd.Index(np.arange(self.num_layers), name="layer"))

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
