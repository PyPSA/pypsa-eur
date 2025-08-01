# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr
from enum import Enum

class TesTemperatureProfile(Enum):
    """
    TES temperature profile assumptions.
    
    CONSTANT: Assumes fixed temperatures at operational limits.
        - Top temperature: constant at max_top_temperature
        - Bottom temperature: constant at min_bottom_temperature
        - Assumes charge-boosting to maintain top temperature.
        - NOTE: Assuming bottom_temperature = min_bottom_temperature ignores that cooling of the return temperature might be necessary in practice.
        
    DYNAMIC: Assumes temperatures follow network conditions.
        - Top temperature: follows forward_temperature (clipped at max_top_temperature)
        - Bottom temperature: follows return_temperature
        - Does not assume charge-boosting.
        - Note: This ignores that the TES temperatures do not match the supply temperatures due to thermal losses or other factors.
    """
    CONSTANT = "constant"
    DYNAMIC = "dynamic"


class PtesTemperatureApproximator:
    """
    A unified class to handle pit thermal energy storage (PTES) temperature-related calculations.

    It calculates top temperature profiles, determines when supplemental heating is needed,
    and approximates storage capacity based on temperature differences.

    Attributes
    ----------
    forward_temperature : xr.DataArray
        The forward temperature profile from the district heating network.
    return_temperature : xr.DataArray
        The return temperature profile from the district heating network.
    max_top_temperature : float
        Maximum operational temperature of top layer in PTES.
    min_bottom_temperature : float
        Minimum operational temperature of bottom layer in PTES.
    temperature_profile : TesTemperatureProfile
        TES temperature profile assumption.
    charge_boosting_required : bool
        Whether charge boosting is required/allowed.
    discharge_boosting_required : bool
        Whether discharge boosting is required/allowed.
    dynamic_capacity : bool
        Whether storage capacity varies with temperature. If False, assumes constant capacity.
    """

    def __init__(
        self,
        forward_temperature: xr.DataArray,
        return_temperature: xr.DataArray,
        max_top_temperature: float,
        min_bottom_temperature: float,
        temperature_profile: TesTemperatureProfile,
        charge_boosting_required: bool,
        discharge_boosting_required: bool,
        dynamic_capacity: bool,
    ):
        """
        Initialize PtesTemperatureApproximator.

        Parameters
        ----------
        forward_temperature : xr.DataArray
            The forward temperature profile from the district heating network.
        return_temperature : xr.DataArray
            The return temperature profile from the district heating network.
        max_top_temperature : float
            Maximum operational temperature of top layer in PTES.
        min_bottom_temperature : float
            Minimum operational temperature of bottom layer in PTES.
        temperature_profile : TesTemperatureProfile
            TES temperature profile assumption.
        charge_boosting_required : bool
            Whether charge boosting is required/allowed.
        discharge_boosting_required : bool
            Whether discharge boosting is required/allowed.
        dynamic_capacity : bool
            Whether storage capacity varies with temperature. If False, assumes constant capacity.
        """
        self.forward_temperature = forward_temperature
        self.return_temperature = return_temperature
        self.max_top_temperature = max_top_temperature
        self.min_bottom_temperature = min_bottom_temperature
        self.temperature_profile = temperature_profile
        self.charge_boosting_required = charge_boosting_required
        self.discharge_boosting_required = discharge_boosting_required
        self.dynamic_capacity = dynamic_capacity

    @property
    def top_temperature(self) -> xr.DataArray:
        """
        Forward temperature clipped at the maximum PTES temperature or constant max temperature.

        Returns
        -------
        xr.DataArray
            The resulting top temperature profile for PTES.
        """
        if self.temperature_profile == TesTemperatureProfile.CONSTANT:
            return xr.full_like(self.forward_temperature, self.max_top_temperature)
        elif self.temperature_profile == TesTemperatureProfile.DYNAMIC:
            return self.forward_temperature.where(
                self.forward_temperature <= self.max_top_temperature,
                self.max_top_temperature,
            )
        else:
            raise NotImplementedError(f"Temperature profile {self.temperature_profile} not implemented")

    @property
    def bottom_temperature(self) -> xr.DataArray:
        """
        Return temperature clipped at the minimum PTES temperature or constant min temperature.

        Returns
        -------
        xr.DataArray
            The resulting bottom temperature profile for PTES.
        """
        if self.temperature_profile == TesTemperatureProfile.CONSTANT:
            return xr.full_like(self.return_temperature, self.min_bottom_temperature)
        elif self.temperature_profile == TesTemperatureProfile.DYNAMIC:
            return self.return_temperature
        else:
            raise NotImplementedError(f"Temperature profile {self.temperature_profile} not implemented")

    @property
    def e_max_pu(self) -> xr.DataArray:
        """
        Calculate the normalized delta T for TES capacity in relation to
        max and min temperature.

        Returns
        -------
        xr.DataArray
            Normalized delta T values between 0 and 1, representing the
            available storage capacity as a percentage of maximum capacity.
            If dynamic_capacity is False, returns constant capacity of 1.0.
        """
        if self.dynamic_capacity:
            delta_t = self.top_temperature - self.bottom_temperature
            normalized_delta_t = delta_t / (
                self.max_top_temperature - self.min_bottom_temperature
            )
            return normalized_delta_t.clip(min=0)  # Ensure non-negative values
        else:
            return xr.ones_like(self.forward_temperature)

    @property
    def boost_per_discharge(self) -> xr.DataArray:
        """
        Calculate the additional lift required between the store's
        current top temperature and the forward temperature with the lift
        already achieved inside the store.

        Notes
        -----
        The total thermal output required to reach the forward temperature is:

            Q_total = Q_source + Q_boost

        The total heat transfer is partitioned into:

            Q_source = Ṽ·ρ·cₚ·(T_top − T_bottom)
            Q_boost  = Ṽ·ρ·cₚ·(T_forward − T_top)

        Solving this to constant Ṽ gives α as the ratio of required boost to available store energy:

            α = Q_boost / Q_source
              = (T_forward − T_top) / (T_top − T_bottom)

        This expression quantifies the share of PTES output that is covered
        by stored energy relative to the additional heating needed to meet
        the desired forward temperature.

        Returns
        -------
        xr.DataArray
            The resulting fraction of PTES charge that must be further heated.
        """
        if self.discharge_boosting_required:
            return ((self.forward_temperature - self.top_temperature) / (
                self.top_temperature - self.bottom_temperature
            )).where(self.forward_temperature > self.top_temperature, 0)
        else:
            return xr.zeros_like(self.forward_temperature)

    @property
    def boost_per_charge(self) -> xr.DataArray:
        """
        Calculate how much of the total energy needed to fill the PTES to its
        maximum capacity has already been delivered by charging up to the forward
        temperature, versus how much extra energy remains to reach the maximum.

        Notes
        -----
        To fill the storage from the return temperature all the way up to its
        maximum top temperature, the total thermal energy required is split into:

            Q_forward   = Ṽ·ρ·cₚ·(T_forward − T_bottom)
            Q_boosting  = Ṽ·ρ·cₚ·(T_top − T_forward)

        - Q_forward is the energy already delivered by charging to the forward setpoint.
        - Q_boosting is the extra boost energy still needed to reach maximum capacity.

        Defining α as the ratio of delivered energy to remaining boost energy:

            α = Q_boosting / Q_forward
              = (T_forward − T_bottom) /
                (T_top − T_forward)

        This ratio quantifies the share of the total charge process that has
        already been completed (via Q_forward) relative to what is still
        required (Q_boosting) to hit the maximum PTES top temperature.

        Wherever the forward temperature meets or exceeds the maximum, α is set
        to zero since no further boost is needed.

        Returns
        -------
        xr.DataArray
            The fraction of the PTES's available storage capacity already used.
        """
        if self.charge_boosting_required:
            return ((self.forward_temperature - self.return_temperature) / (
                self.max_top_temperature - self.forward_temperature
            )).where(self.forward_temperature < self.max_top_temperature, 0)
        else:
            return xr.zeros_like(self.forward_temperature)
