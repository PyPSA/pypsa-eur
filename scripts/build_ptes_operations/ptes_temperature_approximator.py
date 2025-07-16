# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


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
    max_ptes_top_temperature : float
        Maximum operational temperature of top layer in PTES, default 90°C.
    min_ptes_bottom_temperature : float
        Minimum operational temperature of bottom layer in PTES, default 35°C.
    """

    def __init__(
        self,
        forward_temperature: xr.DataArray,
        return_temperature: xr.DataArray,
        max_ptes_top_temperature: float = 90,
        min_ptes_bottom_temperature: float = 35,
    ):
        """
        Initialize PtesTemperatureApproximator.

        Parameters
        ----------
        forward_temperature : xr.DataArray
            The forward temperature profile from the district heating network.
        return_temperature : xr.DataArray
            The return temperature profile from the district heating network.
        max_ptes_top_temperature : float, optional
            Maximum operational temperature of top layer in PTES, default 90°C.
        min_ptes_bottom_temperature : float, optional
            Minimum operational temperature of bottom layer in PTES, default 35°C.
        """
        self.forward_temperature = forward_temperature
        self.return_temperature = return_temperature
        self.max_ptes_top_temperature = max_ptes_top_temperature
        self.min_ptes_bottom_temperature = min_ptes_bottom_temperature

    @property
    def top_temperature(self) -> xr.DataArray:
        """
        Forward temperature clipped at the maximum PTES temperature.

        Returns
        -------
        xr.DataArray
            The resulting top temperature profile for PTES.
        """
        return self.forward_temperature.where(
            self.forward_temperature <= self.max_ptes_top_temperature,
            self.max_ptes_top_temperature,
        )

    @property
    def bottom_temperature(self) -> xr.DataArray:
        """
        Return temperature clipped at the minimum PTES temperature.

        Returns
        -------
        xr.DataArray
            The resulting bottom temperature profile for PTES.
        """
        return self.min_ptes_bottom_temperature

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
        """
        delta_t = self.top_temperature - self.return_temperature
        normalized_delta_t = delta_t / (
            self.max_ptes_top_temperature - self.bottom_temperature
        )
        return normalized_delta_t.clip(min=0)  # Ensure non-negative values

    @property
    def temperature_boost_ratio(self) -> xr.DataArray:
        """
        Calculate the additional lift required between the store's
        current top temperature and the forward temperature with the lift
        already achieved inside the store.

        Notes
        -----
        The total thermal output required to reach the forward temperature is:

            Q_total = Q_source + Q_boost

        The total heat transfer is partitioned into:

            Q_source = Ṽ·ρ·cₚ·(T_max,store − T_return)
            Q_boost  = Ṽ·ρ·cₚ·(T_forward − T_max,store)

        Defining α as the ratio of required boost to available store energy:

            α = Q_boost / Q_source
              = (T_forward − T_max,store) / (T_max,store − T_return)

        This expression quantifies the share of PTES output that needs
        additional heating to meet the desired forward temperature.

        Returns
        -------
        xr.DataArray
            The resulting fraction of PTES charge that must be further heated.
        """
        return (self.forward_temperature - self.top_temperature) / (
            self.top_temperature - self.return_temperature
        )

    def forward_temperature_boost_ratio(self) -> xr.DataArray:
        """
        Calculate the additional lift required between the forward temperature and the maximum possible store's
        tmeperature to increase the maximum storage capacity.
        """
        return (self.forward_temperature - self.return_temperature) / (
            self.max_ptes_top_temperature - self.forward_temperature
        )