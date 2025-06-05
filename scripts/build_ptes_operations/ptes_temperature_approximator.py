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
        Maximum operational temperature of top layer in PTES.
    min_ptes_bottom_temperature : float
        Minimum operational temperature of bottom layer in PTES.
    """

    def __init__(
        self,
        forward_temperature: xr.DataArray,
        return_temperature: xr.DataArray,
        max_ptes_top_temperature: float,
        min_ptes_bottom_temperature: float,
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
            Maximum operational temperature of top layer in PTES.
        min_ptes_bottom_temperature : float, optional
            Minimum operational temperature of bottom layer in PTES.
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
    def direct_utilisation_profile(self) -> xr.DataArray:
        """
        Identify timesteps requiring supplemental heating.

        Returns
        -------
        xr.DataArray
            Array with 1 for direct PTES usage, 0 if supplemental heating is needed.
        """
        return (self.forward_temperature <= self.max_ptes_top_temperature).astype(int)

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

        Returns
        -------
        xr.DataArray
            The resulting fraction of PTES charge that must be further heated.
        """
        return (self.forward_temperature - self.top_temperature) / (
            self.top_temperature - self.return_temperature
        )
