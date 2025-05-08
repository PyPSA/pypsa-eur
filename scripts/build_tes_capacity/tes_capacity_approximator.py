# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class TesCapacityApproximator:
    """
    A class to approximate thermal energy storage (TES) capacity based on temperature data.

    This class calculates the normalized temperature difference between the top and bottom
    layers of a pit thermal energy storage (PTES) system, which determines the available
    energy storage capacity. The top temperature is clipped to a maximum operational
    limit of 90°C.

    Attributes
    ----------
    top_temperature : xr.DataArray
        The top temperature data (forward flow temperature).
    bottom_temperature : xr.DataArray
        The bottom temperature data (return flow temperature).
    max_top_temperature : float
        Maximum operational temperature of top layer in PTES, default 90°C.
    min_bottom_temperature : float
        Minimum operational temperature of bottom layer in PTES, default 35°C.
    """

    def __init__(
        self,
        top_temperature: xr.DataArray,
        bottom_temperature: xr.DataArray,
        max_top_temperature: float = 90,
        min_bottom_temperature: float = 35,
    ):
        """
        Initialize TesCapacityApproximator.

        Parameters
        ----------
        top_temperature : xr.DataArray
            The top temperature data (forward flow temperature).
        bottom_temperature : xr.DataArray
            The bottom temperature data (return flow temperature).
        max_top_temperature : float, optional
            Maximum operational temperature of top layer in PTES, default 90°C.
        min_bottom_temperature : float, optional
            Minimum operational temperature of bottom layer in PTES, default 35°C.
        """
        self.top_temperature = top_temperature
        self.bottom_temperature = bottom_temperature
        self.max_top_temperature = max_top_temperature
        self.min_bottom_temperature = min_bottom_temperature

    @property
    def clipped_top_temperature(self) -> xr.DataArray:
        """
        Clip top temperature to maximum operational limit.

        Returns
        -------
        xr.DataArray
            Top temperature clipped to maximum operational limit.
        """
        return self.top_temperature.where(
            self.top_temperature <= self.max_top_temperature, self.max_top_temperature
        )

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
        delta_t = self.clipped_top_temperature - self.bottom_temperature
        normalized_delta_t = delta_t / (
            self.max_top_temperature - self.min_bottom_temperature
        )
        return normalized_delta_t.clip(min=0)  # Ensure non-negative values

    def calculate_e_max_pu(self, nodes=None) -> xr.DataArray:
        """
        Method to compute e_max_pu based on the top and bottom temperatures.

        Parameters
        ----------
        nodes : list or pd.Index, optional
            Subset of nodes to select from the temperature data.

        Returns
        -------
        xr.DataArray
            Computed e_max_pu values, optionally filtered for specific nodes.
        """
        result = self.e_max_pu
        if nodes is not None:
            result = result.sel(name=nodes)
        return result
