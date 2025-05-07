# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class BuildTesTopTemperature:
    """
    Determines the top temperature profile for Thermal Energy Storage (TES).

    Attributes
    ----------
    forward_temperature_celsius : xr.DataArray
        The forward temperature profile (in Celsius) from the district heating network.
    max_PTES_temperature : float
        The maximum operational temperature (in Celsius) that Pit Thermal Energy Storage (PTES) can directly utilize.

    Methods
    -------
    clipped_top_temperature
        Forward temperature clipped at the maximum PTES temperature.
    """

    def __init__(
        self, forward_temperature_celsius: xr.DataArray, max_PTES_temperature: float
    ):
        """
        Initialize BuildTESTemperature.

        Parameters
        ----------
        forward_temperature_celsius : xr.DataArray
            The forward temperature profile from the district heating network.
        max_PTES_temperature : float
            The maximum operational PTES temperature .
        """
        self.forward_temperature = forward_temperature_celsius
        self.max_PTES_temperature = max_PTES_temperature

    @property
    def clipped_top_temperature(self) -> xr.DataArray:
        """
        Forward temperature clipped at the maximum PTES temperature.

        Returns
        -------
        xr.DataArray
            The resulting top temperature profile for PTES.
        """
        return self.forward_temperature.where(
            self.forward_temperature <= self.max_PTES_temperature,
            self.max_PTES_temperature,
        )
