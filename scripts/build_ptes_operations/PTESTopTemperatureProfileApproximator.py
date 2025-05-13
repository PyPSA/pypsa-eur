# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class PTESTopTemperatureProfileApproximator:
    """
    Determines the top temperature profile for Thermal Energy Storage (TES).

    Attributes
    ----------
    forward_temperature : xr.DataArray
        The forward temperature profile (in Celsius) from the district heating network.
    max_ptes_top_temperature : float
        The maximum operational temperature (in Celsius) that Pit Thermal Energy Storage (PTES) can directly utilize.

    Methods
    -------
    clipped_top_temperature
        Forward temperature clipped at the maximum PTES temperature.
    """

    def __init__(
        self, forward_temperature: xr.DataArray, max_ptes_top_temperature: float
    ):
        """
        Initialize PTESTopTemperatureProfileApproximator.

        Parameters
        ----------
        forward_temperature : xr.DataArray
            The forward temperature profile from the district heating network.
        max_ptes_top_temperature : float
            The maximum operational PTES temperature .
        """
        self.forward_temperature = forward_temperature
        self.max_ptes_top_temperature = max_ptes_top_temperature

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
            self.forward_temperature <= self.max_ptes_top_temperature,
            self.max_ptes_top_temperature,
        )
