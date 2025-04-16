# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build TES Temperature Profiles.

This module provides a class that approximates the top temperature profile for a
Pit Thermal Energy Storage (PTES) system based on the forward temperature profile
from the district heating network. The top temperature is defined as the forward
temperature if it does not exceed the maximum operational PTES temperature. If the
forward temperature is higher than this threshold, the output is clipped to the
maximum PTES temperature, reflecting the operational limits and indicating that any
excess temperature would require additional measures.

Relevant Settings
-----------------
.. code:: yaml
    storage:
        PTES:
            max_temperature: 90  # Maximum operational temperature for PTES

Inputs
------
- forward_temperature_celsius: xarray.DataArray
    Forward temperature profile from the district heating network (in Celsius).
- max_PTES_temperature: float
    Maximum operational temperature (in Celsius) that PTES can directly supply.

Outputs
-------
- clipped_top_temperature: xarray.DataArray
    The top temperature profile for PTES, clipped at the maximum operational limit.
"""

import xarray as xr


class BuildTESTemperature:
    """
    Approximates the top temperature profile for Pit Thermal Energy Storage (PTES).

    This class calculates the PTES top temperature by comparing the district heating
    forward temperature profile with the maximum operational temperature. When the
    forward temperature exceeds the set maximum, it is clipped to that limit.

    Parameters
    ----------
    forward_temperature_celsius : xr.DataArray
        The forward temperature profile (in Celsius) of the district heating network.
    max_PTES_temperature : float
        The maximum direct usage temperature (in Celsius) that the PTES can deliver.

    Attributes
    ----------
    forward_temperature : xr.DataArray
        The input forward temperature profile.
    max_PTES_temperature : float
        The maximum allowable temperature for direct PTES usage.
    """

    def __init__(self, forward_temperature_celsius: xr.DataArray, max_PTES_temperature: float):
        """
        Initialize the TES temperature approximator.

        Parameters
        ----------
        forward_temperature_celsius : xr.DataArray
            The forward temperature profile (in Celsius) of the district heating network.
        max_PTES_temperature : float
            The maximum direct usage temperature (in Celsius) that the PTES can deliver.
        """
        self.forward_temperature = forward_temperature_celsius
        self.max_PTES_temperature = max_PTES_temperature

    @property
    def clipped_top_temperature(self) -> xr.DataArray:
        """
        Calculate the PTES top temperature profile.

        The top temperature is computed by clipping the forward temperature to the maximum
        PTES operational temperature. Values above the maximum are replaced by the maximum,
        ensuring that the output respects the operational limits of the storage.

        Returns
        -------
        xr.DataArray
            Top temperature profile with values exceeding the maximum PTES temperature
            clipped to that maximum.
        """
        return self.forward_temperature.where(
            self.forward_temperature <= self.max_PTES_temperature,
            self.max_PTES_temperature
        )
