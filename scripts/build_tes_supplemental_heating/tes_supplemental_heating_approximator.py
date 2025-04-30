# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class TESSupplementalHeatingApproximator:
    """
    Approximates the need for supplemental heating in TES integration.

    Parameters
    ----------
    forward_temperature_celsius : xr.DataArray
        Forward temperature profile from the district heating network.
    max_PTES_temperature : float
        The maximum operational PTES temperature.
    """

    def __init__(self, forward_temperature_celsius: xr.DataArray, max_PTES_temperature: float):
        """
        Initialize the TESSupplementalHeatingApproximator.

        Parameters
        ----------
        forward_temperature_celsius : xr.DataArray
            Forward temperature profile from the district heating network.
        max_PTES_temperature : float
            The maximum operational PTES temperature.
        """
        self.forward_temperature = forward_temperature_celsius
        self.max_PTES_temperature = max_PTES_temperature

    def determine_ptes_usage(self) -> xr.DataArray:
        """
        Identify timesteps requiring supplemental heating.


        Returns
        -------
        xr.DataArray
            Array with 1 for direct PTES usage, 0 if supplemental heating is needed.
        """
        return (self.forward_temperature < self.max_PTES_temperature).astype(int)
