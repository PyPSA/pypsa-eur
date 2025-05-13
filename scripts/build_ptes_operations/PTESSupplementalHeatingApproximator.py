# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class PTESSupplementalHeatingApproximator:
    """
    Approximates the need for supplemental heating in PTES integration.

    Parameters
    ----------
    forward_temperature : xr.DataArray
        Forward temperature profile from the district heating network.
    max_ptes_top_temperature : float
        The maximum operational PTES temperature.
    """

    def __init__(
        self, forward_temperature: xr.DataArray, max_ptes_top_temperature: float
    ):
        """
        Initialize the PTESSupplementalHeatingApproximator.

        Parameters
        ----------
        forward_temperature : xr.DataArray
            Forward temperature profile from the district heating network.
        max_ptes_top_temperature : float
            The maximum operational PTES temperature.
        """
        self.forward_temperature = forward_temperature
        self.max_ptes_top_temperature = max_ptes_top_temperature

    def determine_ptes_usage(self) -> xr.DataArray:
        """
        Identify timesteps requiring supplemental heating.

        Returns
        -------
        xr.DataArray
            Array with 1 for direct PTES usage, 0 if supplemental heating is needed.
        """
        return (self.forward_temperature < self.max_ptes_top_temperature).astype(int)
