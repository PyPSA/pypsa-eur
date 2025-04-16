# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class BuildTESTemperature:
    """
    Approximates the top temperature profile for Pit Thermal Energy Storage (PTES).

    This class determines the PTES top temperature profile by evaluating the district
    heating forward temperature profile against a maximum operational temperature.
    If the forward temperature exceeds the threshold, the result is clipped to the
    maximum PTES temperature, reflecting the system's operational constraints.

    Parameters
    ----------
    forward_temperature_celsius : xr.DataArray
        The forward temperature profile (in Celsius) from the district heating network.
    max_PTES_temperature : float
        The maximum operational temperature (in Celsius) that PTES can directly utilize.

    Attributes
    ----------
    forward_temperature : xr.DataArray
        The input forward temperature profile from the district heating network.
    max_PTES_temperature : float
        The threshold temperature above which the PTES temperature is limited.
    """

    def __init__(self, forward_temperature_celsius: xr.DataArray, max_PTES_temperature: float):
        """
        Initialize BuildTESTemperature.

        Parameters
        ----------
        forward_temperature_celsius : xr.DataArray
            The forward temperature profile (in Celsius) from the district heating network.
        max_PTES_temperature : float
            The maximum operational temperature (in Celsius) that PTES can deliver.
        """
        self.forward_temperature = forward_temperature_celsius
        self.max_PTES_temperature = max_PTES_temperature

    @property
    def clipped_top_temperature(self) -> xr.DataArray:
        """
        Compute the clipped PTES top temperature profile.

        This property calculates the PTES top temperature by clipping the forward temperature
        profile at the predefined maximum PTES temperature. Values exceeding the maximum
        threshold are replaced by the maximum operational value, ensuring the profile adheres
        to system constraints.

        Returns
        -------
        xr.DataArray
            The resulting top temperature profile for PTES, with temperatures above the maximum
            operational limit clipped accordingly.
        """
        return self.forward_temperature.where(
            self.forward_temperature <= self.max_PTES_temperature,
            self.max_PTES_temperature
        )