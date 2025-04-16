# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class PTESAdditionalHeatingApproximator:
    """
    Determines when additional heating is required for PTES integration.

    This class evaluates whether the thermal energy stored in the PTES system
    can be directly used for district heating or if supplementary heating (for example, via
    a heat pump) is needed. The decision is made by comparing the district heating network's
    forward temperature profile with the maximum direct-usage temperature for PTES.

    Parameters
    ----------
    forward_temperature_celsius : xr.DataArray
        The forward temperature profile (in Celsius) from the district heating network.
    max_PTES_temperature : float
        The maximum operational temperature (in Celsius) that PTES can directly deliver.

    Attributes
    ----------
    forward_temperature : xr.DataArray
        The district heating forward temperature profile.
    max_PTES_temperature : float
        The maximum temperature threshold allowed for direct PTES usage.
    """

    def __init__(self, forward_temperature_celsius: xr.DataArray, max_PTES_temperature: float):
        """
        Initialize the PTESAdditionalHeatingApproximator.

        Parameters
        ----------
        forward_temperature_celsius : xr.DataArray
            The forward temperature profile (in Celsius) from the district heating network.
        max_PTES_temperature : float
            The maximum direct usage temperature (in Celsius) that the PTES can supply.
        """
        self.forward_temperature = forward_temperature_celsius
        self.max_PTES_temperature = max_PTES_temperature

    def determine_ptes_usage(self) -> xr.DataArray:
        """
        Compute when direct PTES usage is possible versus when additional heating is needed.

        This method assesses the forward temperature profile against the maximum allowed
        PTES temperature. It returns an xarray.DataArray containing integer values,
        where:
          - 1 indicates that the forward temperature is within or below the operational limit,
            meaning the PTES can be used directly.
          - 0 indicates that the forward temperature exceeds the limit and additional heating is required.

        Returns
        -------
        xr.DataArray
            An array of integer values (0 or 1) representing the need for extra heating.
        """
        additional_heating_required = (self.forward_temperature < self.max_PTES_temperature).astype(int)
        return additional_heating_required
