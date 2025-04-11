# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class PTESAdditionalHeatingApproximator:
    """
    Determines when additional heating is required for PTES integration.

    This class approximates whether the heat stored in thermal energy storage (TES) can be
    used directly for district heating or if further heating is required (e.g., via a heat pump).
    The decision is made by comparing the forward temperature profile of the district heating
    network with the maximum direct-usage temperature of the storage.

    Attributes
    ----------
    forward_temperature : xr.DataArray
        Temperature profile (in Celsius) of the district heating supply.
    max_PTES_temperature : float
        Maximum temperature (in Celsius) the PTES can deliver directly.

    Returns
    -------
    xr.DataArray
        DataArray with integer values:
            - 0 indicates that the TES can be used directly (no additional heating required).
            - 1 indicates that additional heating is needed.
    """

    def __init__(
            self,
            forward_temperature_celsius: xr.DataArray,
            max_PTES_temperature: float,
    ):
        """
        Initialize the PTESAdditionalHeatingApproximator.

        Parameters
        ----------
        forward_temperature_celsius : xr.DataArray
            The forward temperature profile (in Celsius) of the district heating network.
        max_PTES_temperature : float
            The maximum direct usage temperature (in Celsius) that the PTES can supply.
        """
        self.forward_temperature = forward_temperature_celsius
        self.max_PTES_temperature = max_PTES_temperature

    def determine_ptes_usage(self) -> xr.DataArray:
        """
        Determine when the PTES can directly supply heat versus when additional heating is required.

        This method compares the forward temperature profile with the maximum storage temperature.
        It outputs a DataArray with a value of 0 wherever the stored heat meets the demand directly
        (i.e., forward temperature is less than or equal to the maximum storage temperature) and
        a value of 1 wherever additional heating is needed (i.e., forward temperature exceeds the
        maximum storage temperature).

        Returns
        -------
        xr.DataArray
            DataArray of integers:
            - 0 indicates that direct usage is possible (no extra heating required).
            - 1 indicates that additional heating is necessary.
        """
        additional_heating_required = (self.forward_temperature > self.max_PTES_temperature).astype(int)
        return additional_heating_required
