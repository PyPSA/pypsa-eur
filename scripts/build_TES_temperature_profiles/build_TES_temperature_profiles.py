# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import numpy as np
from typing import Union
import xarray as xr
import pandas as pd


class PTESTemperatureApproximator:
    """
    Class to model the temperature dynamics of a Pit Thermal Energy Storage (PTES).

    This class computes the temperature behavior in the top and bottom layers of the PTES,
    based on the forward and return temperature profiles of the district heating network.

    Attributes
    ----------
    forward_temperature_celsius : Union[xr.DataArray, np.array]
        Forward temperature profile in Celsius, representing the supply temperature of the district heating network.
    return_temperature_celsius : Union[xr.DataArray, np.array]
        Return temperature profile in Celsius, representing the temperature after heat distribution.
    max_PTES_temperature : float
        Maximum allowable storage temperature for PTES.
    snapshots : pd.DatetimeIndex
        Time snapshots used for temperature calculations.

    Methods
    -------
    simplified_top_layer_temperature_model() -> xr.DataArray
        Computes a sinusoidal temperature profile for the top layer of the PTES
        based on the maximum top temperature and minimum bottom temperature.

    ptes_top_temperature() -> xr.DataArray
        Determines the effective maximum temperature for the top layer of the PTES.
        The top temperature is the smaller value between the maximum PTES storage temperature
        and the highest forward temperature in the district heating network.

    ptes_bottom_temperature() -> xr.DataArray
        Retrieves the bottom temperature of the PTES storage, defined as the return
        temperature of the heating system.
    """

    def __init__(
        self,
        forward_temperature_celsius: Union[xr.DataArray, np.array],
        return_temperature_celsius: Union[xr.DataArray, np.array],
        max_PTES_temperature: float,
        snapshots: pd.DatetimeIndex,
    ):
        """
        Initialize the PTES Temperature Approximator.

        Parameters
        ----------
        forward_temperature_celsius : Union[xr.DataArray, np.array]
            Forward temperature profile in Celsius.
        return_temperature_celsius : Union[xr.DataArray, np.array]
            Return temperature profile in Celsius.
        max_PTES_temperature : float
            Maximum storage temperature for PTES.
        snapshots : pd.DatetimeIndex
            Time snapshots used for calculations.
        """
        self.forward_temperature = forward_temperature_celsius
        self.return_temperature = return_temperature_celsius
        self.max_PTES_temperature = max_PTES_temperature
        self.snapshots = snapshots  # Use snapshots instead of deriving time

    def simplified_top_layer_temperature_model(self) -> xr.DataArray:
        """
        Simplified top layer temperature model for PTES.

        Uses the effective maximum PTES top temperature (constant in time) and the
        minimum PTES bottom temperature to compute a sinusoidal base temperature profile.

        Returns
        -------
        xr.DataArray
            Top layer temperature profile for the given timestamps for each PTES,
            with dimensions ("time", "node") and time coordinates from self.snapshots.
        """
        # ptes_top_temperature now returns a DataArray with shape (time: 1, node: 6)
        # Remove the singleton "time" dimension to work with a (node: 6) array.
        max_ptes_top = self.ptes_top_temperature().squeeze("time")

        # Get the bottom temperature (assumed to be a DataArray with dim "node")
        min_ptes_bottom = self.ptes_bottom_temperature()

        # Compute the mean temperature and amplitude for each node.
        t_mean = (max_ptes_top + min_ptes_bottom) / 2
        t_amp = (max_ptes_top - min_ptes_bottom) / 2

        # Create a numeric time array.
        # Here we assume that self.snapshots (a DatetimeIndex) represents 8760 hourly time steps.
        t_hours = np.arange(len(self.snapshots))  # shape: (8760,)

        # Create an xarray DataArray for the time axis using self.snapshots as coordinates.
        time_da = xr.DataArray(t_hours, dims="time", coords={"time": self.snapshots})

        # Compute the sinusoidal term.
        # The sine function is applied on the time axis and is automatically broadcast
        # over the "node" dimension when combined with t_mean and t_amp.
        sin_term = np.sin(2 * np.pi * time_da / 8760 - 2.5)

        # Compute the base temperature by combining the time-constant node values with the time-varying sinusoidal term.
        base_temp = t_mean + t_amp * sin_term

        return base_temp

    def ptes_top_temperature(self) -> xr.DataArray:
        """
        Determines the effective maximum temperature for the top layer of each PTES.

        The maximum top layer temperature is the smaller value between:
          - `max_PTES_temperature` (global PTES limit)
          - The maximum forward temperature within the network (computed per PTES).

        Returns
        -------
        xr.DataArray
            Top layer temperature for each PTES with shape (time: 1, name: 6).
        """
        # Compute the maximum forward temperature for each PTES (across time).
        # This yields a DataArray with dims ("name",) and shape (6,).
        forward_max = self.forward_temperature.max(dim="time")

        # Determine the effective top temperature as the elementwise minimum between:
        # - the global limit, and
        # - the maximum forward temperature for each PTES.
        # Using xr.where preserves the DataArray structure and coordinates.
        effective_top = xr.where(forward_max < self.max_PTES_temperature,
                                 forward_max,
                                 self.max_PTES_temperature)

        # Expand the result to include a "time" dimension of length 1.
        effective_top = effective_top.expand_dims(dim="time", axis=0)

        # (Optional) If you want to check the dimensions, you can print them:

        return effective_top

    def ptes_bottom_temperature(self) -> xr.DataArray:
        """
        Determines the bottom temperature of the PTES storage.

        The bottom temperature is defined as the return temperature of the heating system.

        Returns
        -------
        xr.DataArray
            Bottom layer temperature for each PTES.
        """
        return self.return_temperature
