# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from abc import ABC, abstractmethod
from typing import Union

import numpy as np
import xarray as xr


class BaseCopApproximator(ABC):
    """
    Abstract class for approximating the coefficient of performance (COP) of a
    heat pump.

    Attributes
    ----------
    forward_temperature_celsius : Union[xr.DataArray, np.array]
        The forward temperature in Celsius.
    source_inlet_temperature_celsius : Union[xr.DataArray, np.array]
        The source inlet temperature in Celsius.

    Methods
    -------
    __init__(self, forward_temperature_celsius, source_inlet_temperature_celsius)
        Initialize CopApproximator.
    approximate_cop(self)
        Approximate heat pump coefficient of performance (COP).
    celsius_to_kelvin(t_celsius)
        Convert temperature from Celsius to Kelvin.
    logarithmic_mean(t_hot, t_cold)
        Calculate the logarithmic mean temperature difference.
    """

    def __init__(
        self,
        forward_temperature_celsius: Union[xr.DataArray, np.array],
        source_inlet_temperature_celsius: Union[xr.DataArray, np.array],
    ):
        """
        Initialize CopApproximator.

        Parameters
        ----------
        forward_temperature_celsius : Union[xr.DataArray, np.array]
            The forward temperature in Celsius.
        source_inlet_temperature_celsius : Union[xr.DataArray, np.array]
            The source inlet temperature in Celsius.
        """
        pass

    @abstractmethod
    def approximate_cop(self) -> Union[xr.DataArray, np.array]:
        """
        Approximate heat pump coefficient of performance (COP).

        Returns
        -------
        Union[xr.DataArray, np.array]
            The calculated COP values.
        """
        pass

    @staticmethod
    def celsius_to_kelvin(
        t_celsius: Union[float, xr.DataArray, np.array],
    ) -> Union[float, xr.DataArray, np.array]:
        """
        Convert temperature from Celsius to Kelvin.

        Parameters
        ----------
        t_celsius : Union[float, xr.DataArray, np.array]
            Temperature in Celsius.

        Returns
        -------
        Union[float, xr.DataArray, np.array]
            Temperature in Kelvin.
        """
        if (np.asarray(t_celsius) > 200).any():
            raise ValueError(
                "t_celsius > 200. Are you sure you are using the right units?"
            )
        return t_celsius + 273.15

    @staticmethod
    def logarithmic_mean(
        t_hot: Union[float, xr.DataArray, np.ndarray],
        t_cold: Union[float, xr.DataArray, np.ndarray],
    ) -> Union[float, xr.DataArray, np.ndarray]:
        """
        Calculate the logarithmic mean temperature difference.

        Parameters
        ----------
        t_hot : Union[float, xr.DataArray, np.ndarray]
            Hot temperature.
        t_cold : Union[float, xr.DataArray, np.ndarray]
            Cold temperature.

        Returns
        -------
        Union[float, xr.DataArray, np.ndarray]
            Logarithmic mean temperature difference.
        """
        if (np.asarray(t_hot <= t_cold)).any():
            raise ValueError("t_hot must be greater than t_cold")
        return (t_hot - t_cold) / np.log(t_hot / t_cold)
