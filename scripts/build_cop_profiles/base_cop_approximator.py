# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from abc import ABC, abstractmethod

import numpy as np
import xarray as xr


class BaseCopApproximator(ABC):
    """
    Abstract class for approximating the coefficient of performance (COP) of a
    heat pump.

    Attributes
    ----------
    sink_outlet_temperature_celsius : Union[xr.DataArray, np.array]
        The sink outlet temperature in Celsius.
    source_inlet_temperature_celsius : Union[xr.DataArray, np.array]
        The source inlet temperature in Celsius.

    Methods
    -------
    __init__(self, sink_outlet_temperature_celsius, source_inlet_temperature_celsius)
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
        sink_outlet_temperature_celsius: xr.DataArray | np.ndarray,
        source_inlet_temperature_celsius: xr.DataArray | np.ndarray,
    ):
        """
        Initialize CopApproximator.

        Parameters
        ----------
        sink_outlet_temperature_celsius : Union[xr.DataArray, np.array]
            The sink outlet temperature in Celsius.
        source_inlet_temperature_celsius : Union[xr.DataArray, np.array]
            The source inlet temperature in Celsius.
        """
        pass

    @property
    def cop(self) -> xr.DataArray | np.ndarray:
        """
        Calculate the coefficient of performance (COP) for the system.

        Returns
        -------
            Union[xr.DataArray, np.array]: The calculated COP values.
        """
        ret_val = self._approximate_cop()

        if isinstance(ret_val, xr.DataArray):
            # Use xarray's where method for DataArray
            ret_val = ret_val.where(ret_val >= 1, 0)
        else:
            # Use NumPy indexing for ndarray
            ret_val[ret_val < 1] = 0
        return ret_val

    @abstractmethod
    def _approximate_cop(self) -> xr.DataArray | np.ndarray:
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
        t_celsius: float | xr.DataArray | np.ndarray,
    ) -> float | xr.DataArray | np.ndarray:
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
        t_hot: float | xr.DataArray | np.ndarray,
        t_cold: float | xr.DataArray | np.ndarray,
    ) -> float | xr.DataArray | np.ndarray:
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
        if (np.asarray(t_hot < t_cold)).any():
            raise ValueError("t_hot must be greater than t_cold")
        return xr.where(
            t_hot == t_cold,
            t_hot,
            (t_hot - t_cold) / np.log(t_hot / t_cold),
        )
