# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


from typing import Union

import numpy as np
import xarray as xr
from BaseCopApproximator import BaseCopApproximator


class DecentralHeatingCopApproximator(BaseCopApproximator):
    """
    Approximate the coefficient of performance (COP) for a heat pump in a
    decentral heating system (individual/household heating).

    Uses a quadratic regression on the temperature difference between the source and sink based on empirical data proposed by Staffell et al. 2012.

    Attributes
    ----------
    forward_temperature_celsius : Union[xr.DataArray, np.array]
        The forward temperature in Celsius.
    source_inlet_temperature_celsius : Union[xr.DataArray, np.array]
        The source inlet temperature in Celsius.
    source_type : str
        The source of the heat pump. Must be either 'air' or 'ground'.

    Methods
    -------
    __init__(forward_temperature_celsius, source_inlet_temperature_celsius, source_type)
        Initialize the DecentralHeatingCopApproximator object.
    approximate_cop()
        Compute the COP values using quadratic regression for air-/ground-source heat pumps.
    _approximate_cop_air_source()
        Evaluate quadratic regression for an air-sourced heat pump.
    _approximate_cop_ground_source()
        Evaluate quadratic regression for a ground-sourced heat pump.

    References
    ----------
    [1] Staffell et al., Energy & Environmental Science 11 (2012): A review of domestic heat pumps, https://doi.org/10.1039/C2EE22653G.
    """

    def __init__(
        self,
        forward_temperature_celsius: Union[xr.DataArray, np.array],
        source_inlet_temperature_celsius: Union[xr.DataArray, np.array],
        source_type: str,
    ):
        """
        Initialize the DecentralHeatingCopApproximator object.

        Parameters
        ----------
        forward_temperature_celsius : Union[xr.DataArray, np.array]
            The forward temperature in Celsius.
        source_inlet_temperature_celsius : Union[xr.DataArray, np.array]
            The source inlet temperature in Celsius.
        source_type : str
            The source of the heat pump. Must be either 'air' or 'ground'.
        """

        self.delta_t = forward_temperature_celsius - source_inlet_temperature_celsius
        if source_type not in ["air", "ground"]:
            raise ValueError("'source_type' must be one of ['air', 'ground']")
        else:
            self.source_type = source_type

    def approximate_cop(self) -> Union[xr.DataArray, np.array]:
        """
        Compute the COP values using quadratic regression for air-/ground-
        source heat pumps.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The calculated COP values.
        """
        if self.source_type == "air":
            return self._approximate_cop_air_source()
        elif self.source_type == "ground":
            return self._approximate_cop_ground_source()

    def _approximate_cop_air_source(self) -> Union[xr.DataArray, np.array]:
        """
        Evaluate quadratic regression for an air-sourced heat pump.

        COP = 6.81 - 0.121 * delta_T + 0.000630 * delta_T^2

        Returns
        -------
        Union[xr.DataArray, np.array]
            The calculated COP values.
        """
        return 6.81 - 0.121 * self.delta_t + 0.000630 * self.delta_t**2

    def _approximate_cop_ground_source(self) -> Union[xr.DataArray, np.array]:
        """
        Evaluate quadratic regression for a ground-sourced heat pump.

        COP = 8.77 - 0.150 * delta_T + 0.000734 * delta_T^2

        Returns
        -------
        Union[xr.DataArray, np.array]
            The calculated COP values.
        """
        return 8.77 - 0.150 * self.delta_t + 0.000734 * self.delta_t**2
