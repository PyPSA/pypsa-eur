# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
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

    Uses a quadratic regression on the temperature difference between the source and sink based on empirical data proposed by Staffell et al. 2012 .

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
        Initialize the COPProfileBuilder object.

        Parameters:
        ----------
        forward_temperature_celsius : Union[xr.DataArray, np.array]
            The forward temperature in Celsius.
        return_temperature_celsius : Union[xr.DataArray, np.array]
            The return temperature in Celsius.
        source: str
            The source of the heat pump. Must be either 'air' or 'soil'
        """

        self.delta_t = forward_temperature_celsius - source_inlet_temperature_celsius
        if source_type not in ["air", "soil"]:
            raise ValueError("'source' must be one of ['air', 'soil']")
        else:
            self.source_type = source_type

    def approximate_cop(self) -> Union[xr.DataArray, np.array]:
        """
        Compute output of quadratic regression for air-/ground-source heat
        pumps.

        Calls the appropriate method depending on `source`.
        """
        if self.source_type == "air":
            return self._approximate_cop_air_source()
        elif self.source_type == "soil":
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
