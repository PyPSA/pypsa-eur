# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class PTESOperatingApproximator:
    """
    Determines PTES usage for direct heating or as a heat source.

    Parameters
    ----------
    simplified_temperature : xr.DataArray
        Simplified temperature model of PTES.
    forward_temperature : xr.DataArray
        Forward temperature profile.
    delta_t_pinch_point : float
        Minimum temperature difference for PTES to act as a heat source.

    Returns
    -------
    direct_heating : xr.DataArray
        Boolean array indicating where PTES can directly distribute heat.
    heat_source : xr.DataArray
        Boolean array indicating where PTES can act as a heat source.
    """

    def __init__(
        self,
        forward_temperature_celsius: xr.DataArray,
        max_PTES_temperature: float,
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
        self.max_PTES_temperature = max_PTES_temperature

    def determine_ptes_usage(
        self
    ) -> xr.DataArray:
        """
        Determines PTES usage for direct heating or as a heat source.

        Parameters
        ----------
        simplified_temperature : xr.DataArray
            Simplified temperature model of PTES.
        forward_temperature : xr.DataArray
            Forward temperature profile.
        delta_t_pinch_point : float
            Minimum temperature difference for PTES to act as a heat source.

        Returns
        -------
        direct_heating : xr.DataArray
            Boolean array indicating where PTES can directly distribute heat.
        heat_source : xr.DataArray
            Boolean array indicating where PTES can act as a heat source.
        """
        # Direct heating: PTES simplified temperature > forward temperature
        operating_setting = self.forward_temperature > self.max_PTES_temperature

        return operating_setting

    # gibt mir debke ich nur boolean zurück momentan.

