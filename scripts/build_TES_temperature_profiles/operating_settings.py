# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr

def determine_ptes_usage(
    simplified_temperature: xr.DataArray,
    forward_temperature: xr.DataArray,
    delta_t_pinch_point: float,
) -> (xr.DataArray, xr.DataArray):
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
    direct_heating = simplified_temperature > forward_temperature

    # Heat source: Forward temperature exceeds PTES simplified temperature by at least delta_t_pinch_point
    heat_source = (forward_temperature - simplified_temperature) >= delta_t_pinch_point # schauen, wie das dann noch richtig integriere. Geht um das Delta T, ab wann Wärmepumpe läuft

    return direct_heating, heat_source
