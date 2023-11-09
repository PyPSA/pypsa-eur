# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Adds dynamic line rating timeseries to the base network.

Relevant Settings
-----------------

.. code:: yaml

    lines:
        cutout:
        line_rating:


.. seealso::
    Documentation of the configuration file ``config.yaml`
Inputs
------

- ``data/cutouts``:
- ``networks/base.nc``: confer :ref:`base`

Outputs
-------

- ``resources/line_rating.nc``


Description
-----------

The rule :mod:`build_line_rating` calculates the line rating for transmission lines.
The line rating provides the maximal capacity of a transmission line considering the heat exchange with the environment.

The following heat gains and losses are considered:

- heat gain through resistive losses
- heat gain through solar radiation
- heat loss through radiation of the transmission line
- heat loss through forced convection with wind
- heat loss through natural convection


With a heat balance considering the maximum temperature threshold of the transmission line,
the maximal possible capacity factor "s_max_pu" for each transmission line at each time step is calculated.
"""

import logging
import re

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import configure_logging
from shapely.geometry import LineString as Line
from shapely.geometry import Point


def calculate_resistance(T, R_ref, T_ref=293, alpha=0.00403):
    """
    Calculates the resistance at other temperatures than the reference
    temperature.

    Parameters
    ----------
    T : Temperature at which resistance is calculated in [°C] or [K]
    R_ref : Resistance at reference temperature in [Ohm] or [Ohm/Per Length Unit]
    T_ref : Reference temperature in [°C] or [K]
    alpha: Temperature coefficient in [1/K]
        Defaults are:
            * T_ref : 20 °C
            * alpha : 0.00403 1/K

    Returns
    -------
    Resistance of at given temperature.
    """
    return R_ref * (1 + alpha * (T - T_ref))


def calculate_line_rating(n, cutout):
    """
    Calculates the maximal allowed power flow in each line for each time step
    considering the maximal temperature.

    Parameters
    ----------
    n : pypsa.Network object containing information on grid

    Returns
    -------
    xarray DataArray object with maximal power.
    """
    relevant_lines = n.lines[(n.lines["underground"] == False)]
    buses = relevant_lines[["bus0", "bus1"]].values
    x = n.buses.x
    y = n.buses.y
    shapes = [Line([Point(x[b0], y[b0]), Point(x[b1], y[b1])]) for (b0, b1) in buses]
    shapes = gpd.GeoSeries(shapes, index=relevant_lines.index)
    if relevant_lines.r_pu.eq(0).all():
        # Overwrite standard line resistance with line resistance obtained from line type
        r_per_length = n.line_types["r_per_length"]
        R = (
            relevant_lines.join(r_per_length, on=["type"])["r_per_length"] / 1000
        )  # in meters
        # If line type with bundles is given retrieve number of conductors per bundle
        relevant_lines["n_bundle"] = (
            relevant_lines["type"]
            .where(relevant_lines["type"].str.contains("bundle"))
            .dropna()
            .apply(lambda x: int(re.findall(r"(\d+)-bundle", x)[0]))
        )
        # Set default number of bundles per line
        relevant_lines["n_bundle"].fillna(1, inplace=True)
        R *= relevant_lines["n_bundle"]
        R = calculate_resistance(T=353, R_ref=R)
    Imax = cutout.line_rating(shapes, R, D=0.0218, Ts=353, epsilon=0.8, alpha=0.8)
    line_factor = relevant_lines.eval("v_nom * n_bundle * num_parallel") / 1e3  # in mW
    return xr.DataArray(
        data=np.sqrt(3) * Imax * line_factor.values.reshape(-1, 1),
        attrs=dict(
            description="Maximal possible power in MW for given line considering line rating"
        ),
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_line_rating",
            network="elec",
            simpl="",
            clusters="5",
            ll="v1.0",
            opts="Co2L-4H",
        )
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)
    time = pd.date_range(freq="h", **snakemake.config["snapshots"])
    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=time)

    da = calculate_line_rating(n, cutout)
    da.to_netcdf(snakemake.output[0])
