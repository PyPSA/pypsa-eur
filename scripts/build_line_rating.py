# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Adds dynamic line rating timeseries to the base network.

Relevant Settings
-----------------

.. code:: yaml

    lines_t:
        s_max_pu


.. seealso::
    Documentation of the configuration file ``config.yaml`
Inputs
------

- ``data/cutouts``: 
- ``networks/base.nc``: confer :ref:`base`

Outputs
-------

- ``networks/base_with_line_rating.nc``


Description
-----------

The rule :mod:`build_line_rating` calculates the line rating for transmission lines. 
The line rating provides the maximal capacity of a transmission line considering the heat exchange with the environment. 

The folloing heat gains and losses are considered:

- heat gain through resistive losses
- heat gain trough solar radiation
- heat loss through radiation of the trasnmission line
- heat loss through forced convection with wind
- heat loss through natural convection 


With a heat balance considering the maximum temperature threshold of the tranmission line, 
the maximal possible capacity factor "s_max_pu" for each transmission line at each time step is calculated.
"""

import logging
from _helpers import configure_logging

import pypsa
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, LineString as Line
import atlite
import xarray as xr


def add_line_rating(n):
    buses = n.lines[["bus0", "bus1"]].values
    x = n.buses.x
    y = n.buses.y
    shapes = [Line([Point(x[b0], y[b0]), Point(x[b1], y[b1])]) for (b0, b1) in buses]
    shapes = gpd.GeoSeries(shapes, index=n.lines.index)
    cutout = atlite.Cutout(snakemake.input.cutout)
    da = xr.DataArray(data=np.sqrt(3) * cutout.line_rating(shapes, n.lines.r/n.lines.length) * 1e3,
                      attrs=dict(description="Maximal possible power for given line considering line rating")) # in MW
    return da
    #import netcdf file in add electricity.py
    #n.lines_t.s_max_pu=s.to_pandas().transpose()/n.lines.s_nom
    #n.lines_t.s_max_pu.replace(np.inf, 1.0, inplace=True)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_line_rating', network='elec', simpl='',
                                  clusters='5', ll='copt', opts='Co2L-BAU-CCL-24H')
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)    
    da=add_line_rating(n)

    da.to_netcdf(snakemake.output[0])