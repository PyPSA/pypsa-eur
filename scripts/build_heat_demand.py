# SPDX-FileCopyrightText: The PyPSA-Eur and -Earth authors
#
# SPDX-License-Identifier: MIT

"""

This rule does X.

Relevant Settings
-----------------

.. code:: yaml

    snapshots:

    load:
        interpolate_limit:
        time_shift_for_large_gaps:
        manual_adjustments:


.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`load_cf`

Inputs
------


Outputs
-------

- ``resource/XXX.csv``:


"""
import geopandas as gpd
import atlite
import pandas as pd
import xarray as xr
import numpy as np
import os


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, configure_logging
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            'build_heat_demand',
            simpl='',
            clusters=48,
        )
    configure_logging(snakemake)

    cutout = snakemake.input.cutout
    cutout = atlite.Cutout(cutout)

    grid_cells = cutout.grid_cells()

    # nuts3 has columns country, gdp, pop, geometry
    # population is given in dimensions of 1e3=k
    nuts3 = gpd.read_file(snakemake.input.nuts3_shapes).set_index('index')

    # Indicator matrix NUTS3 -> grid cells
    I = atlite.cutout.compute_indicatormatrix(nuts3.geometry, grid_cells)

    # Indicator matrix grid_cells -> NUTS3; inprinciple Iinv*I is identity
    # but imprecisions mean not perfect
    Iinv = cutout.indicatormatrix(nuts3.geometry)

    countries = np.sort(nuts3.country.unique())

    urban_fraction = pd.read_csv(snakemake.input.urban_percent,
                                header=None, index_col=0,
                                names=['fraction'], squeeze=True) / 100.

    # fill missing Balkans values
    missing = ["AL", "ME", "MK"]
    reference = ["RS", "BA"]
    average = urban_fraction[reference].mean()
    fill_values = pd.Series({ct: average for ct in missing})
    urban_fraction = urban_fraction.append(fill_values)

    # population in each grid cell
    pop_cells = pd.Series(I.dot(nuts3['pop']))

    ycoords = ('y', cutout.coords['y'].data)
    xcoords = ('x', cutout.coords['x'].data)
    values = pop_cells.values.reshape(cutout.shape)
    pop_layout = xr.DataArray(values, [ycoords, xcoords])


    regions_onshore = snakemake.input.regions_onshore
    regions = gpd.read_file(regions_onshore).set_index('name').buffer(0).squeeze()

    I = cutout.indicatormatrix(regions)

    stacked_pop = pop_layout.stack(spatial=('y', 'x'))
    M = I.T.dot(np.diag(I.dot(stacked_pop)))

    heat_demand = cutout.heat_demand(
        matrix=M.T, index=regions.index)

    heat_demand.to_netcdf(snakemake.output.heat_demand)
