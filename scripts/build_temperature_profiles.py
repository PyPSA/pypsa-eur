# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build time series for air and soil temperatures per clustered model region.

Uses ``atlite.Cutout.temperature`` and ``atlite.Cutout.soil_temperature compute temperature ambient air and soil temperature for the respective cutout. The rule is executed in ``build_sector.smk``.


.. seealso::
    `Atlite.Cutout.temperature <https://atlite.readthedocs.io/en/master/ref_api.html#module-atlite.convert>`_
    `Atlite.Cutout.soil_temperature <https://atlite.readthedocs.io/en/master/ref_api.html#module-atlite.convert>`_

Relevant Settings
-----------------

.. code:: yaml

    snapshots:
    drop_leap_day:
    atlite:
        default_cutout:

Inputs
------

- ``resources/<run_name>/pop_layout_total.nc``:
- ``resources/<run_name>/regions_onshore_base_s<simpl>_<clusters>.geojson``:
- ``cutout``: Weather data cutout, as specified in config

Outputs
-------

- ``resources/temp_soil_total_base_s<simpl>_<clusters>.nc``:
- ``resources/temp_air_total_base_s<simpl>_<clusters>.nc`
"""

import atlite
import geopandas as gpd
import numpy as np
import xarray as xr
from _helpers import get_snapshots, set_scenario_config
from dask.distributed import Client, LocalCluster

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_temperature_profiles",
            clusters=48,
        )
    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    cluster = LocalCluster(n_workers=nprocesses, threads_per_worker=1)
    client = Client(cluster, asynchronous=True)

    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)

    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=time)

    clustered_regions = (
        gpd.read_file(snakemake.input.regions_onshore).set_index("name").buffer(0)
    )

    I = cutout.indicatormatrix(clustered_regions)  # noqa: E741

    pop_layout = xr.open_dataarray(snakemake.input.pop_layout)

    stacked_pop = pop_layout.stack(spatial=("y", "x"))
    M = I.T.dot(np.diag(I.dot(stacked_pop)))

    nonzero_sum = M.sum(axis=0, keepdims=True)
    nonzero_sum[nonzero_sum == 0.0] = 1.0
    M_tilde = M / nonzero_sum

    temp_air = cutout.temperature(
        matrix=M_tilde.T,
        index=clustered_regions.index,
        dask_kwargs=dict(scheduler=client),
        show_progress=False,
    )

    temp_air.to_netcdf(snakemake.output.temp_air)

    temp_soil = cutout.soil_temperature(
        matrix=M_tilde.T,
        index=clustered_regions.index,
        dask_kwargs=dict(scheduler=client),
        show_progress=False,
    )

    temp_soil.to_netcdf(snakemake.output.temp_soil)
