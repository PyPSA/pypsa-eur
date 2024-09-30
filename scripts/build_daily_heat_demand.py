# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This rule builds heat demand time series using heating degree day (HDD)
approximation.

Snapshots are resampled to daily time resolution and ``Atlite.convert.heat_demand`` is used to convert ambient temperature from the default weather cutout to heat demand time series for the respective cutout.

Heat demand is distributed by population to clustered onshore regions.

The rule is executed in ``build_sector.smk``.

.. seealso::
    `Atlite.Cutout.heat_demand <https://atlite.readthedocs.io/en/master/ref_api.html#module-atlite.convert>`_

Relevant Settings
-----------------

.. code:: yaml

    snapshots:
    drop_leap_day:

Inputs
------

- ``resources/<run_name>/pop_layout_<scope>.nc``: Population layout (spatial population distribution).
- ``resources/<run_name>/regions_onshore_base_s<simpl>_<clusters>.geojson``: Onshore region shapes.
- ``cutout``: Weather data cutout, as specified in config

Outputs
-------

- ``resources/daily_heat_demand_<scope>_base_s<simpl>_<clusters>.nc``:

Relevant settings
-----------------

.. code:: yaml

    atlite:
        default_cutout``:
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
            "build_daily_heat_demands",
            scope="total",
            clusters=48,
        )
    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    cluster = LocalCluster(n_workers=nprocesses, threads_per_worker=1)
    client = Client(cluster, asynchronous=True)

    cutout_name = snakemake.input.cutout

    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    daily = get_snapshots(
        snakemake.params.snapshots,
        snakemake.params.drop_leap_day,
        freq="D",
    )

    cutout = atlite.Cutout(cutout_name).sel(time=time)

    clustered_regions = (
        gpd.read_file(snakemake.input.regions_onshore).set_index("name").buffer(0)
    )

    I = cutout.indicatormatrix(clustered_regions)  # noqa: E741

    pop_layout = xr.open_dataarray(snakemake.input.pop_layout)

    stacked_pop = pop_layout.stack(spatial=("y", "x"))
    M = I.T.dot(np.diag(I.dot(stacked_pop)))

    heat_demand = cutout.heat_demand(
        matrix=M.T,
        index=clustered_regions.index,
        dask_kwargs=dict(scheduler=client),
        show_progress=False,
    ).sel(time=daily)

    heat_demand.to_netcdf(snakemake.output.heat_demand)
