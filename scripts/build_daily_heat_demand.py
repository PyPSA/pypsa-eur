# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
This rule builds heat demand time series using heating degree day (HDD)
approximation.

Snapshots are resampled to daily time resolution and ``Atlite.convert.heat_demand`` is used to convert ambient temperature from the default weather cutout to heat demand time series for the respective cutout.

Heat demand is distributed by population to clustered onshore regions.

.. seealso::
    `Atlite.Cutout.heat_demand <https://atlite.readthedocs.io/en/master/ref_api.html#module-atlite.convert>`_

"""

import logging

import atlite
import geopandas as gpd
import numpy as np
import xarray as xr
from _helpers import configure_logging, get_snapshots, set_scenario_config
from dask.distributed import Client, LocalCluster

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_daily_heat_demands",
            scope="total",
            clusters=48,
        )
    configure_logging(snakemake)
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
