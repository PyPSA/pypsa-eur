# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build heat demand time series using heating degree day (HDD) approximation.
"""

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from dask.distributed import Client, LocalCluster

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_daily_heat_demands",
            scope="total",
            simpl="",
            clusters=48,
        )

    nprocesses = int(snakemake.threads)
    cluster = LocalCluster(n_workers=nprocesses, threads_per_worker=1)
    client = Client(cluster, asynchronous=True)

    time = pd.date_range(freq="h", **snakemake.params.snapshots)
    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=time)

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
    )

    heat_demand.to_netcdf(snakemake.output.heat_demand)
