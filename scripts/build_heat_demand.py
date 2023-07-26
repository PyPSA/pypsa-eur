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
            "build_heat_demands",
            weather_year="",
            simpl="",
            clusters=48,
        )

    nprocesses = int(snakemake.threads)
    cluster = LocalCluster(n_workers=nprocesses, threads_per_worker=1)
    client = Client(cluster, asynchronous=True)

    cutout_name = snakemake.input.cutout
    year = snakemake.wildcards.weather_year

    if year:
        snapshots = dict(start=year, end=str(int(year) + 1), inclusive="left")
        cutout_name = cutout_name.format(weather_year=year)
    else:
        snapshots = snakemake.params.snapshots

    drop_leap_day = snakemake.config["atlite"].get("drop_leap_day", False)
    time = pd.date_range(freq="h", **snapshots)
    daily = pd.date_range(freq="D", **snapshots)
    if drop_leap_day:
        time = time[~((time.month == 2) & (time.day == 29))]
        daily = daily[~((daily.month == 2) & (daily.day == 29))]

    cutout = atlite.Cutout(cutout_name).sel(time=time)

    clustered_regions = (
        gpd.read_file(snakemake.input.regions_onshore)
        .set_index("name")
        .buffer(0)
        .squeeze()
    )

    I = cutout.indicatormatrix(clustered_regions)

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
