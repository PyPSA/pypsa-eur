# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build solar thermal collector time series.
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
            "build_solar_thermal_profiles",
            simpl="",
            clusters=48,
        )

    nprocesses = int(snakemake.threads)
    cluster = LocalCluster(n_workers=nprocesses, threads_per_worker=1)
    client = Client(cluster, asynchronous=True)

    config = snakemake.params.solar_thermal

    time = pd.date_range(freq="h", **snakemake.params.snapshots)
    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=time)

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

    nonzero_sum = M.sum(axis=0, keepdims=True)
    nonzero_sum[nonzero_sum == 0.0] = 1.0
    M_tilde = M / nonzero_sum

    solar_thermal = cutout.solar_thermal(
        **config,
        matrix=M_tilde.T,
        index=clustered_regions.index,
        dask_kwargs=dict(scheduler=client),
        show_progress=False
    )

    solar_thermal.to_netcdf(snakemake.output.solar_thermal)
