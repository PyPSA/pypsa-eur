# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build solar thermal collector profile time series.

Uses ``atlite.Cutout.solar_thermal` to compute heat generation for clustered onshore regions from population layout and weather data cutout.
The rule is executed in ``build_sector.smk``.

.. seealso::
    `Atlite.Cutout.solar_thermal <https://atlite.readthedocs.io/en/master/ref_api.html#module-atlite.convert>`_
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

        snakemake = mock_snakemake("build_solar_thermal_profiles", clusters=48)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    cluster = LocalCluster(n_workers=nprocesses, threads_per_worker=1)
    client = Client(cluster, asynchronous=True)

    config = snakemake.params.solar_thermal
    config.pop("cutout", None)

    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)

    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=time)

    clustered_regions = (
        gpd.read_file(snakemake.input.regions_onshore).set_index("name").buffer(0)
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
        show_progress=False,
    )

    solar_thermal.to_netcdf(snakemake.output.solar_thermal)
