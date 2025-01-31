# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Aggregate all rastered cutout data to base regions Voronoi cells.
"""

import logging

import atlite
import geopandas as gpd
from _helpers import configure_logging, get_snapshots, set_scenario_config
from atlite.aggregate import aggregate_matrix
from dask.distributed import Client

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_hac_features")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params
    nprocesses = int(snakemake.threads)

    if nprocesses > 1:
        client = Client(n_workers=nprocesses, threads_per_worker=1)
    else:
        client = None

    time = get_snapshots(params.snapshots, params.drop_leap_day)

    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=time)

    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    I = cutout.indicatormatrix(regions)  # noqa: E741

    ds = cutout.data[params.features].map(
        aggregate_matrix, matrix=I, index=regions.index
    )

    ds = ds.load(scheduler=client)

    ds.to_netcdf(snakemake.output[0])
