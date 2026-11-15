# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Aggregates weather cutout time series to the regions of the simplified network as features for hierarchical agglomerative clustering (HAC).

Selected cutout variables, such as wind speed and solar influx, are aggregated
over each onshore region with the cutout's indicator matrix. The resulting
per-region time series are the feature vectors compared when the network is
clustered with the HAC algorithm.
"""

import logging

import geopandas as gpd
from atlite.aggregate import aggregate_matrix

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    load_cutout,
    set_scenario_config,
    setup_dask,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_hac_features")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params
    nprocesses = int(snakemake.threads)

    dask_kwargs = setup_dask(nprocesses)

    time = get_snapshots(params.snapshots, params.drop_leap_day)

    cutout = load_cutout(snakemake.input.cutout, time=time)

    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    I = cutout.indicatormatrix(regions)  # noqa: E741

    ds = cutout.data[params.features].map(
        aggregate_matrix, matrix=I, index=regions.index
    )

    ds = ds.load(**dask_kwargs)

    ds.to_netcdf(snakemake.output[0])
