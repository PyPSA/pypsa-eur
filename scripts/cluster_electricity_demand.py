# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Aggregate simplified electricity demand to clustered resolution.

Maps the per-bus demand of ``electricity_demand_simplified.nc`` onto the
clustered network buses using ``busmap_cluster_network.csv``, so
``compose_network`` can attach the load directly without a clustering step.
"""

import logging

import pandas as pd
import xarray as xr

from scripts._helpers import (
    PYPSA_V1,
    configure_logging,
    sanitize_busmap,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("cluster_electricity_demand")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    load = xr.open_dataarray(snakemake.input.load)

    busmap = pd.read_csv(snakemake.input.busmap, dtype=str)
    index_col = "name" if PYPSA_V1 else "Bus"
    busmap = sanitize_busmap(busmap.set_index(index_col).squeeze())

    demand = load.to_dataframe().squeeze(axis=1).unstack(level="time")
    clustered = demand.groupby(busmap).sum().T
    clustered.index.name = "time"
    clustered.columns.name = "bus"

    out = xr.DataArray(clustered, name="electricity demand (MW)")
    comp = dict(zlib=True, complevel=9, least_significant_digit=5)
    out.to_netcdf(snakemake.output[0], encoding={out.name: comp})

    logger.info(f"Clustered electricity demand to {clustered.shape[1]} buses.")
