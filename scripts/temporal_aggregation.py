# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Applies the time aggregation, using snapshot weightings, on the sector-coupled network.

Description
-----------
Reads the snapshots weightings from the CSV file and applies it on the time-varying
network data prepared in ``prepare_sector_network.py``.
"""

import logging

import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)


def set_temporal_aggregation(
    n: pypsa.Network, resolution: str, snapshot_weightings_fn: str
):
    """
    Aggregate time-varying data to the given snapshots.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network with hourly resolution.
    resolution : str
        Temporal resolution specification.
    snapshot_weightings_fn : str
        Path to CSV file containing snapshot weightings for aggregation.

    Returns
    -------
    pypsa.Network
        Network with aggregated temporal resolution.
    """

    if not resolution:
        logger.info("No temporal aggregation. Using native resolution.")
        return n
    elif "sn" in resolution.lower():
        # Representative snapshots are dealt with directly
        sn = int(resolution[:-2])
        logger.info("Use every %s snapshot as representative", sn)
        n.set_snapshots(n.snapshots[::sn])
        n.snapshot_weightings *= sn
        return n
    else:
        # Otherwise, use the provided snapshots
        snapshot_weightings = pd.read_csv(
            snapshot_weightings_fn, index_col=0, parse_dates=True
        )

        # Define a series used for aggregation, mapping each hour in
        # n.snapshots to the closest previous timestep in
        # snapshot_weightings.index
        aggregation_map = (
            pd.Series(
                snapshot_weightings.index.get_indexer(n.snapshots), index=n.snapshots
            )
            .replace(-1, np.nan)
            .ffill()
            .astype(int)
            .map(lambda i: snapshot_weightings.index[i])
        )

        m = n.copy(snapshots=[])
        m.set_snapshots(snapshot_weightings.index)
        m.snapshot_weightings = snapshot_weightings

        # Aggregation all time-varying data.
        for c in n.iterate_components():
            pnl = getattr(m, c.list_name + "_t")
            for k, df in c.pnl.items():
                if not df.empty:
                    if c.list_name == "stores" and k == "e_max_pu":
                        pnl[k] = df.groupby(aggregation_map).min()
                    elif c.list_name == "stores" and k == "e_min_pu":
                        pnl[k] = df.groupby(aggregation_map).max()
                    else:
                        pnl[k] = df.groupby(aggregation_map).mean()

        return m


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "temporal_aggregation",
            configfiles="config/test/config.overnight.yaml",
            opts="",
            clusters="5",
            sector_opts="",
            planning_horizons="2030",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n_h = pypsa.Network(snakemake.input.network)

    n = set_temporal_aggregation(
        n=n_h,
        resolution=snakemake.params.time_resolution,
        snapshot_weightings_fn=snakemake.input.snapshot_weightings,
    )

    n.export_to_netcdf(snakemake.output[0])
