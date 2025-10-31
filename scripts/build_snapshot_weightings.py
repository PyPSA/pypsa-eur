# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Defines the time aggregation, known as snapshot weightings, to be used for sector-coupled network.

Description
-----------
Computes a time aggregation scheme for the given network, in the form of a CSV
file with the snapshot weightings, indexed by the new subset of snapshots. This
rule only computes said aggregation scheme; aggregation of time-varying network
data is done in ``time_aggregation.py``.
"""

import logging

import numpy as np
import pandas as pd
import pypsa
import tsam.timeseriesaggregation as tsam
import xarray as xr

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_snapshot_weightings",
            configfiles="test/config.overnight.yaml",
            opts="",
            clusters="37",
            sector_opts="Co2L0-24h-T-H-B-I-A-dist1",
            planning_horizons="2030",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n = pypsa.Network(snakemake.input.network)
    resolution = snakemake.params.time_resolution

    if resolution["resolution_elec"] not in (False, 1, "1h", "1H"):
        raise ValueError(
            f"Invalid configuration: expected 'resolution_elec' = False for the "
            f"sector-coupled model, received {resolution['resolution_elec']!r}. "
            "Use 'resolution_sector' to define temporal resolution instead."
        )
    resolution = resolution["resolution_sector"]

    # Representative snapshots
    if not resolution or isinstance(resolution, str) and "sn" in resolution.lower():
        logger.info("Use representative snapshot or no aggregation at all")
        # Output an empty csv; this is taken care of in prepare_sector_network.py
        pd.DataFrame().to_csv(snakemake.output.snapshot_weightings)

    # Plain resampling
    elif isinstance(resolution, str) and "h" in resolution.lower():
        offset = resolution.lower()
        logger.info(f"Averaging every {offset} hours")

        # Resample years separately to handle non-contiguous years
        years = pd.DatetimeIndex(n.snapshots).year.unique()
        snapshot_weightings = []
        for year in years:
            sws_year = n.snapshot_weightings[n.snapshots.year == year]
            sws_year = sws_year.resample(offset).sum()
            snapshot_weightings.append(sws_year)
        snapshot_weightings = pd.concat(snapshot_weightings)

        # The resampling produces a contiguous date range. In case the original
        # index was not contiguous, all rows with zero weight must be dropped
        # (corresponding to time steps not included in the original snapshots).
        zeros_i = snapshot_weightings.query("objective == 0").index
        snapshot_weightings.drop(zeros_i, inplace=True)

        swi = snapshot_weightings.index
        leap_days = swi[(swi.month == 2) & (swi.day == 29)]
        if snakemake.params.drop_leap_day and not leap_days.empty:
            for year in leap_days.year.unique():
                year_leap_days = leap_days[leap_days.year == year]
                leap_weights = snapshot_weightings.loc[year_leap_days].sum()
                march_first = pd.Timestamp(year, 3, 1, 0, 0, 0)
                snapshot_weightings.loc[march_first] = leap_weights
            snapshot_weightings = snapshot_weightings.drop(leap_days).sort_index()

        sns = snapshot_weightings.index
        snapshot_weightings = snapshot_weightings.loc[sns]
        snapshot_weightings.to_csv(snakemake.output.snapshot_weightings)

    # Temporal segmentation
    elif isinstance(resolution, str) and "seg" in resolution.lower():
        segments = int(resolution[:-3])
        logger.info(f"Use temporal segmentation with {segments} segments")

        # Get all time-dependent data
        dfs = [
            pnl
            for c in n.iterate_components()
            for attr, pnl in c.pnl.items()
            if not pnl.empty and attr != "e_min_pu"
        ]
        if snakemake.input.hourly_heat_demand_total:
            dfs.append(
                xr.open_dataset(snakemake.input.hourly_heat_demand_total)
                .to_dataframe()
                .unstack(level=1)
            )
        if snakemake.input.solar_thermal_total:
            dfs.append(
                xr.open_dataset(snakemake.input.solar_thermal_total)
                .to_dataframe()
                .unstack(level=1)
            )
        df = pd.concat(dfs, axis=1)

        # Reset columns to flat index
        df = df.T.reset_index(drop=True).T

        # Normalise all time-dependent data
        annual_max = df.max().replace(0, 1)
        df = df.div(annual_max, level=0)

        # Get representative segments
        agg = tsam.TimeSeriesAggregation(
            df,
            hoursPerPeriod=len(df),
            noTypicalPeriods=1,
            noSegments=segments,
            segmentation=True,
            solver=snakemake.params.solver_name,
        )
        agg = agg.createTypicalPeriods()

        weightings = agg.index.get_level_values("Segment Duration")
        offsets = np.insert(np.cumsum(weightings[:-1]), 0, 0)
        snapshot_weightings = n.snapshot_weightings.loc[n.snapshots[offsets]].mul(
            weightings, axis=0
        )

        logger.info(
            f"Distribution of snapshot durations:\n{snapshot_weightings.objective.value_counts()}"
        )

        snapshot_weightings.to_csv(snakemake.output.snapshot_weightings)
