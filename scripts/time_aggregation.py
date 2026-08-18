# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Defines the time aggregation to be used for sector-coupled network.

Description
-----------
Computes a time aggregation scheme for the given network, in the form of a CSV
file with the snapshot weightings, indexed by the new subset of snapshots. This
rule only computes said aggregation scheme; aggregation of time-varying network
data is done in `prepare_sector_network.py`.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import tsam
import xarray as xr

from scripts._helpers import (
    configure_logging,
    get_temporal_resolution,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def align_basis(data: pd.DataFrame, name: str, snapshots: pd.Index) -> pd.DataFrame:
    """Restrict a time series frame to the snapshots and flatten its columns."""
    data = data.loc[snapshots]
    data.columns = [f"{name}-{i}" for i in range(data.shape[1])]
    return data


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "time_aggregation",
            configfiles="config/test/config.overnight.yaml",
            horizon=2030,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    resolution = get_temporal_resolution(snakemake.params.time_resolution)

    # Representative snapshots and native resolution need no precomputed
    # weightings; they are handled directly in prepare_sector_network.py.
    if resolution is None or resolution[0] == "representative":
        logger.info("Use representative snapshot or no aggregation at all")
        pd.DataFrame().to_csv(snakemake.output.snapshot_weightings)

    # Plain resampling
    elif resolution[0] == "averaging":
        offset = f"{resolution[1]}h"
        logger.info(f"Averaging every {offset}")

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
    elif resolution[0] == "segmentation":
        segments = resolution[1]
        logger.info(f"Use temporal segmentation with {segments} segments")

        # The clustered network carries no time series yet, so the segmentation
        # basis is read from the resources that compose_network.py attaches later.
        sns = n.snapshots
        dfs = [
            align_basis(pnl, f"{c.name}-{attr}", sns)
            for c in n.components
            for attr, pnl in c.dynamic.items()
            if not pnl.empty and attr != "e_min_pu"
        ]

        for fn in snakemake.input.profiles:
            with xr.open_dataset(fn) as ds:
                if ds.indexes["bus"].empty:
                    continue
                if "year" in ds.indexes:
                    ds = ds.sel(year=ds.year.min(), drop=True)
                profile = ds.stack(bus_bin=["bus", "bin"])["profile"].to_pandas()
                dfs.append(align_basis(profile, Path(fn).stem, sns))

        if snakemake.input.hydro_profile:
            inflow = xr.open_dataarray(snakemake.input.hydro_profile).to_pandas()
            dfs.append(align_basis(inflow, "inflow", sns))

        demand = xr.open_dataarray(snakemake.input.electricity_demand).to_pandas()
        dfs.append(align_basis(demand, "electricity-demand", sns))

        if snakemake.input.hourly_heat_demand_total:
            heat_demand = (
                xr.open_dataset(snakemake.input.hourly_heat_demand_total)
                .to_dataframe()
                .unstack(level=1)
            )
            dfs.append(align_basis(heat_demand, "heat-demand", sns))
        if snakemake.input.solar_thermal_total:
            solar_thermal = (
                xr.open_dataset(snakemake.input.solar_thermal_total)
                .to_dataframe()
                .unstack(level=1)
            )
            dfs.append(align_basis(solar_thermal, "solar-thermal", sns))

        if not dfs:
            raise ValueError(
                "No time series available to derive the temporal segmentation from."
            )
        df = pd.concat(dfs, axis=1)

        # Reset columns to flat index
        df = df.T.reset_index(drop=True).T

        # Normalise all time-dependent data
        annual_max = df.max().replace(0, 1)
        df = df.div(annual_max, level=0)

        # Get representative segments
        if hasattr(tsam, "aggregate"):  # tsam >= 3.0
            agg = tsam.aggregate(
                df,
                n_clusters=1,
                period_duration=len(df),
                segments=tsam.SegmentConfig(n_segments=int(segments)),
            )
            agg = agg.cluster_representatives
        else:  # tsam < 3.0
            from tsam import timeseriesaggregation

            agg = timeseriesaggregation.TimeSeriesAggregation(
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
