# -*- coding: utf-8 -*-
import logging

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "time_aggregation",
            configfiles="test/config.overnight.yaml",
            simpl="",
            opts="",
            clusters="37",
            ll="v1.0",
            sector_opts="CO2L0-24h-T-H-B-I-A-dist1",
            planning_horizons="2030",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n = pypsa.Network(snakemake.input.network)
    resolution = snakemake.params.time_resolution

    # Representative snapshots
    if "sn" in resolution.lower():
        logger.info("Use representative snapshots")
        # Output an empty csv; this is taken care of in prepare_sector_network.py
        pd.DataFrame().to_csv(snakemake.output.snapshot_weightings)

    # Plain resampling
    elif "h" in resolution.lower():
        offset = resolution.lower()
        logger.info(f"Averaging every {offset} hours")
        snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
        sns = snapshot_weightings.index
        if snakemake.params.drop_leap_day:
            sns = sns[~((sns.month == 2) & (sns.day == 29))]
        snapshot_weightings = snapshot_weightings.loc[sns]
        snapshot_weightings.to_csv(snakemake.output.snapshot_weightings)

    # Temporal segmentation
    elif "seg" in resolution.lower():
        segments = int(resolution[:-3])
        logger.info(f"Use temporal segmentation with {segments} segments")

        try:
            import tsam.timeseriesaggregation as tsam
        except ImportError:
            raise ModuleNotFoundError(
                "Optional dependency 'tsam' not found." "Install via 'pip install tsam'"
            )

        # Get all time-dependent data
        dfs = [
            pnl
            for c in n.iterate_components()
            for attr, pnl in c.pnl.items()
            if not pnl.empty and attr != "e_min_pu"
        ]
        dfs.append(
            xr.open_dataset(snakemake.input.hourly_heat_demand_total)
            .to_dataframe()
            .unstack(level=1)
        )
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
