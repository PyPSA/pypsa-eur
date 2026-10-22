# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build hourly space and water heating demand profiles from daily heat demand.

The daily heat demand from [build_daily_heat_demand][] is multiplied by
intraday profiles from BDEW, which differ between residential and services
and between weekdays and weekends, to obtain hourly space heating demand per
node. Water heating uses the intraday profile alone since it does not depend
on temperature. A weekly availability profile for residential heat
demand-side management is also written; it is zero at the configured
`restriction_time` hours so that shifted heat demand must be served within
each period.
"""

import logging
from itertools import product

import numpy as np
import pandas as pd
import xarray as xr

from scripts._helpers import (
    configure_logging,
    generate_periodic_profiles,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def heat_dsm_profile(nodes, options):
    """
    Generate heat demand-side management (DSM) availability profile with periodic restrictions.

    Creates a weekly profile that restricts heat storage availability at configured
    checkpoint hours to enforce consumption requirements within 12-hour periods
    (day: 9am-9pm, night: 9pm-9am). This implements building thermal mass flexibility
    based on the smartEn/DNV methodology for residential heat DSM.

    The checkpoint approach operationally enforces the constraint that heat
    consumption requirements must be met within each perio (by default 12-hour periods,
    preventing the building thermal mass from acting as long-term seasonal storage while
    allowing short-term load shifting for demand-side flexibility.

    Parameters
    ----------
    nodes : pd.Index or array-like
        Node identifiers for which to generate profiles.
    options : dict
        Configuration dictionary containing:
        - ``options["residential_heat"]["dsm"]["restriction_time"]``: list of int
            Hours at which storage must be empty (checkpoint hours).

    Returns
    -------
    pd.DataFrame
        DSM availability profile indexed by timestamp with columns for each node.
        Values are 1.0 (storage available) for most hours and 0.0 at checkpoint
        hours to force storage depletion and enforce consumption periods.
    """
    weekly_profile = np.ones(24 * 7)
    for i in options["residential_heat"]["dsm"]["restriction_time"]:
        weekly_profile[(np.arange(0, 7, 1) * 24 + int(i))] = 0

    dsm_profile = generate_periodic_profiles(
        dt_index=pd.date_range(freq="h", **snakemake.params.snapshots, tz="UTC"),
        nodes=nodes,
        weekly_profile=weekly_profile,
    )

    return dsm_profile


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_hourly_heat_demand",
            scope="total",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    sector_options = snakemake.params.sector

    daily_space_heat_demand = (
        xr.open_dataarray(snakemake.input.heat_demand)
        .to_pandas()
        .reindex(index=snapshots, method="ffill")
    )

    intraday_profiles = pd.read_csv(snakemake.input.heat_profile, index_col=0)

    sectors = ["residential", "services"]
    uses = ["water", "space"]

    heat_demand = {}
    dsm_profile = {}
    for sector, use in product(sectors, uses):
        weekday = list(intraday_profiles[f"{sector} {use} weekday"])
        weekend = list(intraday_profiles[f"{sector} {use} weekend"])
        weekly_profile = weekday * 5 + weekend * 2
        intraday_year_profile = generate_periodic_profiles(
            daily_space_heat_demand.index.tz_localize("UTC"),
            nodes=daily_space_heat_demand.columns,
            weekly_profile=weekly_profile,
        )

        if use == "space":
            heat_demand[f"{sector} {use}"] = (
                daily_space_heat_demand * intraday_year_profile
            )
            if sector == "residential":
                dsm_profile[f"{sector} {use}"] = heat_dsm_profile(
                    daily_space_heat_demand.columns, sector_options
                )
        else:
            heat_demand[f"{sector} {use}"] = intraday_year_profile

    heat_demand = pd.concat(heat_demand, axis=1, names=["sector use", "node"])
    dsm_profile = pd.concat(dsm_profile, axis=1, names=["sector use", "node"])

    heat_demand.index.name = "snapshots"

    ds = heat_demand.stack(future_stack=True).to_xarray()

    ds.to_netcdf(snakemake.output.heat_demand)
    dsm_profile.to_csv(snakemake.output.heat_dsm_profile)
