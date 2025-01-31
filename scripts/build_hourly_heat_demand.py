# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build hourly heat demand time series from daily heat demand.

Water and space heating demand profiles are generated using intraday profiles from BDEW. Different profiles are used for the residential and services sectors as well as weekdays and weekend.

The daily heat demand is multiplied by the intraday profile to obtain the hourly heat demand time series. The rule is executed in ``build_sector.smk``.
"""

import logging
from itertools import product

import pandas as pd
import xarray as xr
from _helpers import (
    configure_logging,
    generate_periodic_profiles,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_hourly_heat_demand",
            scope="total",
            clusters=5,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    daily_space_heat_demand = (
        xr.open_dataarray(snakemake.input.heat_demand)
        .to_pandas()
        .reindex(index=snapshots, method="ffill")
    )

    intraday_profiles = pd.read_csv(snakemake.input.heat_profile, index_col=0)

    sectors = ["residential", "services"]
    uses = ["water", "space"]

    heat_demand = {}
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
        else:
            heat_demand[f"{sector} {use}"] = intraday_year_profile

    heat_demand = pd.concat(heat_demand, axis=1, names=["sector use", "node"])

    heat_demand.index.name = "snapshots"

    ds = heat_demand.stack(future_stack=True).to_xarray()

    ds.to_netcdf(snakemake.output.heat_demand)
