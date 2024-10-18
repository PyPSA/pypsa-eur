# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build hourly heat demand time series from daily heat demand.

Water and space heating demand profiles are generated using intraday profiles from BDEW. Different profiles are used for the residential and services sectors as well as weekdays and weekend.

The daily heat demand is multiplied by the intraday profile to obtain the hourly heat demand time series. The rule is executed in ``build_sector.smk``.


Relevant Settings
-----------------

.. code:: yaml

    snapshots:
    drop_leap_day:

Inputs
------

- ``data/heat_load_profile_BDEW.csv``: Intraday heat profile for water and space heating demand for the residential and services sectors for weekends and weekdays.
- ``resources/daily_heat_demand_total_base_s<simpl>_<clusters>.nc``: Daily heat demand per cluster.

Outputs
-------

- ``resources/hourly_heat_demand_total_base_s<simpl>_<clusters>.nc``:
"""

from itertools import product

import numpy as np
import pandas as pd
import xarray as xr
from _helpers import generate_periodic_profiles, get_snapshots, set_scenario_config

def heat_dsm_profile(nodes, options):

    weekly_profile = np.ones((24 * 7))
    for i in options["residential_heat_restriction_time"]:
        weekly_profile[(np.arange(0, 7, 1) * 24 + int(i))] = 0

    dsm_profile = generate_periodic_profiles(
        dt_index=pd.date_range(freq="h", **snakemake.params.snapshots, tz="UTC"),
        nodes=nodes,
        weekly_profile=weekly_profile,
    )

    return dsm_profile

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_hourly_heat_demand",
            scope="total",
            clusters=5,
        )
    set_scenario_config(snakemake)

    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    options = snakemake.params.sector

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
                    daily_space_heat_demand.columns, options
                )
        else:
            heat_demand[f"{sector} {use}"] = intraday_year_profile

    heat_demand = pd.concat(heat_demand, axis=1, names=["sector use", "node"])
    dsm_profile = pd.concat(dsm_profile, axis=1, names=["sector use", "node"])

    heat_demand.index.name = "snapshots"

    ds = heat_demand.stack(future_stack=True).to_xarray()

    ds.to_netcdf(snakemake.output.heat_demand)
    dsm_profile.to_csv(snakemake.output.heat_dsm_profile)
