# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build daily heat demand time series per clustered region from heating degree days.

Ambient temperature from the weather cutout is converted to heat demand with
[atlite.Cutout.heat_demand](https://atlite.readthedocs.io/en/master/ref_api.html#module-atlite.convert),
which counts the degrees by which the daily mean temperature falls below a
threshold. Grid cells are aggregated to clustered onshore regions weighted by
population and the result is kept at daily resolution. The daily profile is
disaggregated to hours in [build_hourly_heat_demand][].
"""

import logging

import geopandas as gpd
import numpy as np
import xarray as xr

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

        snakemake = mock_snakemake(
            "build_daily_heat_demand",
            scope="total",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    dask_kwargs = setup_dask(nprocesses)

    cutout_name = snakemake.input.cutout

    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    daily = get_snapshots(
        snakemake.params.snapshots,
        snakemake.params.drop_leap_day,
        freq="D",
    )

    cutout = load_cutout(cutout_name, time=time)

    clustered_regions = (
        gpd.read_file(snakemake.input.onshore_regions).set_index("name").buffer(0)
    )

    I = cutout.indicatormatrix(clustered_regions)  # noqa: E741

    pop_layout = xr.open_dataarray(snakemake.input.pop_layout)

    stacked_pop = pop_layout.stack(spatial=("y", "x"))
    M = I.T.dot(np.diag(I.dot(stacked_pop)))

    heat_demand = cutout.heat_demand(
        matrix=M.T,
        index=clustered_regions.index,
        dask_kwargs=dask_kwargs,
        show_progress=False,
    ).sel(time=daily)

    heat_demand.to_netcdf(snakemake.output.heat_demand)
