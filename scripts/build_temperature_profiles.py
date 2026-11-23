# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build air and soil temperature time series per clustered region.

[atlite.Cutout.temperature](https://atlite.readthedocs.io/en/master/ref_api.html#module-atlite.convert)
and `atlite.Cutout.soil_temperature` read ambient air and soil temperature
from the weather cutout. Grid cells are aggregated to clustered onshore
regions weighted by population, giving the temperature experienced by the
average inhabitant. The profiles serve, among others, as heat source
temperatures in [build_cop_profiles][] and as ambient temperatures in
[build_central_heating_temperature_profiles][].
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

        snakemake = mock_snakemake("build_temperature_profiles")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    dask_kwargs = setup_dask(nprocesses)

    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)

    cutout = load_cutout(snakemake.input.cutout, time=time)

    clustered_regions = (
        gpd.read_file(snakemake.input.onshore_regions).set_index("name").buffer(0)
    )

    I = cutout.indicatormatrix(clustered_regions)  # noqa: E741

    pop_layout = xr.open_dataarray(snakemake.input.pop_layout)

    stacked_pop = pop_layout.stack(spatial=("y", "x"))
    M = I.T.dot(np.diag(I.dot(stacked_pop)))

    nonzero_sum = M.sum(axis=0, keepdims=True)
    nonzero_sum[nonzero_sum == 0.0] = 1.0
    M_tilde = M / nonzero_sum

    temp_air = cutout.temperature(
        matrix=M_tilde.T,
        index=clustered_regions.index,
        dask_kwargs=dask_kwargs,
        show_progress=False,
    )

    temp_air.to_netcdf(snakemake.output.temp_air)

    temp_soil = cutout.soil_temperature(
        matrix=M_tilde.T,
        index=clustered_regions.index,
        dask_kwargs=dask_kwargs,
        show_progress=False,
    )

    temp_soil.to_netcdf(snakemake.output.temp_soil)
