# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot average capacity factor map for selected sector-coupling technologies.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from plot_choropleth_capacity_factors import plot_choropleth

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_choropleth_capacity_factors_sector",
            simpl="",
            clusters=128,
            configfiles=["../../config/config.test.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    df = pd.DataFrame()

    ds = xr.open_dataset(snakemake.input.cop_soil)
    df["ground-sourced heat pump"] = ds["soil temperature"].mean("time").to_pandas()

    ds = xr.open_dataset(snakemake.input.cop_air)
    df["air-sourced heat pump"] = ds["temperature"].mean("time").to_pandas()

    plot_choropleth(
        df,
        regions_onshore,
        "ground-sourced heat pump",
        cmap="Greens",
        vmax=4,
        vmin=2,
        label="mean COP [-]",
        dir=snakemake.output[0],
    )

    plot_choropleth(
        df,
        regions_onshore,
        "air-sourced heat pump",
        cmap="Greens",
        vmax=4,
        vmin=2,
        label="mean COP [-]",
        dir=snakemake.output[0],
    )
