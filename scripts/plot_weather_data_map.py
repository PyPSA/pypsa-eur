# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot time-averaged cutout weather data on a map.
"""

import logging
import atlite
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from _helpers import configure_logging


logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_weather_data_map",
            configfiles=["../../config/config.test.yaml"]
        )
    configure_logging(snakemake)

    plt.style.use(snakemake.input.rc)

    data = dict(
        irradiation=dict(label=r"Mean Direct Solar Irradiance [W/m$^2$]", cmap="Oranges"),
        runoff=dict(label=r"Total Runoff [m]", cmap="Greens"),
        temperature=dict(label=r"Mean Temperatures [°C]", cmap="Reds"),
        wind=dict(label=r"Mean Wind Speeds [m/s]", cmap="Blues"),
    )

    cutout = atlite.Cutout(snakemake.input.cutout)

    data["irradiation"]["data"] = cutout.data.influx_direct.mean("time")
    data["runoff"]["data"] = cutout.data.runoff.sum("time")
    data["temperature"]["data"] = (cutout.data.temperature.mean("time") - 273.15)
    data["wind"]["data"] = cutout.data.wnd100m.mean("time")

    for k, v in data.items():

        logger.info(f"Plotting weather data map for variable '{k}'")

        fig, ax = plt.subplots(
            figsize=(7, 7), subplot_kw={"projection": ccrs.EqualEarth()}
        )

        v["data"].plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=v["cmap"],
            linewidths=0,
            cbar_kwargs={"label": v["label"], "shrink": 0.6},
        )

        ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
        ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)
        ax.gridlines(linestyle=":")
        ax.axis("off")

        for fn in snakemake.output[k]:
            plt.savefig(fn)
