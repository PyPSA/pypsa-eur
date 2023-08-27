# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot heatmap time series (resources).
"""

import logging
import pypsa
import matplotlib.pyplot as plt
from multiprocessing import Pool

logger = logging.getLogger(__name__)

from plot_heatmap_timeseries import unstack_day_hour, plot_heatmap


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_heatmap_timeseries_resources",
            simpl="",
            clusters=100,
            ll="v1.5",
            opts="",
            sector_opts="Co2L0-2190SEG-T-H-B-I-S-A",
            planning_horizons=2050,
            configfiles="../../config/config.100n-seg.yaml"
        )

    plt.style.use(["bmh", snakemake.input.rc])

    n = pypsa.Network(snakemake.input.network)

    dir = snakemake.output[0]
    snapshots = snakemake.config["snapshots"]
    carriers = n.carriers

    # WWS capacity factors
    data = n.generators_t.p_max_pu.groupby(n.generators.carrier, axis=1).mean() * 100

    def process_capacity_factors(carrier, s):
        df = unstack_day_hour(s, snapshots)
        label = "capacity factor [%]"
        title = carriers.nice_name.get(carrier, carrier)
        plot_heatmap(
            df,
            cmap="Greens",
            vmin=0,
            vmax=100,
            label=label,
            title=title,
            tag="capacity_factor",
            dir=dir,
        )

    with Pool(processes=snakemake.threads) as pool:
        pool.starmap(process_capacity_factors, data.items())
    
    # heat pump COPs
    data = n.links_t.efficiency.groupby(n.links.carrier, axis=1).mean()

    def process_cops(carrier, s):
        df = unstack_day_hour(s, snapshots)
        label = "coefficient of performance [-]"
        title = carriers.nice_name.get(carrier, carrier)
        plot_heatmap(
            df,
            cmap="Greens",
            vmin=2,
            vmax=4,
            label=label,
            title=title,
            cbar_kws=dict(extend='both'),
            tag="cop",
            dir=dir,
        )

    with Pool(processes=snakemake.threads) as pool:
        pool.starmap(process_cops, data.items())

    # EV availability
    data = n.links_t.p_max_pu.groupby(n.links.carrier, axis=1).mean() * 100

    def process_availabilities(carrier, s):
        df = unstack_day_hour(s, snapshots)
        label = "availability [%]"
        title = carriers.nice_name.get(carrier, carrier)
        plot_heatmap(
            df,
            cmap="Greens",
            vmin=60,
            vmax=100,
            label=label,
            title=title,
            tag="availability",
            dir=dir,
        )

    with Pool(processes=snakemake.threads) as pool:
        pool.starmap(process_availabilities, data.items())

