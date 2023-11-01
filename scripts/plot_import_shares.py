# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates stacked bar charts of import shares per carrier.
"""

import logging

logger = logging.getLogger(__name__)

import matplotlib.pyplot as plt
import pypsa
from _helpers import configure_logging

CARRIERS = {
    "AC": "electricity",
    "H2": "hydrogen",
    "NH3": "ammonia",
    "gas": "methane",
    "methanol": "methanol",
    "oil": "Fischer-Tropsch",
    "steel": "steel",
}

THRESHOLD = 1 # MWh

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_import_shares",
            simpl="",
            opts="",
            clusters="100",
            ll="vopt",
            sector_opts="Co2L0-2190SEG-T-H-B-I-S-A-imp",
            planning_horizons="2050",
            configfiles="../../config/config.20230926-zecm.yaml",
        )

    configure_logging(snakemake)

    plt.style.use(['bmh', snakemake.input.rc])

    n = pypsa.Network(snakemake.input.network)

    eb = n.statistics.energy_balance().groupby(["carrier", "bus_carrier"]).sum()
    supply = eb.where(eb > THRESHOLD).dropna().div(1e6)

    ie = supply.unstack(1).groupby(lambda x: "import" if "import" in x or "external" in x else "domestic").sum()
    ie = ie / ie.sum() * 100


    fig, ax = plt.subplots(figsize=(6, 4))
    sel = list(CARRIERS.keys())[::-1]
    ie[sel].rename(columns=CARRIERS).T.plot.barh(
        stacked=True,
        ax=ax,
        color=['lightseagreen', 'coral']
    )

    plt.ylabel("")
    plt.xlabel("domestic share [%]", fontsize=11, color='lightseagreen')
    plt.legend(ncol=2, bbox_to_anchor=(0.35, 1.25))
    plt.grid(axis="y")
    plt.xlim(0,100)

    secax = ax.secondary_xaxis('top', functions=(lambda x: 100 - x, lambda x: 100 -x))
    secax.set_xlabel('import share [%]', fontsize=11, color='coral')

    ticks = range(10, 100, 20)
    ax.set_xticks(ticks, minor=True)
    secax.set_xticks(ticks, minor=True)

    for i in ["top", "right", "left", "bottom"]:
        secax.spines[i].set_visible(False)
        ax.spines[i].set_visible(False)

    def fmt(x):
        return f"{x:.0f}%" if x > 1 else ""

    for container in ax.containers[:2]:
        ax.bar_label(container, label_type='center', color='white', fmt=fmt)

    for fn in snakemake.output:
        plt.savefig(fn, bbox_inches="tight")
