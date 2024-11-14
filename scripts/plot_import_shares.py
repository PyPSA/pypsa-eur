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
import pandas as pd
import pypsa
from _helpers import configure_logging

CARRIERS = {
    "AC": "electricity (HVDC)",
    "H2": "hydrogen (pipeline)",
    "NH3": "ammonia (ship)",
    "gas": "methane (ship)",
    "methanol": "methanol (ship)",
    "oil": "Fischer-Tropsch (ship)",
    "HBI": "HBI (ship)",
    "steel": "steel (ship)",
}

COLOR_MAPPING = {
    "electricity (HVDC)": "import hvdc-to-elec",
    "hydrogen (pipeline)": "import pipeline-h2",
    "ammonia (ship)": "import shipping-lnh3",
    "methane (ship)": "import shipping-lch4",
    "methanol (ship)": "import shipping-meoh",
    "Fischer-Tropsch (ship)": "import shipping-ftfuel",
    "steel (ship)": "import shipping-steel",
    "HBI (ship)": "import shipping-hbi",
}

THRESHOLD = 1  # MWh

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_import_shares",
            opts="",
            clusters="115",
            ll="vopt",
            sector_opts="imp",
            planning_horizons="2050",
            configfiles="config/config.20240826-z1.yaml",
        )

    configure_logging(snakemake)

    plt.style.use(["bmh", snakemake.input.rc])

    tech_colors = snakemake.config["plotting"]["tech_colors"]

    n = pypsa.Network(snakemake.input.network)

    eb = n.statistics.energy_balance(nice_names=False).groupby(["carrier", "bus_carrier"]).sum()
    eb = (
        eb.unstack(1)
        .groupby(lambda x: "import" if "import" in x or "external" in x else x)
        .sum()
    )
    supply = eb.where(eb > THRESHOLD).dropna(how="all", axis=0).div(1e6)

    ie = supply.groupby(lambda x: "import" if x == "import" else "domestic").sum()
    ie_rel = ie / ie.sum() * 100
    ie_sum = ie.sum().astype(int)

    imp_mix = ie.loc["import"].copy()
    if "steel" in imp_mix.index:
        imp_mix.loc[["steel", "HBI"]] *= 2.1  # kWh/kg
    imp_mix = pd.DataFrame(
        imp_mix.where(imp_mix > 1e-3)
        .dropna()
        .rename(index=CARRIERS)
        .sort_values(ascending=False)
    )
    imp_mix.columns = ["import mix"]

    fig, (ax, ax_mix) = plt.subplots(
        2, 1, figsize=(6, 5), gridspec_kw=dict(height_ratios=[5.5, 1])
    )
    sel = list(CARRIERS.keys())[::-1]
    ie_rel = ie_rel[sel].rename(columns=CARRIERS).T
    if ie["HBI"].sum() < 0.5:
        ie_rel.loc["HBI (ship)"] = 0.
    ie_rel.index = ie_rel.index.str.split(" ").str[0]
    ie_rel.plot.barh(stacked=True, ax=ax, width=0.8, color=["lightseagreen", "coral"])

    ax.set_ylabel("")
    ax.set_xlabel("domestic share [%]", fontsize=11, color="lightseagreen")
    ax.legend(ncol=2, bbox_to_anchor=(0.55, 1.35))
    ax.grid(axis="y")
    ax.set_xlim(0, 100)

    imp_mix.T.plot.barh(
        ax=ax_mix,
        stacked=True,
        legend=True,
        width=0.9,
        color=[tech_colors[COLOR_MAPPING[i]] for i in imp_mix.index],
    )

    for i, (carrier, twh) in enumerate(ie_sum[sel].items()):
        if carrier.lower().startswith("steel") or carrier.lower().startswith("hbi"):
            unit = "Mt"
        else:
            unit = "TWh"
        ax.text(119, i, f"{twh} {unit}", va="center", ha="right", color="slateblue")
    ax.text(119, i + 1.5, "total\nsupply", va="center", ha="right", color="slateblue")

    secax = ax.secondary_xaxis("top", functions=(lambda x: 100 - x, lambda x: 100 - x))
    secax.set_xlabel("import share [%]", fontsize=11, color="coral")

    total_imp = imp_mix.sum().sum()
    secax_mix = ax_mix.secondary_xaxis(
        "top", functions=(lambda x: x / total_imp * 100, lambda x: x * total_imp * 100)
    )

    ax_mix.text(total_imp * 1.03, 0, str(total_imp.astype(int)) + " TWh", va="center", ha="left", color="slateblue")
    ax_mix.text(total_imp * 1.1, -1.05, "TWh", va="center")
    ax_mix.text(total_imp * 1.1, 1.05, "%", va="center")
    ax_mix.legend(ncol=3, bbox_to_anchor=(1.2, -0.4), title="")

    ticks = range(10, 100, 20)
    ax.set_xticks(ticks, minor=True)
    secax.set_xticks(ticks, minor=True)
    secax_mix.set_xticks(ticks, minor=True)

    ax_mix.set_xlim(0, total_imp)
    ax_mix.grid(axis="y")

    for i in ["top", "right", "left", "bottom"]:
        secax.spines[i].set_visible(False)
        ax.spines[i].set_visible(False)
        ax_mix.spines[i].set_visible(False)
        secax_mix.spines[i].set_visible(False)

    def fmt(x):
        return f"{x:.0f}%" if x > 1 else ""

    for container in ax.containers[:2]:
        ax.bar_label(container, label_type="center", color="white", fmt=fmt)

    def fmt(x):
        return f"{x:.0f}" if x > 80 else ""

    for container in ax_mix.containers:
        ax_mix.bar_label(container, label_type="center", color="white", fmt=fmt)

    for fn in snakemake.output:
        plt.savefig(fn, bbox_inches="tight")
