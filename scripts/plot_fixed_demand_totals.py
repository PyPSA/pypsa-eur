# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot fixed demand totals.
"""

import matplotlib.pyplot as plt
import pypsa
from pypsa.descriptors import get_switchable_as_dense as as_dense

COLUMNS = {
    'EV battery': "electricity",
    'H2': "hydrogen",
    'NH3': "ammonia",
    'agriculture machinery oil': "oil",
    'gas for industry': "methane",
    'industry methanol': "methanol",
    'kerosene for aviation': "oil",
    'low voltage': "electricity",
    'naphtha for industry': "oil",
    'rural heat': "heat",
    'shipping methanol': "methanol",
    'solid biomass for industry': "biomass",
    'steel': "steel",
    'urban central heat': "heat",
    'urban decentral heat': "heat",
    '': "electricity",
}

INDEX = {
    'H2 for industry': "industry",
    'NH3': "fertilizers",
    'agriculture electricity': "agriculture",
    'agriculture heat': "agriculture",
    'agriculture machinery oil': "agriculture",
    'electricity': "residential/services",
    'gas for industry': "industry",
    'industry electricity': "industry",
    'industry methanol': "industry",
    'kerosene for aviation': "aviation",
    'land transport EV': "road transport",
    'low-temperature heat for industry': "industry",
    'naphtha for industry': "feedstock",
    'rural heat': "rural heat",
    'shipping methanol': "shipping",
    'solid biomass for industry': "industry",
    'steel': "industry",
    'urban central heat': "urban heat",
    'urban decentral heat': "urban heat",
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_fixed_demand_totals",
            simpl="",
            clusters=115,
            ll="vopt",
            opts="",
            sector_opts="",
            planning_horizons=2050,
            configfiles="config/config.20240826-z1.yaml",
        )

    plt.style.use(["bmh", snakemake.input.rc])

    n = pypsa.Network(snakemake.input.network)

    w = n.snapshot_weightings.generators
    loads = w @ as_dense(n, "Load", "p_set")

    stack = (
        loads.groupby([n.loads.carrier, n.loads.bus.map(n.buses.carrier)])
        .sum()
        .div(1e6)
        .unstack()
    )

    stack = stack.drop("process emissions", axis=0).drop("process emissions", axis=1)

    stack = stack.groupby(INDEX).sum().T.groupby(COLUMNS).sum()

    stack.loc["steel"] *= 2.1

    stack = stack.loc[stack.sum(axis=1).sort_values().index]

    fig, ax = plt.subplots(figsize=(8, 2.8))

    stack.plot.barh(ax=ax, stacked=True, width=0.8)

    ax.set_xlabel("Final energy and non-energy demand [TWh/a]")
    ax.set_ylabel("")

    ax.grid(False)

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Used for...")


    def fmt(x):
        return f"{x:.0f}" if x > 200 else ""

    for c in ax.containers:
        ax.bar_label(c, label_type="center", color='white', fmt=fmt)

    plt.tight_layout()

    ax.set_xlim(-50, 5000)
    ax.set_xticks(range(0, 5001, 500))
    ax.set_xticks(range(0, 5001, 250), minor=True)

    
    for i, (idx, total) in enumerate(stack.sum(axis=1).items()):
        ax.text(total + 50, i, f"{total:.0f}", va='center')

    for fn in snakemake.output:
        plt.savefig(fn, bbox_inches="tight")
