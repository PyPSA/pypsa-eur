#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Created on Fri Jun 30 10:50:53 2023.

@author: fabian
"""

import matplotlib.pyplot as plt
import pypsa
import seaborn as sns

sns.set_theme("paper", style="whitegrid")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_statistics",
            simpl="",
            opts="Co2L-3h",
            clusters="37c",
            ll="v1.0",
        )


n = pypsa.Network(snakemake.network)


n.loads.carrier = "load"
n.carriers.loc["load", ["nice_name", "color"]] = "Load", "darkred"
colors = n.carriers.set_index("nice_name").color.where(lambda s: s != "", "lightgrey")


def rename_index(ds):
    return ds.set_axis(ds.index.map(lambda x: f"{x[1]}\n({x[0].lower()})"))


def plot_static_per_carrier(ds, ax, drop_zero=True):
    if drop_zero:
        ds = ds[ds != 0]
    ds = ds.dropna()
    c = colors[ds.index.get_level_values("carrier")]
    ds = ds.pipe(rename_index)
    label = f"{ds.attrs['name']} [{ds.attrs['unit']}]"
    ds.plot.barh(color=c.values, xlabel=label, ax=ax)


fig, ax = plt.subplots()
ds = n.statistics.capacity_factor().dropna()
plot_static_per_carrier(ds, ax)
# fig.savefig("")

fig, ax = plt.subplots()
ds = n.statistics.installed_capacity().dropna()
ds = ds.drop("Line")
ds = ds / 1e3
ds.attrs["unit"] = "GW"
plot_static_per_carrier(ds, ax)
# fig.savefig("")


fig, ax = plt.subplots()
ds = n.statistics.optimal_capacity()
ds = ds.drop("Line")
ds = ds / 1e3
ds.attrs["unit"] = "GW"
plot_static_per_carrier(ds, ax)
# fig.savefig("")


fig, ax = plt.subplots()
ds = n.statistics.capex()
plot_static_per_carrier(ds, ax)
# fig.savefig("")

fig, ax = plt.subplots()
ds = n.statistics.curtailment()
plot_static_per_carrier(ds, ax)
# fig.savefig("")

fig, ax = plt.subplots()
ds = n.statistics.supply()
ds = ds.drop("Line")
ds = ds / 1e6
ds.attrs["unit"] = "TWh"
plot_static_per_carrier(ds, ax)
# fig.savefig("")


fig, ax = plt.subplots()
ds = n.statistics.withdrawal()
ds = ds.drop("Line")
ds = ds / -1e6
ds.attrs["unit"] = "TWh"
plot_static_per_carrier(ds, ax)
# fig.savefig("")
