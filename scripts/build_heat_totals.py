# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Approximate heat demand for all weather years.
"""

from itertools import product

import pandas as pd
from numpy.polynomial import Polynomial

idx = pd.IndexSlice


def approximate_heat_demand(energy_totals, hdd):

    countries = hdd.columns.intersection(energy_totals.index.levels[0])

    demands = {}

    for kind, sector in product(["total", "electricity"], ["services", "residential"]):
        # reduced number years (2007-2021) for regression because it implicitly
        # assumes a constant building stock
        row = idx[:, 2007:2021]
        col = f"{kind} {sector} space"
        demand = energy_totals.loc[row, col].unstack(0)

        # ffill for GB in 2020- and bfill for CH 2007-2009
        # compromise to have more years available for the fit
        demand = demand.ffill(axis=0).bfill(axis=0)

        demand_approx = {}

        for c in countries:
            Y = demand[c].dropna()
            X = hdd.loc[Y.index, c]

            # Sometimes (looking at you, Switzerland) we only have
            # _one_ year of heating data to base the prediction on. In
            # this case we add a point at 0, 0 to make a "polynomial"
            # fit work.
            if len(X) == len(Y) == 1:
                X.loc[-1] = 0
                Y.loc[-1] = 0

            to_predict = hdd.index.difference(Y.index)
            X_pred = hdd.loc[to_predict, c]

            p = Polynomial.fit(X, Y, 1)
            Y_pred = p(X_pred)

            demand_approx[c] = pd.Series(Y_pred, index=to_predict)

        demand_approx = pd.DataFrame(demand_approx)
        demand_approx = pd.concat([demand, demand_approx]).sort_index()
        demands[f"{kind} {sector} space"] = demand_approx.groupby(
            demand_approx.index
        ).sum()

    demands = pd.concat(demands).unstack().T.clip(lower=0)
    demands.index.names = ["country", "year"]

    return demands


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_heat_totals")

    hdd = pd.read_csv(snakemake.input.hdd, index_col=0).T
    hdd.index = hdd.index.astype(int)

    energy_totals = pd.read_csv(snakemake.input.energy_totals, index_col=[0, 1])

    heat_demand = approximate_heat_demand(energy_totals, hdd)

    heat_demand.to_csv(snakemake.output.heat_totals)
