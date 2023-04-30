# -*- coding: utf-8 -*-
"""
Approximate heat demand for all weather years.
"""

from itertools import product

import pandas as pd
from numpy.polynomial import Polynomial

idx = pd.IndexSlice


def approximate_heat_demand(energy_totals, hdd):
    if isinstance(hdd, str):
        hdd = pd.read_csv(hdd, index_col=0).T
    hdd.index = hdd.index.astype(int)

    demands = {}

    for kind, sector in product(["total", "electricity"], ["services", "residential"]):
        row = idx[:, 2007:2015]
        col = f"{kind} {sector} space"
        demand = energy_totals.loc[row, col].unstack(0)

        demand_approx = {}

        for c in countries:
            Y = demand[c].dropna()
            X = hdd.loc[Y.index, c]

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
        from helper import mock_snakemake

        snakemake = mock_snakemake("build_energy_totals")

    hdd = pd.read_csv(snakemake.input.hdd, index_col=0).T

    energy_totals = pd.read_csv(snakemake.input.energy_totals, index_col=[0, 1])

    countries = hdd.columns

    heat_demand = approximate_heat_demand(energy_totals, hdd)

    heat_demand.to_csv(snakemake.output.heat_totals)
