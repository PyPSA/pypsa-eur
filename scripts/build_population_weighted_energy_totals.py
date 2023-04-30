# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Distribute country-level energy demands by population.
"""

import pandas as pd

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_population_weighted_energy_totals",
            weather_year="",
            simpl="",
            clusters=48,
        )

    config = snakemake.config["energy"]
    data_year = int(config["energy_totals_year"])
    if snakemake.wildcards.weather_year and snakemake.wildcards.kind == "heat":
        data_year = int(snakemake.wildcards.weather_year)

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    totals = pd.read_csv(snakemake.input.totals, index_col=[0, 1])
    totals = totals.xs(data_year, level="year")

    nodal_totals = totals.loc[pop_layout.ct].fillna(0.0)
    nodal_totals.index = pop_layout.index
    nodal_totals = nodal_totals.multiply(pop_layout.fraction, axis=0)

    nodal_totals.to_csv(snakemake.output[0])
