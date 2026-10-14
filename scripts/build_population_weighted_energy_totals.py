# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Distribute country-level energy or heat totals to clustered model regions by
population.

Annual totals per country are averaged over the requested years (the
reference year for energy totals, the weather years of the snapshots for heat
totals), mapped to each region via its country and scaled by the region's
share of national population from the clustered population layout. The
regional totals are the basis for building sectoral demand time series.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, get_snapshots, set_scenario_config

idx = pd.IndexSlice

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_population_weighted_energy_totals",
            kind="heat",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    config = snakemake.config["energy"]

    if snakemake.wildcards.kind == "heat":
        snapshots = get_snapshots(
            snakemake.params.snapshots, snakemake.params.drop_leap_day
        )
        data_years = snapshots.year.unique()
    else:
        data_years = int(config["energy_totals_year"])

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    totals = pd.read_csv(snakemake.input.energy_totals, index_col=[0, 1])
    totals = totals.loc[idx[:, data_years], :].groupby("country").mean()

    nodal_totals = totals.loc[pop_layout.ct].fillna(0.0)
    nodal_totals.index = pop_layout.index
    nodal_totals = nodal_totals.multiply(pop_layout.fraction, axis=0)

    nodal_totals.to_csv(snakemake.output[0])
