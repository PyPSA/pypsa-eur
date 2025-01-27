# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Distribute country-level energy demands by population.
"""

import logging

import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_population_weighted_energy_totals",
            kind="heat",
            clusters=60,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    config = snakemake.config["energy"]

    if snakemake.wildcards.kind == "heat":
        years = pd.date_range(freq="h", **snakemake.params.snapshots).year.unique()
        assert len(years) == 1, "Currently only works for single year."
        data_year = years[0]
    else:
        data_year = int(config["energy_totals_year"])

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    totals = pd.read_csv(snakemake.input.energy_totals, index_col=[0, 1])
    totals = totals.xs(data_year, level="year")

    nodal_totals = totals.loc[pop_layout.ct].fillna(0.0)
    nodal_totals.index = pop_layout.index
    nodal_totals = nodal_totals.multiply(pop_layout.fraction, axis=0)

    nodal_totals.to_csv(snakemake.output[0])
