# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build district heat shares at each node, depending on investment year.

Inputs:
-------
- `resources/<run_name>/pop_layout.csv`: Population layout for each node: Total, urban and rural population.
- `resources/<run_name>/district_heat_share.csv`: Historical district heat share at each country. Output of `scripts/build_energy_totals.py`.

Outputs:
--------
- `resources/<run_name>/district_heat_share.csv`: District heat share at each node, potential for each investment year.

Notes
-----
- The district heat share is calculated as the share of urban population at each node, multiplied by the share of district heating in the respective country.
- The `sector.district_heating.potential` setting defines the max. district heating share.
- The max. share of district heating is increased by a progress factor, depending on the investment year (See `sector.district_heating.progress` setting).
"""

import logging

import pandas as pd
from _helpers import configure_logging, set_scenario_config
from prepare_sector_network import get

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_district_heat_share",
            clusters=60,
            planning_horizons="2050",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    investment_year = int(snakemake.wildcards.planning_horizons)

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    year = str(snakemake.params.energy_totals_year)
    district_heat_share = pd.read_csv(snakemake.input.district_heat_share, index_col=0)[
        year
    ]

    # make ct-based share nodal
    district_heat_share = district_heat_share.reindex(pop_layout.ct).fillna(0)
    district_heat_share.index = pop_layout.index

    # total urban population per country
    ct_urban = pop_layout.urban.groupby(pop_layout.ct).sum()

    # distribution of urban population within a country
    pop_layout["urban_ct_fraction"] = pop_layout.urban / pop_layout.ct.map(ct_urban.get)

    # fraction of node that is urban
    urban_fraction = pop_layout.urban / pop_layout[["rural", "urban"]].sum(axis=1)

    # maximum potential of urban demand covered by district heating
    central_fraction = snakemake.config["sector"]["district_heating"]["potential"]

    # district heating share at each node
    dist_fraction_node = (
        district_heat_share * pop_layout["urban_ct_fraction"] / pop_layout["fraction"]
    )

    # if district heating share larger than urban fraction -> set urban
    # fraction to district heating share
    urban_fraction = pd.concat([urban_fraction, dist_fraction_node], axis=1).max(axis=1)

    # difference of max potential and today's share of district heating
    diff = ((urban_fraction * central_fraction) - dist_fraction_node).clip(lower=0)
    progress = get(
        snakemake.config["sector"]["district_heating"]["progress"], investment_year
    )
    dist_fraction_node += diff * progress
    logger.info(
        f"Increase district heating share by a progress factor of {progress:.2%} "
        f"resulting in new average share of {dist_fraction_node.mean():.2%}"
    )

    df = pd.DataFrame(
        {
            "original district heat share": district_heat_share,
            "district fraction of node": dist_fraction_node,
            "urban fraction": urban_fraction,
        },
        dtype=float,
    )

    df.to_csv(snakemake.output.district_heat_share)
