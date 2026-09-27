# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build the district heating share of heat demand at each node for a planning horizon.

Today's country-level district heat share from [build_energy_totals][] is
scaled to each node by the ratio of its share of the country's urban
population to its share of the total population. The maximum share is the
node's urban fraction times the potential set in
`sector.district_heating.potential`, either a single value or one per
country. The gap between today's share and this maximum is closed by the
progress factor of the planning horizon from
`sector.district_heating.progress`. Where today's share exceeds the urban
fraction, the urban fraction is raised to match.

| Column | Description |
|---|---|
| `original district heat share` | Today's share of the country |
| `district fraction of node` | Share of the node's heat demand served by district heating |
| `urban fraction` | Share of the node's population living in urban areas |
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.prepare_sector_network import get

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_district_heat_share",
            horizon="2050",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    investment_year = int(snakemake.wildcards.horizon)

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
    if isinstance(central_fraction, dict):
        # Check if individual district heating shares are given for all countries of the network
        other_countries = set(pop_layout.ct.unique()).difference(
            central_fraction.keys()
        )
        if other_countries:
            default_value = central_fraction.get("default")
            # Default value is required if not all countries are covered
            if default_value is None:
                raise ValueError(
                    "No default district heating potential was provided in the config."
                )
            logger.warning(
                "Some countries do not have a district heating potential defined. "
                f"Using default value {default_value:.2%} for these countries."
            )
            # Fill missing countries with default value from config
            central_fraction = {
                **central_fraction,
                **{ct: default_value for ct in other_countries},
            }
        # Map district heating potentials to bus regions
        central_fraction = pop_layout.ct.map(central_fraction)

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
