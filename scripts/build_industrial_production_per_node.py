# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build industrial production per model region.

Description
-------

This rule maps the industrial production per country from a certain time horizon to each bus region.
The mapping file provides a value between 0 and 1 for each bus and industry subcategory, indicating the share of the country's production of that sector in that bus.
The industrial production per country is multiplied by the mapping value to get the industrial production per bus.
The unit of the production is kt/a.
"""

import logging
from itertools import product

import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

# map JRC/our sectors to hotmaps sector, where mapping exist
sector_mapping = {
    "Electric arc": "EAF",
    "Integrated steelworks": "Integrated steelworks",
    "DRI + Electric arc": "DRI + EAF",
    "Ammonia": "Ammonia",
    "HVC": "Chemical industry",
    "HVC (mechanical recycling)": "Chemical industry",
    "HVC (chemical recycling)": "Chemical industry",
    "Methanol": "Chemical industry",
    "Chlorine": "Chemical industry",
    "Other chemicals": "Chemical industry",
    "Pharmaceutical products etc.": "Chemical industry",
    "Cement": "Cement",
    "Ceramics & other NMM": "Non-metallic mineral products",
    "Glass production": "Glass",
    "Pulp production": "Paper and printing",
    "Paper production": "Paper and printing",
    "Printing and media reproduction": "Paper and printing",
    "Alumina production": "Non-ferrous metals",
    "Aluminium - primary production": "Non-ferrous metals",
    "Aluminium - secondary production": "Non-ferrous metals",
    "Other non-ferrous metals": "Non-ferrous metals",
}


def build_nodal_industrial_production():
    fn = snakemake.input.industrial_production_per_country_tomorrow
    industrial_production = pd.read_csv(fn, index_col=0)

    fn = snakemake.input.industrial_distribution_key
    keys = pd.read_csv(fn, index_col=0)
    keys["country"] = keys.index.str[:2]

    nodal_production = pd.DataFrame(
        index=keys.index, columns=industrial_production.columns, dtype=float
    )

    countries = keys.country.unique()
    sectors = industrial_production.columns

    for country, sector in product(countries, sectors):
        buses = keys.index[keys.country == country]
        mapping = sector_mapping.get(sector, "population")

        key = keys.loc[buses, mapping]
        nodal_production.loc[buses, sector] = (
            industrial_production.at[country, sector] * key
        )

    nodal_production.to_csv(snakemake.output.industrial_production_per_node)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industrial_production_per_node", clusters=48)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    build_nodal_industrial_production()
