# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>>
#
# SPDX-License-Identifier: MIT
"""
Build industrial energy demand per model region.

Description
-------

This rule maps the industrial energy demand per country `industrial_energy_demand_per_country_today.csv` to each bus region.
The energy demand per country is multiplied by the mapping value from the file ``industrial_distribution_key_base_s_{clusters}.csv`` between 0 and 1 to get the industrial energy demand per bus.

The unit of the energy demand is TWh/a.
"""

import logging
from itertools import product

import numpy as np
import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

# map JRC/our sectors to hotmaps sector, where mapping exist
sector_mapping = {
    "Electric arc": "EAF",
    "Integrated steelworks": "Integrated steelworks",
    "DRI + Electric arc": "DRI + EAF",
    "Ammonia": "Ammonia",
    "Basic chemicals (without ammonia)": "Chemical industry",
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


def build_nodal_industrial_energy_demand():
    fn = snakemake.input.industrial_energy_demand_per_country_today
    industrial_demand = pd.read_csv(fn, header=[0, 1], index_col=0)

    fn = snakemake.input.industrial_distribution_key
    keys = pd.read_csv(fn, index_col=0)
    keys["country"] = keys.index.str[:2]

    nodal_demand = pd.DataFrame(
        0.0, dtype=float, index=keys.index, columns=industrial_demand.index
    )

    countries = keys.country.unique()
    sectors = industrial_demand.columns.unique(1)

    for country, sector in product(countries, sectors):
        buses = keys.index[keys.country == country]
        mapping = sector_mapping.get(sector, "population")

        key = keys.loc[buses, mapping]
        demand = industrial_demand[country, sector]

        outer = pd.DataFrame(
            np.outer(key, demand), index=key.index, columns=demand.index
        )

        nodal_demand.loc[buses] += outer

    nodal_demand.index.name = "TWh/a"

    nodal_demand.to_csv(snakemake.output.industrial_energy_demand_per_node_today)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_energy_demand_per_node_today",
            clusters=48,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    build_nodal_industrial_energy_demand()
