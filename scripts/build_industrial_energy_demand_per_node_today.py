"""Build industrial energy demand per node."""

import pandas as pd
import numpy as np
from itertools import product

# map JRC/our sectors to hotmaps sector, where mapping exist
sector_mapping = {
    'Electric arc': 'Iron and steel',
    'Integrated steelworks': 'Iron and steel',
    'DRI + Electric arc': 'Iron and steel',
    'Ammonia': 'Chemical industry',
    'Basic chemicals (without ammonia)': 'Chemical industry',
    'Other chemicals': 'Chemical industry',
    'Pharmaceutical products etc.': 'Chemical industry',
    'Cement': 'Cement',
    'Ceramics & other NMM': 'Non-metallic mineral products',
    'Glass production': 'Glass',
    'Pulp production': 'Paper and printing',
    'Paper production': 'Paper and printing',
    'Printing and media reproduction': 'Paper and printing',
    'Alumina production': 'Non-ferrous metals',
    'Aluminium - primary production': 'Non-ferrous metals',
    'Aluminium - secondary production': 'Non-ferrous metals',
    'Other non-ferrous metals': 'Non-ferrous metals',
}


def build_nodal_industrial_energy_demand():

    fn = snakemake.input.industrial_energy_demand_per_country_today
    industrial_demand = pd.read_csv(fn, header=[0, 1], index_col=0)

    fn = snakemake.input.industrial_distribution_key
    keys = pd.read_csv(fn, index_col=0)
    keys["country"] = keys.index.str[:2]

    nodal_demand = pd.DataFrame(0., dtype=float,
                                index=keys.index,
                                columns=industrial_demand.index)
                                
    countries = keys.country.unique()
    sectors = industrial_demand.columns.levels[1]

    for country, sector in product(countries, sectors):
        
        buses = keys.index[keys.country == country]
        mapping = sector_mapping.get(sector, 'population')

        key = keys.loc[buses, mapping]
        demand = industrial_demand[country, sector]

        outer = pd.DataFrame(np.outer(key, demand),
                            index=key.index,
                            columns=demand.index)

        nodal_demand.loc[buses] += outer

    nodal_demand.index.name = "TWh/a"

    nodal_demand.to_csv(snakemake.output.industrial_energy_demand_per_node_today)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_industrial_energy_demand_per_node_today',
            simpl='',
            clusters=48,
        )

    build_nodal_industrial_energy_demand()
