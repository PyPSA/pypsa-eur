"""Build industrial production per node."""

import pandas as pd
from itertools import product

# map JRC/our sectors to hotmaps sector, where mapping exist
sector_mapping = {
    'Electric arc': 'Iron and steel',
    'Integrated steelworks': 'Iron and steel',
    'DRI + Electric arc': 'Iron and steel',
    'Ammonia': 'Chemical industry',
    'HVC': 'Chemical industry',
    'HVC (mechanical recycling)': 'Chemical industry',
    'HVC (chemical recycling)': 'Chemical industry',
    'Methanol': 'Chemical industry',
    'Chlorine': 'Chemical industry',
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


def build_nodal_industrial_production():

    fn = snakemake.input.industrial_production_per_country_tomorrow
    industrial_production = pd.read_csv(fn, index_col=0)

    fn = snakemake.input.industrial_distribution_key
    keys = pd.read_csv(fn, index_col=0)
    keys["country"] = keys.index.str[:2]

    nodal_production = pd.DataFrame(index=keys.index,
                                    columns=industrial_production.columns,
                                    dtype=float)

    countries = keys.country.unique()
    sectors = industrial_production.columns

    for country, sector in product(countries, sectors):

        buses = keys.index[keys.country == country]
        mapping = sector_mapping.get(sector, "population")

        key = keys.loc[buses, mapping]
        nodal_production.loc[buses, sector] = industrial_production.at[country, sector] * key

    nodal_production.to_csv(snakemake.output.industrial_production_per_node)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('build_industrial_production_per_node',
            simpl='',
            clusters=48,
        )

    build_nodal_industrial_production()
