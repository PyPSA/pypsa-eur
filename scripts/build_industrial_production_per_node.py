
import pandas as pd

def build_nodal_industrial_production():

    industrial_production = pd.read_csv(snakemake.input.industrial_production_per_country_tomorrow,
                                        index_col=0)

    distribution_keys = pd.read_csv(snakemake.input.industrial_distribution_key,
                                        index_col=0)
    distribution_keys["country"] = distribution_keys.index.str[:2]

    nodal_industrial_production = pd.DataFrame(index=distribution_keys.index,
                                               columns=industrial_production.columns,
                                               dtype=float)

    #map JRC/our sectors to hotmaps sector, where mapping exist
    sector_mapping = {'Electric arc' : 'Iron and steel',
                      'Integrated steelworks' : 'Iron and steel',
                      'DRI + Electric arc' : 'Iron and steel',
                      'Ammonia' : 'Chemical industry',
                      'Basic chemicals (without ammonia)' : 'Chemical industry',
                      'Other chemicals' : 'Chemical industry',
                      'Pharmaceutical products etc.' : 'Chemical industry',
                      'Cement' : 'Cement',
                      'Ceramics & other NMM' : 'Non-metallic mineral products',
                      'Glass production' : 'Glass',
                      'Pulp production' : 'Paper and printing',
                      'Paper production' : 'Paper and printing',
                      'Printing and media reproduction' : 'Paper and printing',
                      'Alumina production' : 'Non-ferrous metals',
                      'Aluminium - primary production' : 'Non-ferrous metals',
                      'Aluminium - secondary production' : 'Non-ferrous metals',
                      'Other non-ferrous metals' : 'Non-ferrous metals',
    }

    for c in distribution_keys.country.unique():
        buses = distribution_keys.index[distribution_keys.country == c]
        for sector in industrial_production.columns:
            distribution_key = distribution_keys.loc[buses,sector_mapping.get(sector,"population")]
            nodal_industrial_production.loc[buses,sector] = industrial_production.at[c,sector]*distribution_key

    nodal_industrial_production.to_csv(snakemake.output.industrial_production_per_node)

if __name__ == "__main__":

    build_nodal_industrial_production()
