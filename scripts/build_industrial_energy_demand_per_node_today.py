
import pandas as pd
import numpy as np

def build_nodal_demand():

    industrial_demand = pd.read_csv(snakemake.input.industrial_energy_demand_per_country_today,
                                    header=[0,1],
                                    index_col=0)

    distribution_keys = pd.read_csv(snakemake.input.industrial_distribution_key,
                                        index_col=0)
    distribution_keys["country"] = distribution_keys.index.str[:2]

    nodal_demand = pd.DataFrame(0.,
                                index=distribution_keys.index,
                                columns=industrial_demand.index,
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
        for sector in industrial_demand.columns.levels[1]:
            distribution_key = distribution_keys.loc[buses,sector_mapping.get(sector,"population")]
            demand = industrial_demand[c,sector]
            outer = pd.DataFrame(np.outer(distribution_key,demand),index=distribution_key.index,columns=demand.index)
            nodal_demand.loc[buses] += outer

    nodal_demand.index.name = "TWh/a"

    nodal_demand.to_csv(snakemake.output.industrial_energy_demand_per_node_today)

if __name__ == "__main__":

    build_nodal_demand()
