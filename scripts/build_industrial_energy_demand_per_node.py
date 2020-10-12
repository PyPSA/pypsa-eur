
import pandas as pd
import numpy as np

# import EU ratios df as csv
industry_sector_ratios=pd.read_csv(snakemake.input.industry_sector_ratios,
                                   index_col=0)

#material demand per node and industry (kton/a)
nodal_production = pd.read_csv(snakemake.input.industrial_production_per_node,
                               index_col=0)

#energy demand today to get current electricity
nodal_today = pd.read_csv(snakemake.input.industrial_energy_demand_per_node_today,
                          index_col=0)

#final energy consumption per node and industry (TWh/a)
nodal_df = nodal_production.dot(industry_sector_ratios.T)
nodal_df*= 0.001 #GWh -> TWh (ktCO2 -> MtCO2)


rename_sectors = {'elec':'electricity',
                  'biomass':'solid biomass',
                  'heat':'low-temperature heat'}

nodal_df.rename(columns=rename_sectors,inplace=True)

nodal_df["current electricity"] = nodal_today["electricity"]

nodal_df.index.name = "TWh/a (MtCO2/a)"

nodal_df.to_csv(snakemake.output.industrial_energy_demand_per_node,
                float_format='%.2f')
