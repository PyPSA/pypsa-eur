
import pandas as pd

industrial_production = pd.read_csv(snakemake.input.industrial_production_per_country,
                                    index_col=0)

total_steel = industrial_production[["Integrated steelworks","Electric arc"]].sum(axis=1)

fraction_primary_stays_primary = snakemake.config["industry"]["St_primary_fraction"]*total_steel.sum()/industrial_production["Integrated steelworks"].sum()

industrial_production.insert(2, "DRI + Electric arc",
                             fraction_primary_stays_primary*industrial_production["Integrated steelworks"])

industrial_production["Electric arc"] = total_steel - industrial_production["DRI + Electric arc"]
industrial_production["Integrated steelworks"] = 0.


total_aluminium = industrial_production[["Aluminium - primary production","Aluminium - secondary production"]].sum(axis=1)

fraction_primary_stays_primary = snakemake.config["industry"]["Al_primary_fraction"]*total_aluminium.sum()/industrial_production["Aluminium - primary production"].sum()

industrial_production["Aluminium - primary production"] = fraction_primary_stays_primary*industrial_production["Aluminium - primary production"]
industrial_production["Aluminium - secondary production"] = total_aluminium - industrial_production["Aluminium - primary production"]

industrial_production["Basic chemicals (without ammonia)"] *= snakemake.config["industry"]['HVC_primary_fraction']


industrial_production.to_csv(snakemake.output.industrial_production_per_country_tomorrow,
                             float_format='%.2f')
