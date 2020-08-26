
import pandas as pd

industrial_production = pd.read_csv(snakemake.input.industrial_production_per_country,
                                    index_col=0)

total_steel = industrial_production[["Integrated steelworks","Electric arc"]].sum(axis=1)

industrial_production.insert(2, "DRI + Electric arc",
                             snakemake.config["industry"]["St_primary_fraction"]*total_steel)
industrial_production["Electric arc"] = (1 - snakemake.config["industry"]["St_primary_fraction"])*total_steel
industrial_production["Integrated steelworks"] = 0.


total_aluminium = industrial_production[["Aluminium - primary production","Aluminium - secondary production"]].sum(axis=1)

industrial_production["Aluminium - primary production"] = snakemake.config["industry"]["Al_primary_fraction"]*total_aluminium
industrial_production["Aluminium - secondary production"] = (1 - snakemake.config["industry"]["Al_primary_fraction"])*total_aluminium


industrial_production.to_csv(snakemake.output.industrial_production_per_country_tomorrow,
                             float_format='%.2f')
