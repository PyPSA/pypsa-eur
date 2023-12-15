
import pandas as pd

costs = pd.read_csv(snakemake.input.costs, index_col=[0, 1]).sort_index()

if "modifications" in snakemake.input.keys():
    modifications = pd.read_csv(snakemake.input.modifications, index_col=[0, 1]).sort_index()
    costs.loc[modifications.index] = modifications
    print(modifications)
    print( costs.loc[modifications.index])

costs.to_csv(snakemake.output[0])
