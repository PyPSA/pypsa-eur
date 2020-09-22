
import pandas as pd

idx = pd.IndexSlice

def build_industrial_demand():
    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout,index_col=0)
    pop_layout["ct"] = pop_layout.index.str[:2]
    ct_total = pop_layout.total.groupby(pop_layout["ct"]).sum()
    pop_layout["ct_total"] = pop_layout["ct"].map(ct_total)
    pop_layout["fraction"] = pop_layout["total"]/pop_layout["ct_total"]

    industrial_demand_per_country = pd.read_csv(snakemake.input.industrial_demand_per_country,index_col=0)

    industrial_demand = industrial_demand_per_country.loc[pop_layout.ct].fillna(0.)
    industrial_demand.index = pop_layout.index
    industrial_demand = industrial_demand.multiply(pop_layout.fraction,axis=0)


    industrial_demand.to_csv(snakemake.output.industrial_demand)



if __name__ == "__main__":

    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils import Dict
        import yaml
        snakemake = Dict()
        snakemake.input = Dict()
        snakemake.input['clustered_pop_layout'] = "resources/pop_layout_elec_s_128.csv"
        snakemake.input['industrial_demand_per_country']="resources/industrial_demand_per_country.csv"
        snakemake.output = Dict()
        snakemake.output['industrial_demand'] = "resources/industrial_demand_elec_s_128.csv"
        with open('config.yaml', encoding='utf8') as f:
            snakemake.config = yaml.safe_load(f)

    build_industrial_demand()
