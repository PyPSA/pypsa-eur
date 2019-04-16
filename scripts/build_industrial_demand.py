
import pandas as pd

idx = pd.IndexSlice

def build_industrial_demand():

    population = pd.read_csv(snakemake.input.clustered_pop_layout,
                             index_col=0)

    totals = pd.Series(data=[1100.,1814.,586.,400.,580.,186.],
                       index=["industry new electricity","industry process heat",
                              "naptha feedstock","shipping H2","aviation kerosene","process emissions"])

    industrial_demand = pd.DataFrame({i : population["total"]*totals[i]*1e6/population["total"].sum() for i in totals.index })

    industrial_demand.to_csv(snakemake.output.industrial_demand)



if __name__ == "__main__":

    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils import Dict
        import yaml
        snakemake = Dict()
        snakemake.input = Dict()
        snakemake.input['clustered_pop_layout'] = "resources/pop_layout_elec_s_128.csv"
        snakemake.output = Dict()
        snakemake.output['industrial_demand'] = "resources/industrial_demand_elec_s_128.csv"
        with open('config.yaml') as f:
            snakemake.config = yaml.load(f)

    build_industrial_demand()
