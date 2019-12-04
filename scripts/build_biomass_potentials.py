
import pandas as pd

idx = pd.IndexSlice

def build_biomass_potentials():

    #delete empty column C from this sheet first before reading it in
    df = pd.read_excel(snakemake.input.jrc_potentials,
                       "Potentials (PJ)",
                       index_col=[0,1])

    df.rename(columns={"Unnamed: 18":"Municipal waste"},inplace=True)
    df.drop(columns="Total",inplace=True)
    df.replace("-",0.,inplace=True)

    df_dict = {}

    for i in range(36):
        df_dict[df.iloc[i*16,1]] = df.iloc[1+i*16:(i+1)*16].astype(float)

    df_new = pd.concat(df_dict)

    #  solid biomass includes: Primary agricultural residues (MINBIOAGRW1),
    #  Forestry energy residue (MINBIOFRSF1), 
    #  Secondary forestry residues (MINBIOWOOW1),
    #  Secondary Forestry residues – sawdust (MINBIOWOO1a)',
    #  Forestry residues from landscape care biomass (MINBIOFRSF1a), 
    #  Municipal waste (MINBIOMUN1)',
    
    # biogas includes : Manure biomass potential (MINBIOGAS1),
    # Sludge biomass (MINBIOSLU1)
    
    us_type = pd.Series(index=df_new.columns)
    us_type.iloc[0:7] = "not included"
    us_type.iloc[7:8] = "biogas"
    us_type.iloc[8:9] = "solid biomass"
    us_type.iloc[9:11] = "not included"
    us_type.iloc[11:16] = "solid biomass"
    us_type.iloc[16:17] = "biogas"


    #convert from PJ to MWh
    biomass_potentials = df_new.loc[idx[:,snakemake.config['biomass']['year'],snakemake.config['biomass']['scenario']],:].groupby(us_type,axis=1).sum().groupby(level=0).sum().rename({"UK" : "GB", "BH" : "BA"})/3.6*1e6
    biomass_potentials.to_csv(snakemake.output.biomass_potentials)



if __name__ == "__main__":


    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils import Dict
        import yaml
        snakemake = Dict()
        snakemake.input = Dict()
        snakemake.input['jrc_potentials'] = "data/biomass/JRC Biomass Potentials.xlsx"
        snakemake.output = Dict()
        snakemake.output['biomass_potentials'] = 'data/biomass_potentials.csv'
        with open('config.yaml') as f:
            snakemake.config = yaml.load(f)

    build_biomass_potentials()
