# coding: utf-8 

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

    #convert from PJ to MWh
    df_new = pd.concat(df_dict).rename({"UK" : "GB", "BH" : "BA"})/3.6*1e6
    df_new.index.name = "MWh/a"
    df_new.to_csv(snakemake.output.biomass_potentials_all)

    #  solid biomass includes: Primary agricultural residues (MINBIOAGRW1),
    #  Forestry energy residue (MINBIOFRSF1),
    #  Secondary forestry residues (MINBIOWOOW1),
    #  Secondary Forestry residues – sawdust (MINBIOWOO1a)',
    #  Forestry residues from landscape care biomass (MINBIOFRSF1a),
    #  Municipal waste (MINBIOMUN1)',

    # biogas includes : Manure biomass potential (MINBIOGAS1),
    # Sludge biomass (MINBIOSLU1)

    us_type = pd.Series("", df_new.columns)

    for k,v in snakemake.config['biomass']['classes'].items():
        us_type.loc[v] = k

    biomass_potentials = df_new.swaplevel(0,2).loc[snakemake.config['biomass']['scenario'],snakemake.config['biomass']['year']].groupby(us_type,axis=1).sum()
    biomass_potentials.index.name = "MWh/a"
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
        snakemake.output['biomass_potentials_all']='resources/biomass_potentials_all.csv'
        with open('config.yaml', encoding='utf8') as f:
            snakemake.config = yaml.safe_load(f)
    

    # This is a hack, to be replaced once snakemake is unicode-conform

    if 'Secondary Forestry residues sawdust' in snakemake.config['biomass']['classes']['solid biomass']:
        snakemake.config['biomass']['classes']['solid biomass'].remove('Secondary Forestry residues sawdust')
        snakemake.config['biomass']['classes']['solid biomass'].append('Secondary Forestry residues – sawdust')
    
    build_biomass_potentials()
