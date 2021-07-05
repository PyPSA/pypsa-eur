import pandas as pd

rename = {"UK" : "GB", "BH" : "BA"}


def build_biomass_potentials():

    config = snakemake.config['biomass']
    year = config["year"]
    scenario = config["scenario"]

    df = pd.read_excel(snakemake.input.jrc_potentials,
                    "Potentials (PJ)",
                    index_col=[0,1])

    df.rename(columns={"Unnamed: 18": "Municipal waste"}, inplace=True)
    df.drop(columns="Total", inplace=True)
    df.replace("-", 0., inplace=True)

    column = df.iloc[:,0]
    countries = column.where(column.str.isalpha()).pad()
    countries = [rename.get(ct, ct) for ct in countries]
    countries_i = pd.Index(countries, name='country')
    df.set_index(countries_i, append=True, inplace=True)

    df.drop(index='MS', level=0, inplace=True)

    # convert from PJ to MWh
    df = df / 3.6 * 1e6

    df.to_csv(snakemake.output.biomass_potentials_all)

    # solid biomass includes:
    # Primary agricultural residues (MINBIOAGRW1),
    # Forestry energy residue (MINBIOFRSF1),
    # Secondary forestry residues (MINBIOWOOW1),
    # Secondary Forestry residues – sawdust (MINBIOWOO1a)',
    # Forestry residues from landscape care biomass (MINBIOFRSF1a),
    # Municipal waste (MINBIOMUN1)',

    # biogas includes:
    # Manure biomass potential (MINBIOGAS1),
    # Sludge biomass (MINBIOSLU1),

    df = df.loc[year, scenario, :]

    grouper = {v: k for k, vv in config["classes"].items() for v in vv}
    df = df.groupby(grouper, axis=1).sum()

    df.index.name = "MWh/a"

    df.to_csv(snakemake.output.biomass_potentials)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('build_biomass_potentials')


    # This is a hack, to be replaced once snakemake is unicode-conform

    solid_biomass = snakemake.config['biomass']['classes']['solid biomass']
    if 'Secondary Forestry residues sawdust' in solid_biomass:
        solid_biomass.remove('Secondary Forestry residues sawdust')
        solid_biomass.append('Secondary Forestry residues – sawdust')

    build_biomass_potentials()
