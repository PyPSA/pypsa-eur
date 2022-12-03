"""Build industrial energy demand per node."""

import pandas as pd

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_industrial_energy_demand_per_node',
            simpl='',
            clusters=48,
            planning_horizons=2030,
        )

    # import EU ratios df as csv
    fn = snakemake.input.industry_sector_ratios
    industry_sector_ratios = pd.read_csv(fn, index_col=0)

    # material demand per node and industry (kton/a)
    fn = snakemake.input.industrial_production_per_node
    nodal_production = pd.read_csv(fn, index_col=0)

    # final energy consumption per node, sector and carrier
    nodal_dict = {k: s * industry_sector_ratios for k, s in nodal_production.iterrows()}
    nodal_df = pd.concat(nodal_dict, axis=1).T

    # convert GWh to TWh and ktCO2 to MtCO2
    nodal_df *= 0.001

    rename_sectors = {
        'elec': 'electricity',
        'biomass': 'solid biomass',
        'heat': 'low-temperature heat'
    }
    nodal_df.rename(columns=rename_sectors, inplace=True)

    nodal_df.index.set_names(["node", "sector"], inplace=True)
    nodal_df.index.name = "TWh/a (MtCO2/a)"

    fn = snakemake.output.industrial_energy_demand_per_node
    nodal_df.to_csv(fn, float_format='%.2f')
