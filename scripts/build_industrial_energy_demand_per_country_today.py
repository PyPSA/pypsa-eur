"""Build industrial energy demand per country."""

import pandas as pd
import multiprocessing as mp
from tqdm import tqdm

ktoe_to_twh = 0.011630

# name in JRC-IDEES Energy Balances
sector_sheets = {'Integrated steelworks': 'cisb',
                 'Electric arc': 'cise',
                 'Alumina production': 'cnfa',
                 'Aluminium - primary production': 'cnfp',
                 'Aluminium - secondary production': 'cnfs',
                 'Other non-ferrous metals': 'cnfo',
                 'Basic chemicals': 'cbch',
                 'Other chemicals': 'coch',
                 'Pharmaceutical products etc.': 'cpha',
                 'Basic chemicals feedstock': 'cpch',
                 'Cement': 'ccem',
                 'Ceramics & other NMM': 'ccer',
                 'Glass production': 'cgla',
                 'Pulp production': 'cpul',
                 'Paper production': 'cpap',
                 'Printing and media reproduction': 'cprp',
                 'Food, beverages and tobacco': 'cfbt',
                 'Transport Equipment': 'ctre',
                 'Machinery Equipment': 'cmae',
                 'Textiles and leather': 'ctel',
                 'Wood and wood products': 'cwwp',
                 'Mining and quarrying': 'cmiq',
                 'Construction': 'ccon',
                 'Non-specified': 'cnsi',
                 }


fuels = {'All Products': 'all',
         'Solid Fuels': 'solid',
         'Total petroleum products (without biofuels)': 'liquid',
         'Gases': 'gas',
         'Nuclear heat': 'heat',
         'Derived heat': 'heat',
         'Biomass and Renewable wastes': 'biomass',
         'Wastes (non-renewable)': 'waste',
         'Electricity': 'electricity'
        }

eu28 = ['FR', 'DE', 'GB', 'IT', 'ES', 'PL', 'SE', 'NL', 'BE', 'FI',
        'DK', 'PT', 'RO', 'AT', 'BG', 'EE', 'GR', 'LV', 'CZ',
        'HU', 'IE', 'SK', 'LT', 'HR', 'LU', 'SI', 'CY', 'MT']

jrc_names = {"GR": "EL", "GB": "UK"}


def industrial_energy_demand_per_country(country):

    jrc_dir = snakemake.input.jrc
    jrc_country = jrc_names.get(country, country)
    fn = f'{jrc_dir}/JRC-IDEES-2015_EnergyBalance_{jrc_country}.xlsx'

    sheets = list(sector_sheets.values())
    df_dict = pd.read_excel(fn, sheet_name=sheets, index_col=0)

    def get_subsector_data(sheet):

        df = df_dict[sheet][year].groupby(fuels).sum()

        df['other'] = df['all'] - df.loc[df.index != 'all'].sum()

        return df

    df = pd.concat({sub: get_subsector_data(sheet)
                    for sub, sheet in sector_sheets.items()}, axis=1)

    sel = ['Mining and quarrying', 'Construction', 'Non-specified']
    df['Other Industrial Sectors'] = df[sel].sum(axis=1)
    df['Basic chemicals'] += df['Basic chemicals feedstock']

    df.drop(columns=sel+['Basic chemicals feedstock'], index='all', inplace=True)

    df *= ktoe_to_twh

    return df


def add_ammonia_energy_demand(demand):

    # MtNH3/a
    fn = snakemake.input.ammonia_production
    ammonia = pd.read_csv(fn, index_col=0)[str(year)] / 1e3

    def ammonia_by_fuel(x):

        fuels = {'gas': config['MWh_CH4_per_tNH3_SMR'],
                 'electricity': config['MWh_elec_per_tNH3_SMR']}

        return pd.Series({k: x*v for k,v in fuels.items()})

    ammonia = ammonia.apply(ammonia_by_fuel).T

    demand['Ammonia'] = ammonia.unstack().reindex(index=demand.index, fill_value=0.)

    demand['Basic chemicals (without ammonia)'] = demand["Basic chemicals"] - demand["Ammonia"]

    demand['Basic chemicals (without ammonia)'].clip(lower=0, inplace=True)

    demand.drop(columns='Basic chemicals', inplace=True)

    return demand


def add_non_eu28_industrial_energy_demand(demand):

    # output in MtMaterial/a
    fn = snakemake.input.industrial_production_per_country
    production = pd.read_csv(fn, index_col=0) / 1e3

    #recombine HVC, Chlorine and Methanol to Basic chemicals (without ammonia)
    chemicals = ["HVC", "Chlorine", "Methanol"]
    production["Basic chemicals (without ammonia)"] = production[chemicals].sum(axis=1)
    production.drop(columns=chemicals, inplace=True)

    eu28_production = production.loc[eu28].sum()
    eu28_energy = demand.groupby(level=1).sum()
    eu28_averages = eu28_energy / eu28_production

    non_eu28 = production.index.symmetric_difference(eu28)

    demand_non_eu28 = pd.concat({k: v * eu28_averages
        for k, v in production.loc[non_eu28].iterrows()})

    return pd.concat([demand, demand_non_eu28])


def industrial_energy_demand(countries):

    nprocesses = snakemake.threads
    func = industrial_energy_demand_per_country
    tqdm_kwargs = dict(ascii=False, unit=' country', total=len(countries),
                       desc="Build industrial energy demand")
    with mp.Pool(processes=nprocesses) as pool:
        demand_l = list(tqdm(pool.imap(func, countries), **tqdm_kwargs))

    demand = pd.concat(demand_l, keys=countries)

    return demand


if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('build_industrial_energy_demand_per_country_today')

    config = snakemake.config['industry']
    year = config.get('reference_year', 2015)

    demand = industrial_energy_demand(eu28)

    demand = add_ammonia_energy_demand(demand)

    demand = add_non_eu28_industrial_energy_demand(demand)

    # for format compatibility
    demand = demand.stack(dropna=False).unstack(level=[0,2])

    # style and annotation
    demand.index.name = 'TWh/a'
    demand.sort_index(axis=1, inplace=True)

    fn = snakemake.output.industrial_energy_demand_per_country_today
    demand.to_csv(fn)
