
import pandas as pd

# sub-sectors as used in PyPSA-Eur-Sec and listed in JRC-IDEES industry sheets
sub_sectors = {'Iron and steel' : ['Integrated steelworks','Electric arc'],
               'Non-ferrous metals' : ['Alumina production','Aluminium - primary production','Aluminium - secondary production','Other non-ferrous metals'],
               'Chemicals' : ['Basic chemicals', 'Other chemicals', 'Pharmaceutical products etc.', 'Basic chemicals feedstock'],
               'Non-metalic mineral' : ['Cement','Ceramics & other NMM','Glass production'],
               'Printing' : ['Pulp production','Paper production','Printing and media reproduction'],
               'Food' : ['Food, beverages and tobacco'],
               'Transport equipment' : ['Transport Equipment'],
               'Machinery equipment' : ['Machinery Equipment'],
               'Textiles and leather' : ['Textiles and leather'],
               'Wood and wood products' : ['Wood and wood products'],
               'Other Industrial Sectors' : ['Other Industrial Sectors'],
}


# name in JRC-IDEES Energy Balances
eb_sheet_name = {'Integrated steelworks' : 'cisb',
                 'Electric arc' : 'cise',
                 'Alumina production' : 'cnfa',
                 'Aluminium - primary production' : 'cnfp',
                 'Aluminium - secondary production' : 'cnfs',
                 'Other non-ferrous metals' : 'cnfo',
                 'Basic chemicals' : 'cbch',
                 'Other chemicals' : 'coch',
                 'Pharmaceutical products etc.' : 'cpha',
                 'Basic chemicals feedstock' : 'cpch',
                 'Cement' : 'ccem',
                 'Ceramics & other NMM' : 'ccer',
                 'Glass production' : 'cgla',
                 'Pulp production' : 'cpul',
                 'Paper production' : 'cpap',
                 'Printing and media reproduction' : 'cprp',
                 'Food, beverages and tobacco' : 'cfbt',
                 'Transport Equipment' : 'ctre',
                 'Machinery Equipment' : 'cmae',
                 'Textiles and leather' : 'ctel',
                 'Wood and wood products' : 'cwwp',
                 'Mining and quarrying' : 'cmiq',
                 'Construction' : 'ccon',
                 'Non-specified': 'cnsi',
}



fuels = {'all' : ['All Products'],
         'solid' : ['Solid Fuels'],
         'liquid' : ['Total petroleum products (without biofuels)'],
         'gas' : ['Gases'],
         'heat' : ['Nuclear heat','Derived heat'],
         'biomass' : ['Biomass and Renewable wastes'],
         'waste' : ['Wastes (non-renewable)'],
         'electricity' : ['Electricity'],
}

ktoe_to_twh = 0.011630

eu28 = ['FR', 'DE', 'GB', 'IT', 'ES', 'PL', 'SE', 'NL', 'BE', 'FI',
        'DK', 'PT', 'RO', 'AT', 'BG', 'EE', 'GR', 'LV', 'CZ',
        'HU', 'IE', 'SK', 'LT', 'HR', 'LU', 'SI', 'CY', 'MT']

jrc_names = {"GR" : "EL",
             "GB" : "UK"}

year = 2015
summaries = {}

#for some reason the Energy Balances list Other Industrial Sectors separately
ois_subs = ['Mining and quarrying','Construction','Non-specified']



for ct in eu28:
    print(ct)
    filename = 'data/jrc-idees-2015/JRC-IDEES-2015_EnergyBalance_{}.xlsx'.format(jrc_names.get(ct,ct))

    summary = pd.DataFrame(index=list(fuels.keys()) + ['other'])

    for sector in sub_sectors:
        if sector == 'Other Industrial Sectors':
            subs = ois_subs
        else:
            subs = sub_sectors[sector]

        for sub in subs:
            df = pd.read_excel(filename,
                               sheet_name=eb_sheet_name[sub],
                               index_col=0)

            s = df[year].astype(float)

            for fuel in fuels:
                summary.at[fuel,sub] = s[fuels[fuel]].sum()
                summary.at['other',sub] = summary.at['all',sub] - summary.loc[summary.index^['all','other'],sub].sum()

    summary['Other Industrial Sectors'] = summary[ois_subs].sum(axis=1)
    summary.drop(columns=ois_subs,inplace=True)

    summaries[ct] = summary*ktoe_to_twh

final_summary = pd.concat(summaries,axis=1)

final_summary.index.name = 'TWh/a'

final_summary.to_csv(snakemake.output.industrial_energy_demand_per_country_today)
