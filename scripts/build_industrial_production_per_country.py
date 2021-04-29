
import pandas as pd
import numpy as np


tj_to_ktoe = 0.0238845
ktoe_to_twh = 0.01163

jrc_base_dir = "data/jrc-idees-2015"
eb_base_dir = "data/eurostat-energy_balances-may_2018_edition"

# year for which data is retrieved
raw_year = 2015
year = raw_year-2016

sub_sheet_name_dict = { 'Iron and steel':'ISI',
                        'Chemicals Industry':'CHI',
                        'Non-metallic mineral products': 'NMM',
                        'Pulp, paper and printing': 'PPA',
                        'Food, beverages and tobacco': 'FBT',
                        'Non Ferrous Metals' : 'NFM',
                        'Transport Equipment': 'TRE',
                        'Machinery Equipment': 'MAE',
                        'Textiles and leather':'TEL',
                        'Wood and wood products': 'WWP',
                        'Other Industrial Sectors': 'OIS'}

index = ['elec','biomass','methane','hydrogen','heat','naphtha','process emission','process emission from feedstock']

non_EU = ['NO', 'CH', 'ME', 'MK', 'RS', 'BA', 'AL']

jrc_names = {"GR" : "EL",
             "GB" : "UK"}

eu28 = ['FR', 'DE', 'GB', 'IT', 'ES', 'PL', 'SE', 'NL', 'BE', 'FI',
        'DK', 'PT', 'RO', 'AT', 'BG', 'EE', 'GR', 'LV', 'CZ',
        'HU', 'IE', 'SK', 'LT', 'HR', 'LU', 'SI', 'CY', 'MT']


countries = non_EU + eu28


sectors = ['Iron and steel','Chemicals Industry','Non-metallic mineral products',
           'Pulp, paper and printing', 'Food, beverages and tobacco', 'Non Ferrous Metals',
           'Transport Equipment', 'Machinery Equipment', 'Textiles and leather',
           'Wood and wood products', 'Other Industrial Sectors']

sect2sub = {'Iron and steel':['Electric arc','Integrated steelworks'],
            'Chemicals Industry': ['Basic chemicals', 'Other chemicals', 'Pharmaceutical products etc.'],
            'Non-metallic mineral products': ['Cement','Ceramics & other NMM','Glass production'],
            'Pulp, paper and printing': ['Pulp production','Paper production','Printing and media reproduction'],
            'Food, beverages and tobacco': ['Food, beverages and tobacco'],
            'Non Ferrous Metals': ['Alumina production', 'Aluminium - primary production', 'Aluminium - secondary production', 'Other non-ferrous metals'],
            'Transport Equipment': ['Transport Equipment'],
            'Machinery Equipment': ['Machinery Equipment'],
            'Textiles and leather': ['Textiles and leather'],
            'Wood and wood products' :['Wood and wood products'],
            'Other Industrial Sectors':['Other Industrial Sectors']}

subsectors = [ss for s in sectors for ss in sect2sub[s]]

#material demand per country and industry (kton/a)
countries_demand = pd.DataFrame(index=countries,
                                columns=subsectors,
                                dtype=float)


out_dic ={'Electric arc': 'Electric arc',
          'Integrated steelworks': 'Integrated steelworks',
          'Basic chemicals': 'Basic chemicals (kt ethylene eq.)',
          'Other chemicals':'Other chemicals (kt ethylene eq.)',
          'Pharmaceutical products etc.':'Pharmaceutical products etc. (kt ethylene eq.)',
          'Cement':'Cement (kt)',
          'Ceramics & other NMM':'Ceramics & other NMM (kt bricks eq.)',
          'Glass production':'Glass production  (kt)',
          'Pulp production':'Pulp production (kt)',
          'Paper production':'Paper production  (kt)',
          'Printing and media reproduction':'Printing and media reproduction (kt paper eq.)',
          'Food, beverages and tobacco': 'Physical output (index)',
          'Alumina production':'Alumina production (kt)',
          'Aluminium - primary production': 'Aluminium - primary production',
          'Aluminium - secondary production': 'Aluminium - secondary production',
          'Other non-ferrous metals' : 'Other non-ferrous metals (kt lead eq.)',
          'Transport Equipment': 'Physical output (index)',
          'Machinery Equipment': 'Physical output (index)',
          'Textiles and leather':  'Physical output (index)',
          'Wood and wood products': 'Physical output (index)',
          'Other Industrial Sectors': 'Physical output (index)'}

loc_dic={'Iron and steel':[5,8],
         'Chemicals Industry': [7,11],
         'Non-metallic mineral products': [6,10],
         'Pulp, paper and printing': [7,11],
         'Food, beverages and tobacco': [2,6],
         'Non Ferrous Metals': [9,14],
         'Transport Equipment': [3,5],
         'Machinery Equipment': [3,5],
         'Textiles and leather':  [3,5],
         'Wood and wood products': [3,5],
         'Other Industrial Sectors': [3,5]}

# In the summary sheet (IDEES database) some names include a white space
dic_sec_summary = {'Iron and steel': 'Iron and steel',
                   'Chemicals Industry': 'Chemicals Industry',
                   'Non-metallic mineral products': 'Non-metallic mineral products',
                   'Pulp, paper and printing': 'Pulp, paper and printing',
                   'Food, beverages and tobacco': ' Food, beverages and tobacco',
                   'Non Ferrous Metals': 'Non Ferrous Metals',
                   'Transport Equipment': ' Transport Equipment',
                   'Machinery Equipment': ' Machinery Equipment',
                   'Textiles and leather': ' Textiles and leather',
                   'Wood and wood products': ' Wood and wood products',
                   'Other Industrial Sectors': ' Other Industrial Sectors'}

#countries=['CH']
eb_names={'NO':'Norway', 'AL':'Albania', 'BA':'Bosnia and Herzegovina',
          'MK':'FYR of Macedonia', 'GE':'Georgia', 'IS':'Iceland',
          'KO':'Kosovo', 'MD':'Moldova', 'ME':'Montenegro', 'RS':'Serbia',
          'UA':'Ukraine', 'TR':'Turkey', }

dic_sec ={'Iron and steel':'Iron & steel industry',
          'Chemicals Industry': 'Chemical and Petrochemical industry',
          'Non-metallic mineral products': 'Non-ferrous metal industry',
          'Pulp, paper and printing': 'Paper, Pulp and Print',
          'Food, beverages and tobacco': 'Food and Tabacco',
          'Non Ferrous Metals': 'Non-metallic Minerals (Glass, pottery & building mat. Industry)',
          'Transport Equipment': 'Transport Equipment',
          'Machinery Equipment': 'Machinery',
          'Textiles and leather': 'Textile and Leather',
          'Wood and wood products': 'Wood and Wood Products',
          'Other Industrial Sectors': 'Non-specified (Industry)'}
          # Mining and Quarrying, Construction

#Annual energy consumption in Switzerland by sector in 2015 (in TJ)
#From: Energieverbrauch in der Industrie und im Dienstleistungssektor, Der Bundesrat
#http://www.bfe.admin.ch/themen/00526/00541/00543/index.html?lang=de&dossier_id=00775

dic_Switzerland ={'Iron and steel': 7889.,
          'Chemicals Industry': 26871.,
          'Non-metallic mineral products': 15513.+3820.,
          'Pulp, paper and printing': 12004.,
          'Food, beverages and tobacco': 17728.,
          'Non Ferrous Metals': 3037.,
          'Transport Equipment': 14993.,
          'Machinery Equipment': 4724.,
          'Textiles and leather': 1742.,
          'Wood and wood products': 0.,
          'Other Industrial Sectors': 10825.,
          'current electricity': 53760.}

dic_sec_position={}
for country in countries:
    countries_demand.loc[country] = 0.
    print(country)
    for sector in sectors:
        if country in non_EU:
            if country == 'CH':
                e_country = dic_Switzerland[sector]*tj_to_ktoe
            else:
                # estimate physical output
                #energy consumption in the sector and country
                excel_balances = pd.read_excel('{}/{}.XLSX'.format(eb_base_dir,eb_names[country]),
                                      sheet_name='2016', index_col=2,header=0, skiprows=1 ,squeeze=True)
                e_country = excel_balances.loc[dic_sec[sector], 'Total all products']

            #energy consumption in the sector and EU28
            excel_sum_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_EU28.xlsx'.format(jrc_base_dir),
                                  sheet_name='Ind_Summary', index_col=0,header=0,squeeze=True) # the summary sheet
            s_sum_out = excel_sum_out.iloc[49:76,year]
            e_EU28 = s_sum_out[dic_sec_summary[sector]]

            ratio_country_EU28=e_country/e_EU28

            excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_EU28.xlsx'.format(jrc_base_dir),
                                      sheet_name=sub_sheet_name_dict[sector],index_col=0,header=0,squeeze=True) # the summary sheet

            s_out = excel_out.iloc[loc_dic[sector][0]:loc_dic[sector][1],year]

            for subsector in sect2sub[sector]:
                countries_demand.loc[country,subsector] = ratio_country_EU28*s_out[out_dic[subsector]]

        else:

            # read the input sheets
            excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(jrc_base_dir,jrc_names.get(country,country)), sheet_name=sub_sheet_name_dict[sector],index_col=0,header=0,squeeze=True) # the summary sheet

            s_out = excel_out.iloc[loc_dic[sector][0]:loc_dic[sector][1],year]

            for subsector in sect2sub[sector]:
                countries_demand.loc[country,subsector] = s_out[out_dic[subsector]]


#include ammonia demand separately and remove ammonia from basic chemicals

ammonia = pd.read_csv(snakemake.input.ammonia_production,
                      index_col=0)

there = ammonia.index.intersection(countries_demand.index)
missing = countries_demand.index.symmetric_difference(there)

print("Following countries have no ammonia demand:", missing)

countries_demand.insert(2,"Ammonia",0.)

countries_demand.loc[there,"Ammonia"] = ammonia.loc[there, str(raw_year)]

countries_demand["Basic chemicals"] -= countries_demand["Ammonia"]

#EE, HR and LT got negative demand through subtraction - poor data
countries_demand.loc[countries_demand["Basic chemicals"] < 0.,"Basic chemicals"] = 0.

countries_demand.rename(columns={"Basic chemicals" : "Basic chemicals (without ammonia)"},
                        inplace=True)

countries_demand.index.name = "kton/a"

countries_demand.to_csv(snakemake.output.industrial_production_per_country,
                        float_format='%.2f')
