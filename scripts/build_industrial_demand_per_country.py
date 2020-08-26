

#%matplotlib inline
import pandas as pd
import numpy as np



jrc_base_dir = "data/jrc-idees-2015"
eb_base_dir = "data/eurostat-energy_balances-may_2018_edition"



tj_to_ktoe = 0.0238845

ktoe_to_twh = 0.01163

# import EU ratios df as csv
df=pd.read_csv(snakemake.input.industry_sector_ratios, sep=';', index_col=0)





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

rename = {"GR" : "EL",
          "GB" : "UK"}

eu28 = ['FR', 'DE', 'GB', 'IT', 'ES', 'PL', 'SE', 'NL', 'BE', 'FI', 'CZ',
        'DK', 'PT', 'RO', 'AT', 'BG', 'EE', 'GR', 'LV',
        'HU', 'IE', 'SK', 'LT', 'HR', 'LU', 'SI'] + ['CY','MT']


countries = non_EU + [rename.get(eu,eu) for eu in eu28]


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

#final energy consumption per country and industry (TWh/a)
countries_df = pd.DataFrame(index=countries,
                            columns=index,
                            dtype=float)

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
dic_countries={'NO':'Norway', 'AL':'Albania', 'BA':'Bosnia and Herzegovina',
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
    countries_df.loc[country] = 0.
    countries_demand.loc[country] = 0.
    print(country)
    for sector in sectors:
        if country in non_EU:
            if country == 'CH':
                e_country = dic_Switzerland[sector]*tj_to_ktoe
            else:
                # estimate physical output
                #energy consumption in the sector and country
                excel_balances = pd.read_excel('{}/{}.XLSX'.format(eb_base_dir,dic_countries[country]),
                                      sheet_name='2016', index_col=2,header=0, skiprows=1 ,squeeze=True)
                e_country = excel_balances.loc[dic_sec[sector], 'Total all products']

            #energy consumption in the sector and EU28
            excel_sum_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_EU28.xlsx'.format(jrc_base_dir),
                                  sheet_name='Ind_Summary', index_col=0,header=0,squeeze=True) # the summary sheet
            s_sum_out = excel_sum_out.iloc[49:76,-1]
            e_EU28 = s_sum_out[dic_sec_summary[sector]]

            ratio_country_EU28=e_country/e_EU28

            excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_EU28.xlsx'.format(jrc_base_dir),
                                      sheet_name=sub_sheet_name_dict[sector],index_col=0,header=0,squeeze=True) # the summary sheet

            s_out = excel_out.iloc[loc_dic[sector][0]:loc_dic[sector][1],-1]

            for subsector in sect2sub[sector]:
                output = ratio_country_EU28*s_out[out_dic[subsector]]
                countries_demand.loc[country,subsector] = output
                for ind in index:
                    countries_df.loc[country, ind] += float(output*df.loc[ind, subsector]) # kton * MWh = GWh (# kton * tCO2 = ktCO2)

        else:

            # read the input sheets
            excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(jrc_base_dir,country), sheet_name=sub_sheet_name_dict[sector],index_col=0,header=0,squeeze=True) # the summary sheet

            s_out = excel_out.iloc[loc_dic[sector][0]:loc_dic[sector][1],-1]

            for subsector in sect2sub[sector]:
                output = s_out[out_dic[subsector]]
                countries_demand.loc[country,subsector] = output
                for ind in index:
                    countries_df.loc[country, ind] += output*df.loc[ind, subsector] #kton * MWh = GWh (# kton * tCO2 = ktCO2)

countries_df*= 0.001 #GWh -> TWh (ktCO2 -> MtCO2)

# save current electricity consumption
for country in countries:
    if country in non_EU:
        if country == 'CH':
            countries_df.loc[country, 'current electricity']=dic_Switzerland['current electricity']*tj_to_ktoe*ktoe_to_twh
        else:
            excel_balances = pd.read_excel('{}/{}.XLSX'.format(eb_base_dir,dic_countries[country]),
                                      sheet_name='2016', index_col=1,header=0, skiprows=1 ,squeeze=True)

            countries_df.loc[country, 'current electricity'] = excel_balances.loc['Industry', 'Electricity']*ktoe_to_twh

    else:

        excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(jrc_base_dir,country),
                                  sheet_name='Ind_Summary',index_col=0,header=0,squeeze=True) # the summary sheet

        s_out = excel_out.iloc[27:48,-1]
        countries_df.loc[country, 'current electricity'] = s_out['Electricity']*ktoe_to_twh
        print(countries_df.loc[country, 'current electricity'])



# save df as csv
for ind in index:
    countries_df[ind]=countries_df[ind].astype('float')
countries_df = countries_df.round(3)

countries_df.rename(index={value : key for key,value in rename.items()},inplace=True)

rename_sectors = {'elec':'electricity',
                  'biomass':'solid biomass',
                  'heat':'low-temperature heat'}

countries_df.rename(columns=rename_sectors,inplace=True)

countries_df.to_csv(snakemake.output.industrial_energy_demand_per_country,
                    float_format='%.2f')
countries_demand.to_csv(snakemake.output.industrial_demand_per_country,
                        float_format='%.2f')
