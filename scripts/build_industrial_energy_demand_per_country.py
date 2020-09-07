
import pandas as pd
import numpy as np


tj_to_ktoe = 0.0238845
ktoe_to_twh = 0.01163

eb_base_dir = "data/eurostat-energy_balances-may_2018_edition"
jrc_base_dir = "data/jrc-idees-2015"

# import EU ratios df as csv
industry_sector_ratios=pd.read_csv(snakemake.input.industry_sector_ratios,
                                   index_col=0)

#material demand per country and industry (kton/a)
countries_production = pd.read_csv(snakemake.input.industrial_production_per_country, index_col=0)

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


eb_names={'NO':'Norway', 'AL':'Albania', 'BA':'Bosnia and Herzegovina',
          'MK':'FYR of Macedonia', 'GE':'Georgia', 'IS':'Iceland',
          'KO':'Kosovo', 'MD':'Moldova', 'ME':'Montenegro', 'RS':'Serbia',
          'UA':'Ukraine', 'TR':'Turkey', }

jrc_names = {"GR" : "EL",
             "GB" : "UK"}

#final energy consumption per country and industry (TWh/a)
countries_df = countries_production.dot(industry_sector_ratios.T)
countries_df*= 0.001 #GWh -> TWh (ktCO2 -> MtCO2)



non_EU = ['NO', 'CH', 'ME', 'MK', 'RS', 'BA', 'AL']


# save current electricity consumption
for country in countries_df.index:
    if country in non_EU:
        if country == 'CH':
            countries_df.loc[country, 'current electricity']=dic_Switzerland['current electricity']*tj_to_ktoe*ktoe_to_twh
        else:
            excel_balances = pd.read_excel('{}/{}.XLSX'.format(eb_base_dir,eb_names[country]),
                                      sheet_name='2016', index_col=1,header=0, skiprows=1 ,squeeze=True)

            countries_df.loc[country, 'current electricity'] = excel_balances.loc['Industry', 'Electricity']*ktoe_to_twh

    else:

        excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(jrc_base_dir,jrc_names.get(country,country)),
                                  sheet_name='Ind_Summary',index_col=0,header=0,squeeze=True) # the summary sheet

        s_out = excel_out.iloc[27:48,-1]
        countries_df.loc[country, 'current electricity'] = s_out['Electricity']*ktoe_to_twh


rename_sectors = {'elec':'electricity',
                  'biomass':'solid biomass',
                  'heat':'low-temperature heat'}

countries_df.rename(columns=rename_sectors,inplace=True)

countries_df.index.name = "TWh/a (MtCO2/a)"

countries_df.to_csv(snakemake.output.industrial_energy_demand_per_country,
                    float_format='%.2f')
