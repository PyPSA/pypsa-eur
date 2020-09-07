

import pandas as pd
import numpy as np

base_dir = "data/jrc-idees-2015"

# year for which data is retrieved
raw_year = 2015
year = raw_year-2016

conv_factor=11.630 #GWh/ktoe OR MWh/toe

country = 'EU28'


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

index = ['elec','coal','coke','biomass','methane','hydrogen','heat','naphtha','process emission','process emission from feedstock']

df = pd.DataFrame(index=index)


## Iron and steel
#
#> There are two different approaches to produce iron and steel: i.e., integrated steelworks and electric arc.
#
#> Electric arc approach has higher efficiency and relies more on electricity.
#
#> We assume that integrated steelworks will be replaced by electric arc entirely.

sector = 'Iron and steel'

# read the input sheets
excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True) # the summary sheet

excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True) # the final energy consumption sheet

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True) # the used energy sheet

excel_emi = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_emi',
                      index_col=0,header=0,squeeze=True) # the emission sheet

### Electric arc

sector = 'Electric arc'

df[sector] = 0.

# read the corresponding lines
s_fec = excel_fec.iloc[51:57,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec',sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat',sector] += s_fec['Low enthalpy heat']

#### Steel: Smelters

subsector = 'Steel: Smelters'

# read the corresponding lines
s_fec = excel_fec.iloc[61:67,year]

s_ued = excel_ued.iloc[61:67,year]

assert s_fec.index[0] == subsector

# Efficiency changes due to transforming all the smelters into methane
eff_met=s_ued['Natural gas (incl. biogas)']/s_fec['Natural gas (incl. biogas)']

df.loc['methane', sector] += s_ued[subsector]/eff_met

#### Steel: Electric arc

subsector = 'Steel: Electric arc'

# read the corresponding lines
s_fec = excel_fec.iloc[67:68,year]

assert s_fec.index[0] == subsector

# only electricity
df.loc['elec',sector] += s_fec[subsector]

#### Steel: Furnaces, Refining and Rolling
#> assume fully electrified
#
#> other processes are scaled by the used energy

subsector = 'Steel: Furnaces, Refining and Rolling'

# read the corresponding lines
s_fec = excel_fec.iloc[68:75,year]

s_ued = excel_ued.iloc[68:75,year]

assert s_fec.index[0] == subsector

# this process can be electrified
eff = s_ued['Steel: Furnaces, Refining and Rolling - Electric']/s_fec['Steel: Furnaces, Refining and Rolling - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff

#### Steel: Products finishing
#> assume fully electrified

subsector = 'Steel: Products finishing'

# read the corresponding lines
s_fec = excel_fec.iloc[75:92,year]

s_ued = excel_ued.iloc[75:92,year]

assert s_fec.index[0] == subsector

# this process can be electrified
eff = s_ued['Steel: Products finishing - Electric']/s_fec['Steel: Products finishing - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff

#### Process emissions (per physical output)

s_emi = excel_emi.iloc[51:93,year]

assert s_emi.index[0] == sector

s_out = excel_out.iloc[7:8,year]

assert sector in str(s_out.index)

df.loc['process emission',sector] = s_emi['Process emissions']/s_out[sector] # unit tCO2/t material

# final energy consumption per t
df.loc[['elec','heat','methane'],sector] = df.loc[['elec','heat','methane'],sector]*conv_factor/s_out[sector] # unit MWh/t material

### For primary route: DRI with H2 + EAF

df['DRI + Electric arc'] = df['Electric arc']

# adding the Hydrogen necessary for the Direct Reduction of Iron. consumption 1.7 MWh H2 /ton steel
df.loc['hydrogen', 'DRI + Electric arc'] = snakemake.config["industry"]["H2_DRI"]
# add electricity consumption in DRI shaft (0.322 MWh/tSl)
df.loc['elec', 'DRI + Electric arc'] += snakemake.config["industry"]["elec_DRI"]


### Integrated steelworks (could be used in combination with CCS)
### Assume existing fuels are kept, except for furnaces, refining, rolling, finishing
### Ignore 'derived gases' since these are top gases from furnaces

sector = 'Integrated steelworks'

df['Integrated steelworks']= 0.

# read the corresponding lines
s_fec = excel_fec.iloc[3:9,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec',sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat',sector] += s_fec['Low enthalpy heat']


#### Steel: Sinter/Pellet making

subsector = 'Steel: Sinter/Pellet making'

# read the corresponding lines
s_fec = excel_fec.iloc[13:19,year]

s_ued = excel_ued.iloc[13:19,year]

assert s_fec.index[0] == subsector

df.loc['elec',sector] += s_fec['Electricity']
df.loc['methane',sector] += s_fec['Natural gas (incl. biogas)']
df.loc['methane',sector] += s_fec['Residual fuel oil']
df.loc['coal',sector] += s_fec['Solids']


#### Steel: Blast / Basic Oxygen Furnace

subsector = 'Steel: Blast /Basic oxygen furnace'

# read the corresponding lines
s_fec = excel_fec.iloc[19:25,year]

s_ued = excel_ued.iloc[19:25,year]

assert s_fec.index[0] == subsector

df.loc['methane',sector] += s_fec['Natural gas (incl. biogas)']
df.loc['methane',sector] += s_fec['Residual fuel oil']
df.loc['coal',sector] += s_fec['Solids']
df.loc['coke',sector] += s_fec['Coke']


#### Steel: Furnaces, Refining and Rolling
#> assume fully electrified
#
#> other processes are scaled by the used energy

subsector = 'Steel: Furnaces, Refining and Rolling'

# read the corresponding lines
s_fec = excel_fec.iloc[25:32,year]

s_ued = excel_ued.iloc[25:32,year]

assert s_fec.index[0] == subsector

# this process can be electrified
eff = s_ued['Steel: Furnaces, Refining and Rolling - Electric']/s_fec['Steel: Furnaces, Refining and Rolling - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff

#### Steel: Products finishing
#> assume fully electrified

subsector = 'Steel: Products finishing'

# read the corresponding lines
s_fec = excel_fec.iloc[32:49,year]

s_ued = excel_ued.iloc[32:49,year]

assert s_fec.index[0] == subsector

# this process can be electrified
eff = s_ued['Steel: Products finishing - Electric']/s_fec['Steel: Products finishing - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff


#### Process emissions (per physical output)

s_emi = excel_emi.iloc[3:50,year]

assert s_emi.index[0] == sector

s_out = excel_out.iloc[6:7,year]

assert sector in str(s_out.index)

df.loc['process emission',sector] = s_emi['Process emissions']/s_out[sector] # unit tCO2/t material

# final energy consumption per t
df.loc[['elec','heat','methane','coke','coal'],sector] = df.loc[['elec','heat','methane','coke','coal'],sector]*conv_factor/s_out[sector] # unit MWh/t material



## Chemicals Industry

sector = 'Chemicals Industry'

# read the input sheets
excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True) # the summary sheet

excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True) # the final energy consumption sheet

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True) # the used energy sheet

excel_emi = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_emi',
                      index_col=0,header=0,squeeze=True) # the emission sheet

### Basic chemicals

## Ammonia is separated afterwards

sector = 'Basic chemicals'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:9,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec',sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat',sector] += s_fec['Low enthalpy heat']

#### Chemicals: Feedstock (energy used as raw material)
#> There are Solids, Refinery gas, LPG, Diesel oil, Residual fuel oil, Other liquids, Naphtha, Natural gas for feedstock.
#
#>  Naphta represents 47%, methane 17%. LPG (18%) solids, refinery gas, diesel oil, residual fuel oils and other liquids are asimilated to Naphtha


subsector = 'Chemicals: Feedstock (energy used as raw material)'

# read the corresponding lines
s_fec = excel_fec.iloc[13:22,year]

assert s_fec.index[0] == subsector

# naphtha
df.loc['naphtha',sector] += s_fec['Naphtha']

# natural gas
df.loc['methane',sector] += s_fec['Natural gas']

# LPG and other feedstock materials are assimilated to naphtha since they will be produced trough Fischer-Tropsh process
df.loc['naphtha',sector] += (s_fec['Solids'] + s_fec['Refinery gas'] + s_fec['LPG'] + s_fec['Diesel oil']
                            + s_fec['Residual fuel oil'] + s_fec['Other liquids'])

#### Chemicals: Steam processing
#> All the final energy consumption in the Steam processing is converted to methane, since we need >1000 C temperatures here.
#
#> The current efficiency of methane is assumed in the conversion.

subsector = 'Chemicals: Steam processing'

# read the corresponding lines
s_fec = excel_fec.iloc[22:33,year]

s_ued = excel_ued.iloc[22:33,year]

assert s_fec.index[0] == subsector

# efficiency of natural gas
eff_ch4 = s_ued['Natural gas (incl. biogas)']/s_fec['Natural gas (incl. biogas)']

# replace all fec by methane
df.loc['methane',sector] += s_ued[subsector]/eff_ch4

#### Chemicals: Furnaces
#> assume fully electrified

subsector = 'Chemicals: Furnaces'

# read the corresponding lines
s_fec = excel_fec.iloc[33:41,year]

s_ued = excel_ued.iloc[33:41,year]

assert s_fec.index[0] == subsector

#efficiency of electrification
eff_elec = s_ued['Chemicals: Furnaces - Electric']/s_fec['Chemicals: Furnaces - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff_elec

#### Chemicals: Process cooling
#> assume fully electrified

subsector = 'Chemicals: Process cooling'

# read the corresponding lines
s_fec = excel_fec.iloc[41:55,year]

s_ued = excel_ued.iloc[41:55,year]

assert s_fec.index[0] == subsector

eff_elec = s_ued['Chemicals: Process cooling - Electric']/s_fec['Chemicals: Process cooling - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff_elec

#### Chemicals: Generic electric process

subsector = 'Chemicals: Generic electric process'

# read the corresponding lines
s_fec = excel_fec.iloc[55:56,year]

assert s_fec.index[0] == subsector

df.loc['elec',sector] += s_fec[subsector]

#### Process emissions

s_emi = excel_emi.iloc[3:57,year]

assert s_emi.index[0] == sector


## Correct everything by subtracting 2015's ammonia demand and putting in ammonia demand for H2 and electricity separately

s_out = excel_out.iloc[8:9,year]

assert sector in str(s_out.index)

ammonia = pd.read_csv(snakemake.input.ammonia_production,
                      index_col=0)

eu28 = ['FR', 'DE', 'GB', 'IT', 'ES', 'PL', 'SE', 'NL', 'BE', 'FI',
        'DK', 'PT', 'RO', 'AT', 'BG', 'EE', 'GR', 'LV', 'CZ',
        'HU', 'IE', 'SK', 'LT', 'HR', 'LU', 'SI', 'CY', 'MT']

#ktNH3/a
total_ammonia = ammonia.loc[ammonia.index.intersection(eu28),str(raw_year)].sum()

s_out -= total_ammonia

df.loc['process emission',sector] += (s_emi['Process emissions'] - snakemake.config["industry"]['petrochemical_process_emissions']*1e3 - snakemake.config["industry"]['NH3_process_emissions']*1e3)/s_out.values # unit tCO2/t material

#these are emissions originating from feedstock, i.e. could be non-fossil origin
df.loc['process emission from feedstock',sector] += (snakemake.config["industry"]['petrochemical_process_emissions']*1e3)/s_out.values # unit tCO2/t material


# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']

#convert from ktoe/a to GWh/a
df.loc[sources,sector] *= conv_factor

df.loc['methane',sector] -= total_ammonia*snakemake.config['industry']['MWh_CH4_per_tNH3_SMR']
df.loc['elec',sector] -= total_ammonia*snakemake.config['industry']['MWh_elec_per_tNH3_SMR']

df.loc[sources,sector] = df.loc[sources,sector]/s_out.values # unit MWh/t material

df.rename(columns={sector : sector + " (without ammonia)"},
          inplace=True)

sector = 'Ammonia'

df[sector] = 0.

df.loc['hydrogen',sector] = snakemake.config['industry']['MWh_H2_per_tNH3_electrolysis']
df.loc['elec',sector] = snakemake.config['industry']['MWh_elec_per_tNH3_electrolysis']


### Other chemicals

sector = 'Other chemicals'

df[sector] = 0
# read the corresponding lines
s_fec = excel_fec.iloc[58:64,year]

# check the position
assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec',sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat',sector] += s_fec['Low enthalpy heat']

#### Chemicals: High enthalpy heat  processing
#> assume fully electrified

subsector = 'Chemicals: High enthalpy heat  processing'

# read the corresponding lines
s_fec = excel_fec.iloc[68:81,year]

s_ued = excel_ued.iloc[68:81,year]

assert s_fec.index[0] == subsector

eff_elec = s_ued['High enthalpy heat  processing - Electric (microwave)']/s_fec['High enthalpy heat  processing - Electric (microwave)']

df.loc['elec',sector] += s_ued[subsector]/eff_elec

#### Chemicals: Furnaces
#> assume fully electrified

subsector = 'Chemicals: Furnaces'

# read the corresponding lines
s_fec = excel_fec.iloc[81:89,year]

s_ued = excel_ued.iloc[81:89,year]

assert s_fec.index[0] == subsector

eff_elec = s_ued['Chemicals: Furnaces - Electric']/s_fec['Chemicals: Furnaces - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff_elec

#### Chemicals: Process cooling
#> assume fully electrified

subsector = 'Chemicals: Process cooling'

# read the corresponding lines
s_fec = excel_fec.iloc[89:103,year]

s_ued = excel_ued.iloc[89:103,year]

assert s_fec.index[0] == subsector

eff = s_ued['Chemicals: Process cooling - Electric']/s_fec['Chemicals: Process cooling - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff

#### Chemicals: Generic electric process

subsector = 'Chemicals: Generic electric process'

# read the corresponding lines
s_fec = excel_fec.iloc[103:104,year]

assert s_fec.index[0] == subsector

df.loc['elec',sector] += s_fec[subsector]

#### Process emissions

s_emi = excel_emi.iloc[58:105,year]

assert s_emi.index[0] == sector

s_out = excel_out.iloc[9:10,year]

assert sector in str(s_out.index)

df.loc['process emission',sector] += s_emi['Process emissions']/s_out.values # unit tCO2/t material

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']

df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out.values # unit MWh/t material
# 1 ktoe = 11630 MWh

### Pharmaceutical products etc.

sector = 'Pharmaceutical products etc.'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[106:112,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec',sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat',sector] += s_fec['Low enthalpy heat']

#### Chemicals: High enthalpy heat  processing
#> assume fully electrified

subsector = 'Chemicals: High enthalpy heat  processing'

# read the corresponding lines
s_fec = excel_fec.iloc[116:129,year]

s_ued = excel_ued.iloc[116:129,year]

assert s_fec.index[0] == subsector

eff_elec = s_ued['High enthalpy heat  processing - Electric (microwave)']/s_fec['High enthalpy heat  processing - Electric (microwave)']

df.loc['elec',sector] += s_ued[subsector]/eff_elec

#### Chemicals: Furnaces
#> assume fully electrified

subsector = 'Chemicals: Furnaces'

# read the corresponding lines
s_fec = excel_fec.iloc[129:137,year]

s_ued = excel_ued.iloc[129:137,year]

assert s_fec.index[0] == subsector

eff = s_ued['Chemicals: Furnaces - Electric']/s_fec['Chemicals: Furnaces - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff

#### Chemicals: Process cooling
#> assume fully electrified

subsector = 'Chemicals: Process cooling'

# read the corresponding lines
s_fec = excel_fec.iloc[137:151,year]

s_ued = excel_ued.iloc[137:151,year]

assert s_fec.index[0] == subsector

eff_elec = s_ued['Chemicals: Process cooling - Electric']/s_fec['Chemicals: Process cooling - Electric']

df.loc['elec',sector] += s_ued[subsector]/eff_elec

#### Chemicals: Generic electric process

subsector = 'Chemicals: Generic electric process'

# read the corresponding lines
s_fec = excel_fec.iloc[151:152,year]

assert s_fec.index[0] == subsector

df.loc['elec',sector] += s_fec[subsector]

# read the corresponding lines
s_out = excel_out.iloc[10:11,year]

# check the position
assert sector in str(s_out.index)

df.loc['process emission',sector] += 0 # unit tCO2/t material

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat', 'naphtha']

df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out.values # unit MWh/t material
# 1 ktoe = 11630 MWh

## Non-metallic mineral products
#
#> This includes cement, ceramic and glass production.
#
#> This sector includes process-emissions related to the fabrication of clinker.

sector = 'Non-metallic mineral products'

# read the input sheets
excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True)

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True)

excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True)

excel_emi = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_emi',
                      index_col=0,header=0,squeeze=True)

### Cement
#
#> This sector has process-emissions.
#
#> Includes three subcategories: (a) Grinding, milling of raw material, (b) Pre-heating and pre-calcination, (c) clinker production (kilns), (d) Grinding, packaging. (b)+(c) represent 94% of fec. So (a) is joined to (b) and (d) is joined to (c).
#
#> Temperatures above 1400C are required for procesing limestone and sand into clinker.
#
#> Everything (except current electricity and heat consumption and existing biomass) is transformed into methane for high T.

sector = 'Cement'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:25,year]

s_ued = excel_ued.iloc[3:25,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec',sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat',sector] += s_fec['Low enthalpy heat']

# pre-processing: keep existing elec and biomass, rest to methane
df.loc['elec', sector] += s_fec['Cement: Grinding, milling of raw material']
df.loc['biomass', sector] += s_fec['Biomass']
df.loc['methane', sector] += s_fec['Cement: Pre-heating and pre-calcination'] - s_fec['Biomass']


#### Cement: Clinker production (kilns)

subsector = 'Cement: Clinker production (kilns)'

# read the corresponding lines
s_fec = excel_fec.iloc[34:43,year]

s_ued = excel_ued.iloc[34:43,year]

assert s_fec.index[0] == subsector

df.loc['biomass', sector] += s_fec['Biomass']
df.loc['methane', sector] += s_fec['Cement: Clinker production (kilns)'] - s_fec['Biomass']
df.loc['elec', sector] += s_fec['Cement: Grinding, packaging']


#### Process-emission came from the calcination of limestone to chemically reactive calcium oxide (lime).
#> Calcium carbonate -> lime + CO2
#
#> CaCO3  -> CaO + CO2

s_emi = excel_emi.iloc[3:44,year]

assert s_emi.index[0] == sector

s_out = excel_out.iloc[7:8,year]

assert sector in str(s_out.index)

df.loc['process emission',sector] +=s_emi['Process emissions']/s_out.values # unit tCO2/t material

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']

df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out.values # unit MWh/t material

### Ceramics & other NMM
#
#> This sector has process emissions.
#
#> Includes four subcategories: (a) Mixing of raw material, (b) Drying and sintering of raw material, (c) Primary production process, (d) Product finishing. (b)represents 65% of fec and (a) 4%. So (a) is joined to (b).
#
#> Everything is electrified

sector = 'Ceramics & other NMM'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[45:94,year]

s_ued = excel_ued.iloc[45:94,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Ceramics: Microwave drying and sintering']/s_fec['Ceramics: Microwave drying and sintering']
df.loc['elec', sector] += s_ued[['Ceramics: Mixing of raw material','Ceramics: Drying and sintering of raw material']].sum()/eff_elec

eff_elec=s_ued['Ceramics: Electric kiln']/s_fec['Ceramics: Electric kiln']
df.loc['elec', sector] += s_ued['Ceramics: Primary production process']/eff_elec

eff_elec=s_ued['Ceramics: Electric furnace']/s_fec['Ceramics: Electric furnace']
df.loc['elec', sector] += s_ued['Ceramics: Product finishing']/eff_elec

s_emi = excel_emi.iloc[45:94,year]

assert s_emi.index[0] == sector

s_out = excel_out.iloc[8:9,year]

assert sector in str(s_out.index)

df.loc['process emission',sector] += s_emi['Process emissions']/s_out.values # unit tCO2/t material

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']

df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out.values # unit MWh/t material
# 1 ktoe = 11630 MWh

### Glass production
#
#> This sector has process emissions.
#
#> Includes four subcategories: (a) Melting tank, (b) Forming, (c) Annealing, (d) Finishing processes. (a)represents 73%. (b), (d) are joined to (c).
#
#> Everything is electrified.

sector = 'Glass production'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[95:123,year]

s_ued = excel_ued.iloc[95:123,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Glass: Electric melting tank']/s_fec['Glass: Electric melting tank']
df.loc['elec', sector] += s_ued['Glass: Melting tank']/eff_elec

eff_elec=s_ued['Glass: Annealing - electric']/s_fec['Glass: Annealing - electric']
df.loc['elec', sector] += s_ued[['Glass: Forming','Glass: Annealing','Glass: Finishing processes']].sum()/eff_elec

s_emi = excel_emi.iloc[95:124,year]

assert s_emi.index[0] == sector

s_out = excel_out.iloc[9:10,year]

assert sector in str(s_out.index)

df.loc['process emission',sector] += s_emi['Process emissions']/s_out.values # unit tCO2/t material

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']

df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out.values # unit MWh/t material

## Pulp, paper and printing
#
#> Pulp, paper and printing can be completely electrified.
#
#> There are no process emissions associated to this sector.

sector = 'Pulp, paper and printing'

# read the input sheets
excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True)

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True)

excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True)

### Pulp production
#
#> Includes three subcategories: (a) Wood preparation, grinding; (b) Pulping;  (c) Cleaning.
#
#> (b) Pulping is either biomass or electric; left like this (dominated by biomass).
#
#> (a) Wood preparation, grinding and (c) Cleaning represent only 10% their current energy consumption is assumed to be electrified without any change in efficiency

sector = 'Pulp production'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:28,year]

s_ued = excel_ued.iloc[3:28,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Industry-specific
df.loc['elec', sector] += s_fec[['Pulp: Wood preparation, grinding', 'Pulp: Cleaning', 'Pulp: Pulping electric']].sum()

# Efficiency changes due to biomass
eff_bio=s_ued['Biomass']/s_fec['Biomass']
df.loc['biomass', sector] += s_ued['Pulp: Pulping thermal']/eff_bio

s_out = excel_out.iloc[8:9,year]

assert sector in str(s_out.index)

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']

df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Pulp production (kt)'] # unit MWh/t material

### Paper production
#
#> Includes three subcategories: (a) Stock preparation; (b) Paper machine;  (c) Product finishing.
#
#> (b) Paper machine and (c) Product finishing are left electric and thermal is moved to biomass. The efficiency is calculated from the pulping process that is already biomass.
#
#> (a) Stock preparation represents only 7% and its current energy consumption is assumed to be electrified without any change in efficiency.

sector = 'Paper production'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[29:78,year]

s_ued = excel_ued.iloc[29:78,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Industry-specific
df.loc['elec', sector] += s_fec['Paper: Stock preparation']

# add electricity from process that is already electrified
df.loc['elec', sector] += s_fec['Paper: Paper machine - Electricity']

# add electricity from process that is already electrified
df.loc['elec', sector] += s_fec['Paper: Product finishing - Electricity']


s_fec = excel_fec.iloc[53:64,year]

s_ued = excel_ued.iloc[53:64,year]

assert s_fec.index[0] == 'Paper: Paper machine - Steam use'

# Efficiency changes due to biomass
eff_bio=s_ued['Biomass']/s_fec['Biomass']
df.loc['biomass', sector] += s_ued['Paper: Paper machine - Steam use']/eff_bio


s_fec = excel_fec.iloc[66:77,year]

s_ued = excel_ued.iloc[66:77,year]

assert s_fec.index[0] == 'Paper: Product finishing - Steam use'

# Efficiency changes due to biomass
eff_bio=s_ued['Biomass']/s_fec['Biomass']
df.loc['biomass', sector] += s_ued['Paper: Product finishing - Steam use']/eff_bio


# read the corresponding lines
s_out = excel_out.iloc[9:10,year]

assert sector in str(s_out.index)

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']

df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out.values # unit MWh/t material\

### Printing and media reproduction
#
#> (a) Printing and publishing is assumed to be electrified without any change in efficiency.

sector='Printing and media reproduction'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[79:90,year]

s_ued = excel_ued.iloc[79:90,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec',sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()
df.loc['elec',sector] += s_ued[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat',sector] += s_fec['Low enthalpy heat']
df.loc['heat',sector] += s_ued['Low enthalpy heat']

# Industry-specific
df.loc['elec', sector] += s_fec['Printing and publishing']
df.loc['elec', sector] += s_ued['Printing and publishing']

# read the corresponding lines
s_out = excel_out.iloc[10:11,year]

assert sector in str(s_out.index)

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']

df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out.values # unit MWh/t material

## Food, beverages and tobaco
#
#> Food, beverages and tobaco can be completely electrified.
#
#> There are no process emissions associated to this sector.

sector = 'Food, beverages and tobacco'

# read the input sheets
excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True)

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True)

excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True)

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:78,year]

s_ued = excel_ued.iloc[3:78,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Food: Direct Heat - Electric']/s_fec['Food: Direct Heat - Electric']
df.loc['elec', sector] += s_ued['Food: Oven (direct heat)']/eff_elec

eff_elec=s_ued['Food: Process Heat - Electric']/s_fec['Food: Process Heat - Electric']
df.loc['elec', sector] += s_ued['Food: Specific process heat']/eff_elec

eff_elec=s_ued['Food: Electric drying']/s_fec['Food: Electric drying']
df.loc['elec', sector] += s_ued['Food: Drying']/eff_elec

eff_elec=s_ued['Food: Electric cooling']/s_fec['Food: Electric cooling']
df.loc['elec', sector] += s_ued['Food: Process cooling and refrigeration']/eff_elec

# Steam processing goes all to biomass without change in efficiency
df.loc['biomass', sector] += s_fec['Food: Steam processing']

# add electricity from process that is already electrified
df.loc['elec', sector] += s_fec['Food: Electric machinery']

# read the corresponding lines
s_out = excel_out.iloc[3:4,year]

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']

df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Physical output (index)'] # unit MWh/t material

## Non Ferrous Metals

sector = 'Non Ferrous Metals'

# read the input sheets
excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True)

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True)

excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True)

excel_emi = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_emi',
                      index_col=0,header=0,squeeze=True) # the emission sheet

### Alumina
#
#> High enthalpy heat is converted to methane. Process heat at T>500ºC is required here.
#
#> Refining is electrified.
#
#> There are no process emissions associated to Alumina manufacturing

sector = 'Alumina production'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:31,year]

s_ued = excel_ued.iloc[3:31,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# High-enthalpy heat is transformed into methane
s_fec = excel_fec.iloc[13:24,year]

s_ued = excel_ued.iloc[13:24,year]

assert s_fec.index[0] == 'Alumina production: High enthalpy heat'

eff_met=s_ued['Natural gas (incl. biogas)']/s_fec['Natural gas (incl. biogas)']
df.loc['methane', sector] += s_fec['Alumina production: High enthalpy heat']/eff_met

# Efficiency changes due to electrification
s_fec = excel_fec.iloc[24:30,year]

s_ued = excel_ued.iloc[24:30,year]

assert s_fec.index[0] == 'Alumina production: Refining'

eff_elec=s_ued['Electricity']/s_fec['Electricity']
df.loc['elec', sector] += s_ued['Alumina production: Refining']/eff_elec

# read the corresponding lines
s_out = excel_out.iloc[9:10,year]

assert sector in str(s_out.index)

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']
df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Alumina production (kt)'] # unit MWh/t material

### Aluminium primary route
#
#> Production through the primary route is divided into 50% remains as today and 50% is transformed into secondary route

sector = 'Aluminium - primary production'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[31:66,year]

s_ued = excel_ued.iloc[31:66,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Add aluminium  electrolysis (smelting
df.loc['elec', sector] += s_fec['Aluminium electrolysis (smelting)']

# Efficiency changes due to electrification
eff_elec=s_ued['Aluminium processing - Electric']/s_fec['Aluminium processing - Electric']
df.loc['elec', sector] += s_ued['Aluminium processing  (metallurgy e.g. cast house, reheating)']/eff_elec

# Efficiency changes due to electrification
eff_elec=s_ued['Aluminium finishing - Electric']/s_fec['Aluminium finishing - Electric']
df.loc['elec', sector] += s_ued['Aluminium finishing']/eff_elec

s_emi = excel_emi.iloc[31:67,year]

assert s_emi.index[0] == sector

s_out = excel_out.iloc[11:12,year]

assert sector in str(s_out.index)

df.loc['process emission',sector] = s_emi['Process emissions']/s_out['Aluminium - primary production'] # unit tCO2/t material

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']
df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Aluminium - primary production'] # unit MWh/t material

### Aluminium secondary route
#
#> All is coverted into secondary route fully electrified

sector = 'Aluminium - secondary production'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[68:109,year]

s_ued = excel_ued.iloc[68:109,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Secondary aluminium - Electric']/s_fec['Secondary aluminium - Electric']
df.loc['elec', sector] += s_ued['Secondary aluminium (incl. pre-treatment, remelting)']/eff_elec

# Efficiency changes due to electrification
eff_elec=s_ued['Aluminium processing - Electric']/s_fec['Aluminium processing - Electric']
df.loc['elec', sector] += s_ued['Aluminium processing  (metallurgy e.g. cast house, reheating)']/eff_elec

# Efficiency changes due to electrification
eff_elec=s_ued['Aluminium finishing - Electric']/s_fec['Aluminium finishing - Electric']
df.loc['elec', sector] += s_ued['Aluminium finishing']/eff_elec

# read the corresponding lines
s_out = excel_out.iloc[12:13,year]

assert sector in str(s_out.index)

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']
df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Aluminium - secondary production'] # unit MWh/t material
# 1 ktoe = 11630 MWh


### Other non-ferrous metals

sector = 'Other non-ferrous metals'

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[110:152,year]

s_ued = excel_ued.iloc[110:152,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Metal production - Electric']/s_fec['Metal production - Electric']
df.loc['elec', sector] += s_ued['Other Metals: production']/eff_elec

# Efficiency changes due to electrification
eff_elec=s_ued['Metal processing - Electric']/s_fec['Metal processing - Electric']
df.loc['elec', sector] += s_ued['Metal processing  (metallurgy e.g. cast house, reheating)']/eff_elec

# Efficiency changes due to electrification
eff_elec=s_ued['Metal finishing - Electric']/s_fec['Metal finishing - Electric']
df.loc['elec', sector] += s_ued['Metal finishing']/eff_elec

s_emi = excel_emi.iloc[110:153,year]

assert s_emi.index[0] == sector

s_out = excel_out.iloc[13:14,year]

assert sector in str(s_out.index)

df.loc['process emission',sector] = s_emi['Process emissions']/s_out['Other non-ferrous metals (kt lead eq.)'] # unit tCO2/t material

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']
df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Other non-ferrous metals (kt lead eq.)'] # unit MWh/t material

## Transport Equipment

sector = 'Transport Equipment'
# read the input sheets
excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True)

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True)

excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True)

excel_emi = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_emi',
                      index_col=0,header=0,squeeze=True) # the emission sheet

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:45,year]

s_ued = excel_ued.iloc[3:45,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Trans. Eq.: Electric Foundries']/s_fec['Trans. Eq.: Electric Foundries']
df.loc['elec', sector] += s_ued['Trans. Eq.: Foundries']/eff_elec

# Efficiency changes due to electrification
eff_elec=s_ued['Trans. Eq.: Electric connection']/s_fec['Trans. Eq.: Electric connection']
df.loc['elec', sector] += s_ued['Trans. Eq.: Connection techniques']/eff_elec

# Efficiency changes due to electrification
eff_elec=s_ued['Trans. Eq.: Heat treatment - Electric']/s_fec['Trans. Eq.: Heat treatment - Electric']
df.loc['elec', sector] += s_ued['Trans. Eq.: Heat treatment']/eff_elec

df.loc['elec', sector] += s_fec['Trans. Eq.: General machinery']
df.loc['elec', sector] += s_fec['Trans. Eq.: Product finishing']

# Steam processing is supplied with biomass
eff_biomass=s_ued['Biomass']/s_fec['Biomass']
df.loc['biomass', sector] += s_ued['Trans. Eq.: Steam processing']/eff_biomass

# read the corresponding lines
s_out = excel_out.iloc[3:4,year]
# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']
df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Physical output (index)'] # unit MWh/t material
# 1 ktoe = 11630 MWh

## Machinery Equipment

sector = 'Machinery Equipment'

# read the input sheets
excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True)

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True)

excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True)

excel_emi = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_emi',
                      index_col=0,header=0,squeeze=True) # the emission sheet

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:45,year]

s_ued = excel_ued.iloc[3:45,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Mach. Eq.: Electric Foundries']/s_fec['Mach. Eq.: Electric Foundries']
df.loc['elec', sector] += s_ued['Mach. Eq.: Foundries']/eff_elec

# Efficiency changes due to electrification
eff_elec=s_ued['Mach. Eq.: Electric connection']/s_fec['Mach. Eq.: Electric connection']
df.loc['elec', sector] += s_ued['Mach. Eq.: Connection techniques']/eff_elec

# Efficiency changes due to electrification
eff_elec=s_ued['Mach. Eq.: Heat treatment - Electric']/s_fec['Mach. Eq.: Heat treatment - Electric']
df.loc['elec', sector] += s_ued['Mach. Eq.: Heat treatment']/eff_elec

df.loc['elec', sector] += s_fec['Mach. Eq.: General machinery']
df.loc['elec', sector] += s_fec['Mach. Eq.: Product finishing']

# Steam processing is supplied with biomass
eff_biomass=s_ued['Biomass']/s_fec['Biomass']
df.loc['biomass', sector] += s_ued['Mach. Eq.: Steam processing']/eff_biomass

# read the corresponding lines
s_out = excel_out.iloc[3:4,year]

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']
df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Physical output (index)'] # unit MWh/t material

## Textiles and leather

sector = 'Textiles and leather'
# read the input sheets
excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True)

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True)

excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True)

excel_emi = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_emi',
                      index_col=0,header=0,squeeze=True) # the emission sheet

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:57,year]

s_ued = excel_ued.iloc[3:57,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Textiles: Electric drying']/s_fec['Textiles: Electric drying']
df.loc['elec', sector] += s_ued['Textiles: Drying']/eff_elec

df.loc['elec', sector] += s_fec['Textiles: Electric general machinery']
df.loc['elec', sector] += s_fec['Textiles: Finishing Electric']

# Steam processing is supplied with biomass
eff_biomass=s_ued[15:26]['Biomass']/s_fec[15:26]['Biomass']
df.loc['biomass', sector] += s_ued['Textiles: Pretreatment with steam']/eff_biomass
df.loc['biomass', sector] += s_ued['Textiles: Wet processing with steam']/eff_biomass

# read the corresponding lines
s_out = excel_out.iloc[3:4,year]

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']
df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Physical output (index)'] # unit MWh/t material

## Wood and wood products

sector = 'Wood and wood products'
# read the input sheets
excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True)

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True)

excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True)

excel_emi = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_emi',
                      index_col=0,header=0,squeeze=True) # the emission sheet

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:46,year]

s_ued = excel_ued.iloc[3:46,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Wood: Electric drying']/s_fec['Wood: Electric drying']
df.loc['elec', sector] += s_ued['Wood: Drying']/eff_elec

df.loc['elec', sector] += s_fec['Wood: Electric mechanical processes']
df.loc['elec', sector] += s_fec['Wood: Finishing Electric']

# Steam processing is supplied with biomass
eff_biomass=s_ued[15:25]['Biomass']/s_fec[15:25]['Biomass']
df.loc['biomass', sector] += s_ued['Wood: Specific processes with steam']/eff_biomass

# read the corresponding lines
s_out = excel_out.iloc[3:4,year]

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']
df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Physical output (index)'] # unit MWh/t material

##  Other Industrial Sectors

sector = 'Other Industrial Sectors'
# read the input sheets
excel_fec = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_fec',
                      index_col=0,header=0,squeeze=True)

excel_ued = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_ued',
                      index_col=0,header=0,squeeze=True)

excel_out = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector],
                      index_col=0,header=0,squeeze=True)

excel_emi = pd.read_excel('{}/JRC-IDEES-2015_Industry_{}.xlsx'.format(base_dir,country), sheet_name=sub_sheet_name_dict[sector]+'_emi',
                      index_col=0,header=0,squeeze=True) # the emission sheet

df[sector] = 0

# read the corresponding lines
s_fec = excel_fec.iloc[3:67,year]

s_ued = excel_ued.iloc[3:67,year]

assert s_fec.index[0] == sector

# Lighting, Air compressors, Motor drives, Fans and pumps
df.loc['elec', sector] += s_fec[['Lighting','Air compressors','Motor drives','Fans and pumps']].sum()

# Low enthalpy heat
df.loc['heat', sector] += s_fec['Low enthalpy heat']

# Efficiency changes due to electrification
eff_elec=s_ued['Other Industrial sectors: Electric processing']/s_fec['Other Industrial sectors: Electric processing']
df.loc['elec', sector] += s_ued['Other Industrial sectors: Process heating']/eff_elec

eff_elec=s_ued['Other Industries: Electric drying']/s_fec['Other Industries: Electric drying']
df.loc['elec', sector] += s_ued['Other Industrial sectors: Drying']/eff_elec

eff_elec=s_ued['Other Industries: Electric cooling']/s_fec['Other Industries: Electric cooling']
df.loc['elec', sector] += s_ued['Other Industrial sectors: Process Cooling']/eff_elec

# Diesel motors are electrified
df.loc['elec', sector] += s_fec['Other Industrial sectors: Diesel motors (incl. biofuels)']
df.loc['elec', sector] += s_fec['Other Industrial sectors: Electric machinery']

# Steam processing is supplied with biomass
eff_biomass=s_ued[15:25]['Biomass']/s_fec[15:25]['Biomass']
df.loc['biomass', sector] += s_ued['Other Industrial sectors: Steam processing']/eff_biomass

# read the corresponding lines
s_out = excel_out.iloc[3:4,year]

# final energy consumption per t
sources=['elec','biomass', 'methane', 'hydrogen', 'heat','naphtha']
df.loc[sources,sector] = df.loc[sources,sector]*conv_factor/s_out['Physical output (index)'] # unit MWh/t material


df.index.name = "MWh/tMaterial"
df.to_csv('resources/industry_sector_ratios.csv')
