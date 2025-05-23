# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:58:15 2024

@author: alice
"""

import pypsa 
import pandas as pd

def country_values(df):
    df.index = df.index.str[:2]
    grouped = df.groupby(df.index).sum()
    grouped = grouped[abs(grouped) > 1e-4]
    return grouped

n = pypsa.Network("../results/baseline_eu_dem/postnetworks/base_s_39_lvopt___2030.nc")

timestep = n.snapshot_weightings.iloc[0,0]
alinks = n.links

# %%
# Check if cement demand is met by cement supply
cement_load = n.loads_t.p.filter(like="cement", axis=1).sum()*timestep #kton cement per year per node
cement_load = country_values(cement_load)
cement_supply = -n.links_t.p1.filter(like="Cement",axis=1).sum()*timestep #kton cement per year per node 
cement_supply = country_values(cement_supply)

# Check if energy consumption and emission balances are respected
cement_links_p0 = n.links_t.p0.filter(like="Cement",axis=1).sum()*timestep #kton limestone per year per node
cement_links_p0 = country_values(cement_links_p0)
cement_links_p2 = n.links_t.p2.filter(like="Cement",axis=1).sum()*timestep #MWh industry heat per year per node
cement_links_p2 = country_values(cement_links_p2)
cement_links_p3 = n.links_t.p3.filter(like="Cement",axis=1).sum()*timestep #tCO2 per yaer per node
cement_links_p3 = country_values(cement_links_p3)

efficiency1 = n.links.efficiency.loc[n.links.index.str.contains('Cement', case=False)].iloc[0]
efficiency2 = n.links.efficiency2.loc[n.links.index.str.contains('Cement', case=False)].iloc[0]
efficiency3 = n.links.efficiency3.loc[n.links.index.str.contains('Cement', case=False)].iloc[0]

check_load = (cement_load - cement_supply).sum()
check0_1 = (cement_links_p0 * efficiency1 - cement_supply).sum()
check0_2 = (cement_links_p0 * (- efficiency2) - cement_links_p2).sum()
check0_3 = (cement_links_p0 * (-efficiency3) - cement_links_p3).sum()

# Check conditions and print message if any are above 0.01
if any(check > 0.01 for check in [check_load, check0_1, check0_2, check0_3]):
    print("Check cement")
    
# %%
    
# Check if emissions balances are okay
cement_process_cc_p0 = n.links_t.p0.filter(like="cement CC",axis=1).sum()*timestep #tCO2 entering CC

cement_process_cc_p1 = n.links_t.p1.filter(like="cement CC",axis=1).sum()*timestep #tCO2 to atmosphere
cement_process_cc_p1 = country_values(cement_process_cc_p1)

cement_process_cc_p2 = n.links_t.p2.filter(like="cement CC",axis=1).sum()*timestep #tCO2 to store
cement_process_cc_p2 = country_values(cement_process_cc_p2)

cement_process_cc_p3 = n.links_t.p3.filter(like="cement CC",axis=1).sum()*timestep #MWh heat required for carbon capture
cement_process_cc_p3 = country_values(cement_process_cc_p3)

efficiency1 = n.links.efficiency.loc[n.links.index.str.contains('cement CC', case=False)].iloc[0]
efficiency2 = n.links.efficiency2.loc[n.links.index.str.contains('cement CC', case=False)].iloc[0]
efficiency3 = n.links.efficiency3.loc[n.links.index.str.contains('cement CC', case=False)].iloc[0]

check0_1 = (cement_process_cc_p0 * efficiency1 - (-cement_process_cc_p1)).sum()
check0_2 = (cement_process_cc_p0 * (- efficiency2) - cement_process_cc_p2).sum()
check0_3 = (cement_process_cc_p0 * (-efficiency3) - cement_process_cc_p3).sum()

cement_process_p0 = n.links_t.p0.filter(like="cement process emis to atmosphere",axis=1).sum()*timestep #kton limestone per year per node
cement_process_p0 = country_values(cement_process_p0)

cement_process_p1 = n.links_t.p1.filter(like="cement process emis to atmosphere",axis=1).sum()*timestep #kton limestone per year per node
cement_process_p1 = country_values(cement_process_p1)


tot_cement_process = -(cement_process_p1 + cement_process_cc_p1).fillna(cement_process_p1).sum()/1e6 #Mton CO2 output of cement industry
# https://www.umweltbundesamt.de/sites/default/files/medien/1410/publikationen/2022-01-04_climate-change_02-2022_decomposition_of_co2_emissions_in_the_european_cement_sector_0.pdf page 5 120Mt in 2018

share_noncc = cement_process_p1.sum() / (cement_process_p1.sum() + cement_process_cc_p1.sum())
dac = n.links_t.p2.filter(like="DAC",axis = 1).sum()*timestep
dac = country_values(dac)
dac = dac[abs(dac)> 0.5]

# Check conditions and print message if any are above 0.01
if any(check > 0.01 for check in [check_load, check0_1, check0_2, check0_3]):
    print("Check cement emissions")


# %%

# Check if steel demand is met by steel supply
steel_load = n.loads_t.p.filter(like="steel", axis=1).sum()*timestep #kton steel per year per node
steel_load = country_values(steel_load)

steel_supply_bof = -n.links_t.p1.filter(like="BOF", axis=1).sum() * timestep  # kton steel per year per node
steel_supply_bof = country_values(steel_supply_bof)

steel_supply_eaf = -n.links_t.p1.filter(like="DRI-EAF", axis=1).sum() * timestep  # kton steel per year per node
steel_supply_eaf = country_values(steel_supply_eaf)


# Check if energy consumption and emission balances are respected for BOF
bof_links_p0 = n.links_t.p0.filter(like="BOF",axis=1).sum()*timestep #kton iron per year per node
bof_links_p0 = country_values(bof_links_p0)
bof_links_p2 = n.links_t.p2.filter(like="BOF",axis=1).sum()*timestep #kton coal per year per node
bof_links_p2 = country_values(bof_links_p2)
bof_links_p3 = n.links_t.p3.filter(like="BOF",axis=1).sum()*timestep #MWh industry heat per year per node
bof_links_p3 = country_values(bof_links_p3)
bof_links_p4 = n.links_t.p4.filter(like="BOF",axis=1).sum()*timestep #MWh electricity per year per node
bof_links_p4 = country_values(bof_links_p4)
bof_links_p5 = n.links_t.p5.filter(like="BOF",axis=1).sum()*timestep #tCO2 per year per node
bof_links_p5 = country_values(bof_links_p5)

efficiency1 = n.links.efficiency.loc[n.links.index.str.contains('BOF', case=False)].iloc[0]
efficiency2 = n.links.efficiency2.loc[n.links.index.str.contains('BOF', case=False)].iloc[0]
efficiency3 = n.links.efficiency3.loc[n.links.index.str.contains('BOF', case=False)].iloc[0]
efficiency4 = n.links.efficiency4.loc[n.links.index.str.contains('BOF', case=False)].iloc[0]
efficiency5 = n.links.efficiency5.loc[n.links.index.str.contains('BOF', case=False)].iloc[0]

check0_1 = (bof_links_p0 * efficiency1 - steel_supply_bof).sum() # Iron to steel
check0_2 = (bof_links_p0 * (-efficiency2) - bof_links_p2).sum()  # Coal to steel
check0_3 = (bof_links_p0 * (-efficiency3) - bof_links_p3).sum()  # Heat to steel
check0_4 = (bof_links_p0 * (-efficiency4) - bof_links_p4).sum()  # Electricity to steel
check0_5 = (bof_links_p0 * (-efficiency5) - bof_links_p5).sum()  # Emissions to steel

# Check conditions and print message if any are above 0.01
if any(check > 0.01 for check in [check0_1, check0_2, check0_3, check0_4, check0_5]):
    print("Check BOF")

# Check if energy consumption and emission balances are respected for NG EAF
eaf_links_p0 = n.links_t.p0.filter(like="EAF",axis=1).sum()*timestep #kton iron per year per node
eaf_links_p0 = country_values(eaf_links_p0)
eaf_links_p2 = n.links_t.p2.filter(like="EAF",axis=1).sum()*timestep #kton coal per year per node
eaf_links_p2 = country_values(eaf_links_p2)
eaf_links_p3 = n.links_t.p3.filter(like="EAF",axis=1).sum()*timestep #MWh industry heat per year per node
eaf_links_p3 = country_values(eaf_links_p3)
eaf_links_p4 = n.links_t.p4.filter(like="EAF",axis=1).sum()*timestep #MWh electricity per year per node
eaf_links_p4 = country_values(eaf_links_p4)

efficiency1 = n.links.efficiency.loc[n.links.index.str.contains('EAF', case=False)].iloc[0]
efficiency2 = n.links.efficiency2.loc[n.links.index.str.contains('EAF', case=False)].iloc[0]
efficiency3 = n.links.efficiency3.loc[n.links.index.str.contains('EAF', case=False)].iloc[0]
efficiency4 = n.links.efficiency4.loc[n.links.index.str.contains('EAF', case=False)].iloc[0]

check_load = steel_supply_eaf.sum() + steel_supply_bof.sum() - steel_load.sum()
check0_1 = (eaf_links_p0 * efficiency1 - steel_supply_eaf).sum() # Iron to steel
check0_2 = (eaf_links_p0 * (-efficiency2) - eaf_links_p2).sum()  # Coal to steel
check0_3 = (eaf_links_p0 * (-efficiency3) - eaf_links_p3).sum()  # Heat to steel
check0_4 = (eaf_links_p0 * (-efficiency4) - eaf_links_p4).sum()  # Electricity to steel


# Check conditions and print message if any are above 0.01
if any(check > 0.01 for check in [check0_1, check0_2, check0_3, check0_4, check_load]):
    print("Check EAF or load")
    
# Chech DRI gas
ng_drigas = n.links_t.p0.filter(like="CH4 to DRI", axis = 1).sum()*timestep
h2_drigas = n.links_t.p0.filter(like="H2 to DRI", axis = 1).sum()*timestep
ng_efficiency = n.links.efficiency.loc[n.links.index.str.contains('CH4 to DRI', case=False)].iloc[0]
h2_efficiency = n.links.efficiency.loc[n.links.index.str.contains('H2 to DRI', case=False)].iloc[0]

# They should match with bus2 in EAF
check_dri = ng_drigas.sum()*ng_efficiency + h2_drigas.sum()*h2_efficiency - eaf_links_p2.sum()

# Check if emissions balances are okay
dri_process_p0 = n.links_t.p0.filter(like="steel dri process",axis=1).sum()*timestep #kton limestone per year per node
dri_process_p0 = country_values(dri_process_p0)
dri_process_p1 = n.links_t.p1.filter(like="steel dri process",axis=1).sum()*timestep #kton limestone per year per node
dri_process_p1 = country_values(dri_process_p1)
dri_efficiency = n.links.efficiency.loc[n.links.index.str.contains('steel dri process', case=False)].iloc[0]

check_dri_emis = dri_process_p0.sum() * dri_efficiency - dri_process_p1.sum()

# Check BOF emissions
bof_process_p0 = n.links_t.p0.filter(like="steel bof process",axis=1).sum()*timestep #kton limestone per year per node
bof_process_p0 = country_values(bof_process_p0)
bof_process_p1 = n.links_t.p1.filter(like="steel bof process",axis=1).sum()*timestep #kton limestone per year per node
bof_process_p1 = country_values(bof_process_p1)
bof_efficiency = n.links.efficiency.loc[n.links.index.str.contains('steel bof process', case=False)].iloc[0]

check_bof_emis = bof_process_p0.sum() * bof_efficiency - bof_process_p1.sum()

bof_cc_p0 = n.links_t.p0.filter(like="steel bof CC",axis=1).sum()*timestep #kton limestone per year per node
bof_cc_p0 = country_values(bof_cc_p0)
bof_cc_p1 = n.links_t.p1.filter(like="steel bof CC",axis=1).sum()*timestep #kton limestone per year per node
bof_cc_p1 = country_values(bof_cc_p1)
bof_cc_p2 = n.links_t.p2.filter(like="steel bof CC",axis=1).sum()*timestep #kton limestone per year per node
bof_cc_p2 = country_values(bof_cc_p2)
bof_cc_efficiency1 = n.links.efficiency.loc[n.links.index.str.contains('steel bof CC', case=False)].iloc[0]
bof_cc_efficiency2 = n.links.efficiency2.loc[n.links.index.str.contains('steel bof CC', case=False)].iloc[0]

check_bof_cc0_1 = bof_cc_p0.sum() * bof_cc_efficiency1 - bof_cc_p1.sum()
check_bof_cc0_2 = bof_cc_p0.sum() * bof_cc_efficiency2 - bof_cc_p2.sum()



# %% Ammonia

ammonia_prod = n.links[n.links.index.str.contains("Haber") ]
ammonia_load = n.loads_t.p.filter(like="NH3", axis=1)
