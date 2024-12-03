# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:10:30 2024

@author: alice
"""

import pypsa
import pandas as pd
import os

def country_values(df):
    df.index = df.index.str[:2]
    grouped = df.groupby(df.index).sum()
    grouped = grouped[abs(grouped) > 1e-4]
    return grouped

def retrieve_loads(df, suffix):

    filtered_df = df[[col for col in df.columns if col.endswith(suffix)]]
    filtered_df.columns = filtered_df.columns.str.split(" 0").str[0] + " 0"
    return filtered_df


def retrieve_links(df, suffix):

    filtered_df = df.loc[
        :, n.links[(n.links.index.str.contains(suffix, case=False))].index
    ]  # .sum()
    filtered_df.columns = filtered_df.columns.str.split(" 0").str[0] + " 0"
    filtered_df = filtered_df.T.groupby(filtered_df.columns).sum().T
    filtered_df = filtered_df[
        [col for col in filtered_df.columns if not col.startswith("EU")]
    ]
    return filtered_df


n = pypsa.Network("../results_24h/baseline_eu_dem/postnetworks/base_s_39_lvopt___2030.nc")

timestep = n.snapshot_weightings.iloc[0,0]

# %% STEEL
countries = n.buses[n.buses.index.str.endswith(' 0')]
countries.index = countries.index.str[:2]
countries = countries.index
countries = countries[~countries.duplicated()]
years = [2030,2040,2050]
df = pd.DataFrame(0, index=countries,columns=years)

# %%

# Initialize Excel writer
scenarios = ["baseline_eu_dem", "policy_eu_dem", "baseline_regional_dem", "policy_regional_dem"]
os.makedirs("excels_24h", exist_ok=True)

for scenario in scenarios:
    
    output_filename = scenario + ".xlsx"
    with pd.ExcelWriter("excels_24h/" + output_filename, engine="xlsxwriter") as writer:
        
        # Create empty DataFrames for each type of data
        eaf_capacity_df = df.copy()
        bof_capacity_df = df.copy()
        h2_eaf_df = df.copy()
        ng_eaf_df = df.copy()
        tgr_capacity_df = df.copy()
        eaf_prod_df = df.copy()
        bof_prod_df = df.copy()
        elec_prices_df = df.copy()
        ammonia_prod_df = df.copy()
        steel_proc_emis_df = df.copy()
        bofcc_emis_capt_df = df.copy()
        cement_proc_emis_df = df.copy()
        
        # Hydrogen
        elec_h2prod_df = df.copy()
        smrcc_h2prod_df = df.copy()
        smr_h2prod_df = df.copy()
        nh3crack_h2prod_df = df.copy()
        
        if 'eu' in scenario:
            steel_load_df = pd.DataFrame(0, index=['EU'], columns = years)
            ammonia_load_df = pd.DataFrame(0, index=['EU'], columns = years)
        else:
            steel_load_df = df.copy()
            ammonia_load_df = df.copy()
            
        for year in years:
            
            n = pypsa.Network(f"../results_24h/{scenario}/postnetworks/base_s_39_lvopt___{year}.nc")
            timestep = n.snapshot_weightings.iloc[0,0]
            
            # EAF capacity
            eaf_efficiency = n.links.efficiency.loc[n.links.index.str.contains('EAF', case=False)].iloc[0]
            eaf_capacity = n.links.p_nom_opt[n.links.index.str.contains("EAF")] * 8760 * eaf_efficiency  # kton steel / year
            eaf_capacity = country_values(eaf_capacity)
            eaf_capacity_df[year] = eaf_capacity.reindex(eaf_capacity_df.index).fillna(0)
            
            # BOF capacity
            bof_efficiency = n.links.efficiency.loc[n.links.index.str.contains('BOF', case=False)].iloc[0]
            bof_capacity = n.links.p_nom_opt[n.links.index.str.contains("BOF")] * 8760 * bof_efficiency  # kton steel / year
            bof_capacity = country_values(bof_capacity)
            bof_capacity_df[year] = bof_capacity.reindex(bof_capacity_df.index).fillna(0)
            
            # H2 EAF
            h2_eaf = n.links_t.p0.filter(like="H2 to DRI", axis=1).sum() * timestep  # MWh H2 / year
            h2_eaf = country_values(h2_eaf)
            h2_eaf_df[year] = h2_eaf.reindex(h2_eaf_df.index).fillna(0)
            
            # NG EAF
            ng_eaf = n.links_t.p0.filter(like="CH4 to DRI", axis=1).sum() * timestep  # MWh CH4 / year
            ng_eaf = country_values(ng_eaf)
            ng_eaf_df[year] = ng_eaf.reindex(ng_eaf_df.index).fillna(0)
            
            # Steel TGR capacity
            tgr_capacity = n.links.p_nom_opt[n.links.index.str.contains("steel bof CC")] * 8760  # kton steel / year
            tgr_capacity = country_values(tgr_capacity)
            tgr_capacity_df[year] = tgr_capacity.reindex(tgr_capacity_df.index).fillna(0)
            
            # EAF production
            eaf_prod = -n.links_t.p1.filter(like="EAF", axis=1).sum() * timestep  # kt steel / year
            eaf_prod = country_values(eaf_prod)
            eaf_prod_df[year] = eaf_prod.reindex(eaf_prod_df.index).fillna(0)
            
            # BOF production
            bof_prod = -n.links_t.p1.filter(like="BOF", axis=1).sum() * timestep  # kt steel / year
            bof_prod = country_values(bof_prod)
            bof_prod_df[year] = bof_prod.reindex(bof_prod_df.index).fillna(0)
            
            # Steel load
            steel_load = n.loads_t.p.filter(like="steel", axis=1).sum()*timestep #kton steel per year per node
            steel_load = country_values(steel_load)
            steel_load_df[year] = steel_load.reindex(steel_load_df.index).fillna(0)
            
            # Electricity price     
            mprice = n.buses_t.marginal_price
            elec_prices = mprice[[col for col in mprice.columns if col.endswith(" 0")]]
    
            elec_loads_resi = retrieve_loads(n.loads_t.p, "0")
            elec_loads_ind = retrieve_loads(n.loads_t.p, "industry electricity") 
            elec_loads_agri = retrieve_loads(n.loads_t.p, "agriculture electricity")
            elec_loads_tra = retrieve_loads(n.loads_t.p, "EV")
            elec_links_th = retrieve_links(n.links_t.p0, "thermal|heat")
            elec_links_dac = retrieve_links(n.links_t.p1, "DAC")
            elec_links_meth = retrieve_links(n.links_t.p2, "methanolisation")
    
            elec_loads = (elec_loads_resi + elec_loads_ind + elec_loads_agri + elec_loads_tra + elec_links_th + elec_links_dac + elec_links_meth)
    
            # Get the average of prices per node based on power load
            total_elec_expen = elec_prices * elec_loads
            total_elec_expen = country_values(total_elec_expen.sum())
            elec_loads = country_values(elec_loads.sum())
            weighted_average = total_elec_expen / elec_loads  # €/MWh
            elec_prices_df[year] = weighted_average.reindex(elec_prices_df.index).fillna(0)
            
            # Ammonia production
            ammonia_prod = -n.links_t.p1.filter(like="Haber", axis=1).sum() * timestep  # kt NH3 / year
            ammonia_prod = country_values(ammonia_prod)
            ammonia_prod_df[year] = ammonia_prod.reindex(ammonia_prod_df.index).fillna(0)
            
            # Ammonia load
            ammonia_load = n.loads_t.p.filter(like="NH3", axis=1).sum() * timestep  # kt NH3 / year
            ammonia_load = country_values(ammonia_load)
            ammonia_load_df[year] = ammonia_load.reindex(ammonia_load_df.index).fillna(0)
            
            # H2 production (electrolysis, SMR, SMR + CC, ammonia cracker)
            # Electrolysis
            elec_h2prod = -n.links_t.p1.filter(like="Electrolysis", axis=1).sum() * timestep
            elec_h2prod = country_values(elec_h2prod)
            elec_h2prod_df[year] = elec_h2prod.reindex(elec_h2prod_df.index).fillna(0)
            
            # SMR + CC
            smrcc_h2prod = -n.links_t.p1.filter(like="SMR CC", axis=1).sum() * timestep
            smrcc_h2prod = country_values(smrcc_h2prod)
            smrcc_h2prod_df[year] = smrcc_h2prod.reindex(smrcc_h2prod_df.index).fillna(0)
            
            # SMR
            smr_h2prod = -n.links_t.p1.loc[:, n.links_t.p1.columns.str.contains("SMR") & ~n.links_t.p1.columns.str.contains("SMR CC")].sum() * timestep
            smr_h2prod = country_values(smr_h2prod)
            smr_h2prod_df[year] = smr_h2prod.reindex(smr_h2prod_df.index).fillna(0)
            
            # Ammonia
            nh3crack_h2prod = -n.links_t.p1.filter(like="ammonia cracker", axis=1).sum() * timestep
            nh3crack_h2prod = country_values(nh3crack_h2prod)
            nh3crack_h2prod_df[year] = nh3crack_h2prod.reindex(nh3crack_h2prod_df.index).fillna(0)
            
            # Emissions
            dri_emis = -n.links_t.p1.filter(like="dri process", axis=1).sum() * timestep / 1e3 # from t to kt per year
            dri_emis = country_values(dri_emis)
            bof_emis = -n.links_t.p1.filter(like="bof process", axis=1).sum() * timestep / 1e3 # from t to kt per year
            bof_emis = country_values(bof_emis)
            bofcc_emis = -n.links_t.p1.filter(like="bof CC", axis=1).sum() * timestep / 1e3 # from t to kt per year
            bofcc_emis = country_values(bofcc_emis)
            steel_proc_emis = dri_emis + bof_emis + bofcc_emis
            steel_proc_emis_df[year] = steel_proc_emis.reindex(steel_proc_emis_df.index).fillna(0)
            
            bofcc_emis_capt = -n.links_t.p2.filter(like="bof CC", axis=1).sum() * timestep / 1e3 # from t to kt per year
            bofcc_emis_capt = country_values(bofcc_emis_capt)
            bofcc_emis_capt_df[year] = bofcc_emis_capt.reindex(bofcc_emis_capt_df.index).fillna(0)
            
            cement_emis = -n.links_t.p1.filter(like="cement process", axis=1).sum() * timestep / 1e3 # from t to kt per year
            cement_emis = country_values(cement_emis)
            cementcc_emis = -n.links_t.p1.filter(like="cement CC", axis=1).sum() * timestep / 1e3 # from t to kt per year
            cementcc_emis = country_values(cementcc_emis)
            cement_proc_emis = cement_emis + cementcc_emis
            cement_proc_emis_df[year] = cement_proc_emis.reindex(cement_proc_emis_df.index).fillna(0)
            
            
            
            
        # Write all DataFrames to separate sheets
        eaf_capacity_df.to_excel(writer, sheet_name="EAF_Capacity")
        bof_capacity_df.to_excel(writer, sheet_name="BOF_Capacity")
        h2_eaf_df.to_excel(writer, sheet_name="H2_EAF")
        ng_eaf_df.to_excel(writer, sheet_name="NG_EAF")
        tgr_capacity_df.to_excel(writer, sheet_name="TGR_Capacity")
        eaf_prod_df.to_excel(writer, sheet_name="EAF_Prod")
        bof_prod_df.to_excel(writer, sheet_name="BOF_Prod")
        steel_load_df.to_excel(writer, sheet_name="Steel_Load")
        elec_prices_df.to_excel(writer, sheet_name="Elec_Prices")
        ammonia_prod_df.to_excel(writer, sheet_name="Ammonia_Prod")
        ammonia_load_df.to_excel(writer, sheet_name="Ammonia_Load")
        
        # Hydrogen
        elec_h2prod_df.to_excel(writer, sheet_name="Elec_H2Prod")
        smrcc_h2prod_df.to_excel(writer, sheet_name="SMRCC_H2Prod")
        smr_h2prod_df.to_excel(writer, sheet_name="SMR_H2Prod")
        nh3crack_h2prod_df.to_excel(writer, sheet_name="NH3crack_H2Prod")
        
        # Emissions
        steel_proc_emis_df.to_excel(writer, sheet_name="Steel_Proc_Emis")
        bofcc_emis_capt_df.to_excel(writer, sheet_name="Steel_Capt_Emis")
        cement_proc_emis_df.to_excel(writer, sheet_name="Cement_Proc_Emis")
    
    print(f"Data saved to {output_filename}")



