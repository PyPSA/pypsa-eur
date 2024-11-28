# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 14:55:16 2024

@author: alice
"""

import pypsa
import pandas as pd


root_dir = "C:/Users/alice/Desktop/CMCC/pypsa-adb-industry/"
years = [2030, 2040, 2050]
scenarios = ["climate_policy_eu_dem", "climate_policy_regional_dem"]
n_eu = dict()
n_reg = dict()

costs = pd.DataFrame(0, index = years, columns = scenarios)

for scenario in scenarios:
    for year in years:
        
        fn = (root_dir + "results/" + scenario + f"/postnetworks/base_s_39_lvopt___{year}.nc")
        if scenario == scenarios[0]:
            n_eu[year] = pypsa.Network(fn)
            costs.loc[year,scenario] = n_eu[year].objective/1e6 # M€
        else:
            n_reg[year] = pypsa.Network(fn)
            costs.loc[year,scenario] = n_reg[year].objective/1e6 # M€
        
    
    
share_change = (costs["climate_policy_regional_dem"] / costs["climate_policy_eu_dem"] -1 ) *100

# %% Let's try to get costs by country

