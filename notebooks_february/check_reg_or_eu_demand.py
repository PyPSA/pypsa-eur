# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 12:14:55 2025

@author: Dibella
"""


import pypsa
import matplotlib.pyplot as plt
import os
import pandas as pd

# Cement

scenarios = ["base_eu_regain", "policy_eu_regain", "base_eu_deindustrial","policy_eu_deindustrial",
             "base_reg_regain", "policy_reg_regain", "base_reg_deindustrial","policy_reg_deindustrial"]
# Replace ISO2 country codes with full names
country_names = {
    "AL": "Albania", "AT": "Austria", "BA": "Bosnia", "BE": "Belgium", "BG": "Bulgaria",
    "CH": "Switzerland", "CZ": "Czechia", "DE": "Germany", "DK": "Denmark", "EE": "Estonia", "ES": "Spain",
    "FI": "Finland", "FR": "France", "GB": "UK", "GR": "Greece", "HR": "Croatia", "HU": "Hungary",
    "IE": "Ireland", "IT": "Italy", "LT": "Lithuania", "LU": "Luxembourg", "LV": "Latvia", "ME": "Montenegro",
    "MK": "North Macedonia", "NL": "Netherlands", "NO": "Norway", "PL": "Poland", "PT": "Portugal", "RO": "Romania",
    "RS": "Serbia", "SE": "Sweden", "SI": "Slovenia", "SK": "Slovakia", "XK": "Kosovo"
}
#scenarios = ["base_eu_regain", "policy_eu_regain"]

commodities = ['steel','cement','NH3','HVC','industry methanol']
df = pd.DataFrame(index=commodities,columns = scenarios)

# First pass: Collect data across all scenarios
for scenario in scenarios:
    cwd = os.getcwd()
    parent_dir = os.path.dirname(cwd)
    file_path = os.path.join(parent_dir, "results_march", scenario, "networks", "base_s_39___2050.nc")
    n = pypsa.Network(file_path)
    timestep = n.snapshot_weightings.iloc[0, 0]
    
    loads = n.loads
    
    # Check commodities
    for com in commodities:
        com_load = loads[loads.index.str.contains(com)]
        df.loc[com,scenario] = len(com_load)
    
