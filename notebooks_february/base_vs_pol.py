# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 13:41:53 2025

@author: Dibella
"""

import pypsa
import matplotlib.pyplot as plt
import os
import pandas as pd


scenarios = ["base_eu_regain", "policy_eu_regain"]

# Replace ISO2 country codes with full names
country_names = {
    "AL": "Albania", "AT": "Austria", "BA": "Bosnia", "BE": "Belgium", "BG": "Bulgaria",
    "CH": "Switzerland", "CZ": "Czechia", "DE": "Germany", "DK": "Denmark", "EE": "Estonia", "ES": "Spain",
    "FI": "Finland", "FR": "France", "GB": "UK", "GR": "Greece", "HR": "Croatia", "HU": "Hungary",
    "IE": "Ireland", "IT": "Italy", "LT": "Lithuania", "LU": "Luxembourg", "LV": "Latvia", "ME": "Montenegro",
    "MK": "North Macedonia", "NL": "Netherlands", "NO": "Norway", "PL": "Poland", "PT": "Portugal", "RO": "Romania",
    "RS": "Serbia", "SE": "Sweden", "SI": "Slovenia", "SK": "Slovakia", "XK": "Kosovo"
}

commodities = ['steel','cement','NH3','HVC','industry methanol']
df = pd.DataFrame(index=commodities,columns = scenarios)
df = pd.DataFrame(columns=scenarios)

# First pass: Collect data across all scenarios
for scenario in scenarios:
    cwd = os.getcwd()
    parent_dir = os.path.dirname(cwd)
    file_path = os.path.join(parent_dir, "results", scenario, "networks", "base_s_39___2050.nc")
    n = pypsa.Network(file_path)
    timestep = n.snapshot_weightings.iloc[0, 0]
    df.loc['Annual system cost [b€/yr]',scenario] = n.objective/1e9
    capex = n.statistics.capex()
    opex = n.statistics.opex()
    threshold = 0.1
    capex = capex[capex>threshold]
    capex.index = capex.reset_index()
    
    
    capex = capex[~capex.index.str.contains('StorageUnit')]
    opex = opex[opex>threshold]
