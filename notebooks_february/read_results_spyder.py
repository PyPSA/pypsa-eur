# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:04:08 2025

@author: alice
"""

import pypsa

n = pypsa.Network(r"C:\Users\alice\Desktop\CMCC\pypsa-adb-industry\results\base_eu_regain\networks\base_s_39___2050.nc")

timestep = n.snapshot_weightings.iloc[0,0]
alinks = -n.links_t.p1.filter(regex='H2 to syn gas DRI', axis = 1).sum() * timestep

cement_not_captured = -n.links_t.p1.filter(like='cement process emis to atmosphere', axis=1).sum() * timestep
cement_ccs = -n.links_t.p1.filter(like='cement CC', axis=1).sum() * timestep

cement_not_captured.index = cement_not_captured.index.str[:2]
cement_not_captured = cement_not_captured.groupby(cement_not_captured.index).sum()

cement_ccs.index = cement_ccs.index.str[:2]
cement_ccs = cement_ccs.groupby(cement_ccs.index).sum()

# Calculate the share of green H2 and dirty H2
share_ccs = round(cement_ccs / (cement_ccs + cement_not_captured), 2)
