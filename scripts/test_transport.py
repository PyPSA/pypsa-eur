#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 09:46:12 2024

@author: lisa
"""
import pypsa

car = ['land transport EV heavy', 
       'land transport fuel cell heavy', 'land transport oil heavy']
def get_load(n, car):
    l_i = n.links[n.links.carrier.isin(car) & ~n.links.p_nom_extendable].index
    load = n.loads_t.p_set.loc[:, n.loads_t.p_set.columns.str.contains("heavy")]
    p=n.links.loc[l_i, "p_nom"].mul(n.links_t.p_max_pu[l_i]).mul(n.links_t.efficiency[l_i])
    p_grouped = p.groupby(n.links.bus1, axis=1).sum()
    return l_i, load, p, p_grouped



n = pypsa.Network("/home/lisa/Documents/playground/pypsa-eur/results/test-new-transport/postnetworks/elec_s_37_lv1.0___2040.nc")
l_i, load_2040, p_2040, p_grouped_2040 = get_load(n, car)

m = pypsa.Network("/home/lisa/Documents/playground/pypsa-eur/results/test-new-transport/prenetworks-brownfield/elec_s_37_lv1.0___2050.nc")
l_i, load_2050, p_2050, p_grouped_2050 = get_load(m, car)
