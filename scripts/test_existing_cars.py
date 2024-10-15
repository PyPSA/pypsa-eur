#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:08:27 2024

@author: lisa
"""
changed_d = {}
registrations = pd.read_csv(snakemake.input.car_registration, index_col=[0,1])
for year in np.arange(2025,2055, 5):
    
    transport_type="heavy"
    filter_links = (~n.links.p_nom_extendable & (n.links.lifetime==np.inf)
                    & (n.links.carrier == f"land transport oil {transport_type}"))
    links_i = n.links[filter_links].index
    
    factor = options["car_reg_factor"]
    reg = registrations.loc[transport_type].iloc[:,0] * factor
    
    unchanged_fleet = (1-(reg*(year-ref_year))).clip(lower=0)
    previous_year = n_p.links.build_year.max()
    already_reduced =  (1-(reg*(previous_year-ref_year))).clip(lower=0)
    changed_d[year] = (unchanged_fleet/already_reduced.replace(0,1)).rename(index= lambda x: x + f" land transport oil {transport_type}-existing")