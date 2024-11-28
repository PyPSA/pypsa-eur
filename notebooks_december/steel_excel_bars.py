# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:46:12 2024

@author: alice
"""

def load_projection(plotting_params):
    proj_kwargs = plotting_params.get("projection", dict(name="EqualEarth"))
    proj_func = getattr(ccrs, proj_kwargs.pop("name"))
    return proj_func(**proj_kwargs)

def steel_preprocessing(excel_dir, scenario, config):
    
    eaf_prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "EAF_Prod", index_col="Bus")
    bof_prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "BOF_Prod", index_col="Bus")
    ng_eaf_prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "NG_EAF", index_col="Bus")
    h2_eaf_prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "H2_EAF", index_col="Bus")
    
    eaf_prod = eaf_prod[2050]
    bof_prod = bof_prod[2050]
    ng_eaf_prod = ng_eaf_prod[2050]
    h2_eaf_prod = h2_eaf_prod[2050]
    
    ng_eaf_prod = ((ng_eaf_prod / (ng_eaf_prod + h2_eaf_prod)) * eaf_prod).fillna(0)
    h2_eaf_prod = ((h2_eaf_prod / (ng_eaf_prod + h2_eaf_prod)) * eaf_prod).fillna(0)
    
    steel_prod = pd.concat([bof_prod, ng_eaf_prod, h2_eaf_prod], axis=1)
    steel_prod.columns = ['BOF','NG EAF','H2 EAF']
    
    steel_prod = steel_prod/ 1e3
    steel_prod = steel_prod.where(steel_prod >= 1, 0)
    #steel_prod = steel_prod[steel_prod.sum(axis=1) > 1]
    
    row_sums = steel_prod.sum(axis=1)

    steel_shares = steel_prod.div(row_sums, axis=0).fillna(0) * 100
    
    
    return steel_prod, steel_shares


def plot_bar(n, excel_data, regions, year, title,i,j,ncol, ax=None):

    
    # Plot the stacked bar chart
       if ax is None:
           fig, ax = plt.subplots(figsize=(10, 6))
       
       # Plot the stacked bar chart
       excel_data.plot(kind='bar', stacked=True, ax=ax, legend=False, width=1.0)
    
       # Customize the plot
       #ax.set_title(title, fontsize=14)
       ax.set_xlabel('Country', fontsize=12)
       ax.set_ylabel('Percentage (%)', fontsize=12)
       #ax.legend(title='Categories', fontsize=10)
    
       # Display the plot
       plt.tight_layout()


# %%


root_dir = "C:/Users/alice/Desktop/CMCC/pypsa-adb-industry/"
scenario = "baseline/"
res_dir = "results/"
regions_fn = root_dir + "resources/" + scenario + "regions_onshore_base_s_39.geojson"

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
import yaml
import numpy as np

with open(
    root_dir + res_dir + "baseline_eu_dem/configs/config.base_s_39_lvopt___2050.yaml"
) as config_file:
    config = yaml.safe_load(config_file)

regions = gpd.read_file(regions_fn).set_index("name")
map_opts = config["plotting"]["map"]

if map_opts["boundaries"] is None:
    map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

config["plotting"]["projection"]["name"] = "EqualEarth"
proj = load_projection(config["plotting"])

year = 2050

excel_dir = "excels_24h/"
# Plotting the green share of steel production
title = "Steel prod [Mt/yr]"

scenarios = ["baseline_eu_dem", "policy_eu_dem", "baseline_regional_dem", "policy_regional_dem"]
cool_names = ["Baseline European", "Policy European", "Baseline Regional", "Policy Regional"]
ncol = int(len(scenarios)/2)

fig, axes = plt.subplots(
    ncol, ncol,
    figsize=(3 * ncol, 3 * ncol),
    constrained_layout=False,
    subplot_kw={"projection": proj},
    gridspec_kw={'width_ratios': [0.805] + [1] * (ncol - 1) }
)

fn = (root_dir + "results/" + scenarios[0] + "/postnetworks/base_s_39_lvopt___2050.nc")
n = pypsa.Network(fn)

j = 0  # row index
i = 0  # column index

for scenario in scenarios:
    
    steel_prod, steel_shares = steel_preprocessing(excel_dir, scenario, config)

    ax = axes[j, i]
    plot_bar(n, steel_shares,  regions, year, title, i,j,ncol, ax=ax)
    
    if i == 0:
        if j == 0:
            ax.annotate(
                "Centralized", xy=(-0.05, 0.5), xycoords='axes fraction',
                ha='center', va='center', fontsize=12, rotation=90
            )
            ax.annotate(
                "Baseline", xy=(0.5,1.05), xycoords='axes fraction',
                ha='center', va='center', fontsize=12, rotation=0
            )
        elif j == 1:
            ax.annotate(
                "Regional", xy=(-0.05, 0.5), xycoords='axes fraction',
                ha='center', va='center', fontsize=12, rotation=90
            )
    elif i == 1:
        if j == 0:
            ax.annotate(
                "Climate Policy", xy=(0.5, 1.05), xycoords='axes fraction',
                ha='center', va='center', fontsize=12, rotation=0
            )
        
    # If we reach the last column, reset i and increment j (move to next row)
    if i == ncol - 1:  
        i = 0
        j += 1
    else:        
        i += 1
        
            

plt.tight_layout()
fig.suptitle('Year 2050', fontsize=16, y=1.02)
plt.savefig("graphs/steel_shares_2050.png", bbox_inches="tight")
plt.show()

# %%

steel_prod, steel_shares = steel_preprocessing(excel_dir, scenario, config)





