# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:25:49 2024

@author: alice
"""

def load_projection(plotting_params):
    proj_kwargs = plotting_params.get("projection", dict(name="EqualEarth"))
    proj_func = getattr(ccrs, proj_kwargs.pop("name"))
    return proj_func(**proj_kwargs)

def read_year(excel_dir, scenario, sheet, year):
    df = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = sheet, index_col="Bus")
    df = df[year]/1e3 # from kt to Mt of steel
    
    return df

def steel_pie_preprocessing(excel_dir, scenario, config, year):
    
    eaf_prod = read_year(excel_dir, scenario, 'EAF_Prod', year)
    bof_prod = read_year(excel_dir, scenario, "BOF_Prod", year)
    steel_prod = eaf_prod + bof_prod
    
    bof_cc_prod = read_year(excel_dir, scenario, "BOF_CC_Prod", year)
    bof_prod = bof_prod - bof_cc_prod
    ng_eaf_prod = read_year(excel_dir, scenario, "NG_EAF", year)
    h2_eaf_prod = read_year(excel_dir, scenario, "H2_EAF", year)
    
    
    # Hydrogen
    elec_h2prod = read_year(excel_dir, scenario, "Elec_H2Prod", year)
    smrcc_h2prod = read_year(excel_dir, scenario, "SMRCC_H2Prod", year)
    smr_h2prod = read_year(excel_dir, scenario, "SMR_H2Prod", year)
    nh3crack_h2prod = read_year(excel_dir, scenario, "NH3crack_H2Prod", year)
    
    green_h2_share = ((elec_h2prod + smrcc_h2prod) / ( elec_h2prod + smrcc_h2prod + smr_h2prod + nh3crack_h2prod)).fillna(0)
    green_h2_eaf_prod = green_h2_share * h2_eaf_prod
    grey_h2_eaf_prod = (1-green_h2_eaf_prod) * h2_eaf_prod

    steel_tech = pd.concat([bof_prod, ng_eaf_prod, grey_h2_eaf_prod, bof_cc_prod, green_h2_eaf_prod ], axis = 1)
    steel_tech.columns = ['BOF','NG EAF', 'Grey H2 EAF','BOF CC','Green H2 EAF']
    
    return steel_prod, steel_tech


def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue
            names = ifind.index[ifind == i]
            c.df.loc[names, "location"] = names.str[:i]


def plot_pie_map(n, steel_prod, steel_tech, regions, year, i, j, ax=None): 
    
    assign_location(n)
    
    regions["data"] = regions.index.str[:2].map(steel_prod) 
    regions["data"] = regions["data"].where(regions["data"] > 0.1, 0)
    regions["data"] = regions["data"].fillna(0)
    
    
    row_sums = steel_tech.sum(axis=1)
    steel_prod_shares = round(steel_tech.div(row_sums, axis=0),1)

    regions = regions.to_crs(proj.proj4_init)
    max_value = steel_prod.values.max()
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": proj})
    
       
    regions.plot(
        ax=ax,
        column="data",
        cmap="Blues",
        linewidths=0.5,
        edgecolor='black',
        legend=(i == ncol - 1),
        vmax=max_value,
        vmin=0,
        legend_kwds={
            "label": "Steel production [Mt steel/yr]",
            #"fontsize": 12,
            "shrink": 0.5,
            "extend": "max",
        },
    )
    
    # Adjust the colors for the five items
    piecolors = ['#000000', '#636161', '#A09F9F', '#090C9B', '#42CF05']
    idx_prefix = None
    
    for idx, region in regions.iterrows():
        # Check if the first two characters of the current idx match the last processed one
        if idx[:2] == idx_prefix:
            idx_prefix = idx[:2]
            continue
        idx_prefix = idx[:2]
        # Get the centroid of each region
        centroid = region['geometry'].centroid
    
        if idx_prefix in steel_prod_shares.index:
            pie_values = [
                steel_prod_shares.loc[idx_prefix, 'BOF'],
                steel_prod_shares.loc[idx_prefix, 'NG EAF'],
                steel_prod_shares.loc[idx_prefix, 'Grey H2 EAF'],
                steel_prod_shares.loc[idx_prefix, 'BOF CC'],
                steel_prod_shares.loc[idx_prefix, 'Green H2 EAF']
            ]
    
            # If all values are 0, skip plotting
            if sum(pie_values) == 0:
                continue
    
            # Normalize the size for the pie chart
            size = 100  # Adjust size scaling as needed
    
            # If one value is 1 (or 100% share), plot a full circle
            if max(pie_values) == 1:
                full_circle_color = piecolors[pie_values.index(max(pie_values))]
                ax.scatter(
                    centroid.x, centroid.y,
                    s=size,
                    color=full_circle_color,
                    edgecolor='black',  # Black border for the full circle
                    linewidth=0.3,
                    marker='o'  # Full circle
                )
                continue
    
            # Initialize the starting angle for the pie slices
            previous = 0
            markers = []
    
            # Generate pie chart segments as markers
            for color, ratio in zip(piecolors, pie_values):
                if ratio == 0:
                    continue
    
                this = 2 * np.pi * ratio + previous
                x = [0] + np.cos(np.linspace(previous, this, 10)).tolist() + [0]
                y = [0] + np.sin(np.linspace(previous, this, 10)).tolist() + [0]
                xy = np.column_stack([x, y])
                previous = this
    
                # Append marker settings for this pie slice
                markers.append({
                    'marker': xy,
                    's': size,
                    'facecolor': color,
                    'edgecolor': 'black',
                    'linewidth': 0.3
                })
    
            # Scatter each pie slice
            for marker in markers:
                ax.scatter(centroid.x, centroid.y, **marker)

            
    ax.set_facecolor("white")
    
    # Add a title and subtitle (if provided)
    ax.set_title(year, fontsize=14, loc="center")  # Main title

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

years = [2030, 2040, 2050]

excel_dir = "excels_24h/"
# Plotting the green share of steel production
title = "Steel prod [Mt/yr]"

scenarios = ["baseline_eu_dem", "policy_eu_dem", "baseline_regional_dem", "policy_regional_dem"]
cool_names = ["Baseline European", "Policy European", "Baseline Regional", "Policy Regional"]
nrows = int(len(scenarios))
ncol = int(len(years))

fig, axes = plt.subplots(
    nrows, ncol,
    figsize=(3 * nrows, 3 * ncol),
    constrained_layout=False,
    subplot_kw={"projection": proj},
    gridspec_kw={'width_ratios': [0.805] + [1] * (ncol - 1) }
)

fn = (root_dir + "results/" + scenarios[0] + "/postnetworks/base_s_39_lvopt___2050.nc")
n = pypsa.Network(fn)

j = 0  # row index
i = 0  # column index

for j, scenario in enumerate(scenarios):
    for i, year in enumerate(years):
        
        steel_prod, steel_shares = steel_pie_preprocessing(excel_dir, scenario, config, year)
    
        ax = axes[j, i]
        plot_pie_map(n, steel_prod, steel_shares, regions, year, i, j, ax=ax)
        
        if i == 0:
            if j == 0:
                ax.annotate(
                    "BASELINE\nEU DEM", xy=(-0.1, 0.5), xycoords='axes fraction',
                    ha='center', va='center', fontsize=12, rotation=90
                )
            elif j == 1:
                ax.annotate(
                    "POLICY\nEU DEM", xy=(-0.1, 0.5), xycoords='axes fraction',
                    ha='center', va='center', fontsize=12, rotation=90
                )
            elif j == 2:
                ax.annotate(
                    "BASELINE\nREG. DEM", xy=(-0.1, 0.5), xycoords='axes fraction',
                    ha='center', va='center', fontsize=12, rotation=90
                )
            elif j == 3:
                ax.annotate(
                    "POLICY\nREG. DEM", xy=(-0.1, 0.5), xycoords='axes fraction',
                    ha='center', va='center', fontsize=12, rotation=90
                )
   
            
            
plt.tight_layout()
fig.suptitle('Year 2050', fontsize=16, y=1.02)
plt.savefig("graphs/steel_pie_charts_years.png", bbox_inches="tight")
plt.show()

