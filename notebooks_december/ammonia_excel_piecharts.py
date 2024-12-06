# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:03:24 2024

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

def ammonia_pie_preprocessing(excel_dir, scenario, config, year):
    
    ammonia_prod = read_year(excel_dir, scenario, 'Ammonia_Prod', year)
    nh3_lhv = 5.3 #MWh/t ammonia
    ammonia_prod = ammonia_prod / nh3_lhv / 1e3 # kt NH3
    # Hydrogen
    elec_h2prod = read_year(excel_dir, scenario, "Elec_H2Prod", year)
    smrcc_h2prod = read_year(excel_dir, scenario, "SMRCC_H2Prod", year)
    smr_h2prod = read_year(excel_dir, scenario, "SMR_H2Prod", year)
    nh3crack_h2prod = read_year(excel_dir, scenario, "NH3crack_H2Prod", year)
    
    green_h2_share = ((elec_h2prod + smrcc_h2prod) / ( elec_h2prod + smrcc_h2prod + smr_h2prod + nh3crack_h2prod)).fillna(0)
    green_hb_share = green_h2_share * ammonia_prod
    nongreen_hb_share = (1- green_h2_share) * ammonia_prod

    ammonia_tech = pd.concat([green_hb_share, nongreen_hb_share], axis=1)
    ammonia_tech.columns = ['Green H2 HB', 'Grey H2 HB']
    ammonia_tech = ammonia_tech[ammonia_tech.sum(axis=1) >= 1]
    
    return ammonia_prod, ammonia_tech


def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue
            names = ifind.index[ifind == i]
            c.df.loc[names, "location"] = names.str[:i]


def plot_pie_map(n, ammonia_prod, ammonia_tech, regions, year, i, j, piecolors, ax=None): 
    
    assign_location(n)
    
    regions["data"] = regions.index.str[:2].map(ammonia_prod) 
    regions["data"] = regions["data"].where(regions["data"] > 0.1, 0)
    regions["data"] = regions["data"].fillna(0)
    
    
    row_sums = ammonia_tech.sum(axis=1)
    ammonia_prod_shares = round(ammonia_tech.div(row_sums, axis=0),2)

    regions = regions.to_crs(proj.proj4_init)
    #max_value = steel_prod.values.max()
    max_value=10
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": proj})
    
       
    regions.plot(
        ax=ax,
        column="data",
        cmap="Reds",
        linewidths=0.5,
        edgecolor='black',
        legend=(i == ncol - 1),
        vmax=max_value,
        vmin=0,
        legend_kwds={
            "label": "Ammonia prod\n[kt NH3/yr]",
            #"fontsize": 12,
            "shrink": 0.5,
            "extend": "max",
        },
    )
    
    # Adjust the colors for the five items
    idx_prefix = None
    
    for idx, region in regions.iterrows():
        # Check if the first two characters of the current idx match the last processed one
        if idx[:2] == idx_prefix:
            idx_prefix = idx[:2]
            continue
        idx_prefix = idx[:2]
        # Get the centroid of each region
        centroid = region['geometry'].centroid
    
        if idx_prefix in ammonia_prod_shares.index:
            pie_values = [
                ammonia_prod_shares.loc[idx_prefix, 'Green H2 HB'],
                ammonia_prod_shares.loc[idx_prefix, 'Grey H2 HB'],
            ]
    
            # If all values are 0, skip plotting
            if sum(pie_values) == 0:
                continue
    
            # Normalize the size for the pie chart
            size = 75  # Adjust size scaling as needed
    
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
    if j == 0:
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
piecolors = ['#42CF05', '#636161']
pielabels = ['Green H2 HB', 'Grey H2 HB']

excel_dir = "excels_24h/"
# Plotting the green share of steel production
title = "Steel prod [Mt/yr]"

scenarios = ["baseline_eu_dem", "policy_eu_dem", "baseline_regional_dem", "policy_regional_dem"]
nrows = int(len(scenarios))
ncol = int(len(years))

fig, axes = plt.subplots(
    nrows, ncol,
    figsize=(2.5 * nrows, 2.5 * ncol),
    constrained_layout=False,
    subplot_kw={"projection": proj},
    #gridspec_kw={'width_ratios': [0.805] + [1] * (ncol - 1) }
)

fig.subplots_adjust(wspace=-0.40, hspace=0)

fn = (root_dir + "results/" + scenarios[0] + "/postnetworks/base_s_39_lvopt___2050.nc")
n = pypsa.Network(fn)

j = 0  # row index
i = 0  # column index

for j, scenario in enumerate(scenarios):
    for i, year in enumerate(years):
        
        ammonia_prod, ammonia_tech = ammonia_pie_preprocessing(excel_dir, scenario, config, year)
    
        ax = axes[j, i]
        plot_pie_map(n, ammonia_prod, ammonia_tech, regions, year, i, j, piecolors, ax=ax)
        ax.axis('off')
        
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

patches = [plt.Line2D([0], [0], marker='o', color=color, markersize=10, linestyle='None') for color in piecolors]
fig.legend(patches, pielabels, loc="lower center", ncol=len(piecolors), fontsize=10, frameon=False,
           bbox_to_anchor=(0.5, 0.05))

"""
fig.legend(
    handles=[
        plt.Line2D([0], [0], marker='o', color=color, markersize=10, linestyle='None')
        for color in piecolors
    ],
    labels=pielabels,
    loc='center right',  # Position legend on the right
    bbox_to_anchor=(1.1, 0.5),  # Adjust position (x, y)
    title='Steel Production Technologies',
    title_fontsize=12,
    fontsize=10
)
"""          
#plt.constrained_layout()
#fig.suptitle('Year 2050', fontsize=16, y=1.02)
plt.savefig("graphs/ammonia_pie_charts_years.png", bbox_inches="tight")
plt.show()

