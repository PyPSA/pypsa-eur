# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:39:20 2024

@author: alice
"""

import pypsa
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import numpy as np

root_dir = "C:/Users/alice/Desktop/CMCC/pypsa-adb-industry/"
scenario = "baseline/"
res_dir = "results/"
regions_fn = root_dir + "resources/" + scenario + "regions_onshore_base_s_39.geojson"


with open(
    root_dir + res_dir + scenario + "configs/config.base_s_39_lvopt___2030.yaml"
) as config_file:
    config = yaml.safe_load(config_file)

def country_values(df):
    df.index = df.index.str[:2]
    grouped = df.groupby(df.index).sum()
    grouped = grouped[abs(grouped) > 1e-4]
    return grouped

def nodal_values(df):
    df.index = df.index.str[:5]
    grouped = df.groupby(df.index).sum()
    grouped = grouped[abs(grouped) > 1e-2]
    return grouped

def load_projection(plotting_params):
    proj_kwargs = plotting_params.get("projection", dict(name="EqualEarth"))
    proj_func = getattr(ccrs, proj_kwargs.pop("name"))
    return proj_func(**proj_kwargs)


def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue
            names = ifind.index[ifind == i]
            c.df.loc[names, "location"] = names.str[:i]


def plot_steel_map(n, regions, year, ax=None): 
    
    assign_location(n)
    timestep = n.snapshot_weightings.iloc[0,0]
    
    steel_load = n.loads_t.p.filter(like="steel", axis=1).sum()*timestep #kton steel per year per node
    steel_load = nodal_values(steel_load)
    
    steel_supply_bof = -n.links_t.p1.filter(like="BOF", axis=1).sum() * timestep  # kton steel per year per node
    steel_supply_bof = nodal_values(steel_supply_bof)
    
    steel_supply_h2_eaf = -n.links_t.p1.filter(like="H2-DRI-EAF", axis=1).sum() * timestep  # kton steel per year per node
    steel_supply_h2_eaf = nodal_values(steel_supply_h2_eaf)
    
    steel_supply_ng_eaf = -n.links_t.p1.filter(like="NG-DRI-EAF", axis=1).sum() * timestep  # kton steel per year per node
    steel_supply_ng_eaf = nodal_values(steel_supply_ng_eaf)
    
    steel_prod_tech = pd.concat([steel_supply_h2_eaf, steel_supply_ng_eaf, steel_supply_bof], axis=1)
    steel_prod_tech.columns = ['H2-EAF','NG-EAF','BOF']
    steel_prod_tech[steel_prod_tech < 0.001] = 0
    row_sums = steel_prod_tech.sum(axis=1)
    steel_prod_shares = steel_prod_tech.div(row_sums, axis=0).fillna(0)
    
    steel_prod_shares = (steel_prod_shares * 100).round()
    round(row_sums.sum()) == round(steel_load.sum())
     
    regions["steel"] = (row_sums.div(1e3)) #Mton
    
    regions["steel"] = regions["steel"].where(regions["steel"] > 0.1, 0)
    regions["steel"] = regions["steel"].fillna(0)
    
    
    regions = regions.to_crs(proj.proj4_init)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": proj})
    
    
        
    regions.plot(
        ax=ax,
        column="steel",
        cmap="Blues",
        linewidths=0,
        legend=True,
        vmax=10,
        vmin=0,
        legend_kwds={
            "label": "Steel production [Mt steel/yr]",
            #"fontsize": 12,
            "shrink": 0.5,
            "extend": "max",
        },
    )
    
    piecolors = ['#5DCB4C', '#94958C', '#39393A']
    
    for idx, region in regions.iterrows():
        # Get the centroid of each region
        centroid = region['geometry'].centroid
        
        # Get the EAF and BOF shares for this region
        if idx in steel_prod_shares.index:
            h2_eaf_percent = steel_prod_shares.loc[idx, 'H2-EAF']
            ng_eaf_percent = steel_prod_shares.loc[idx, 'NG-EAF']
            bof_percent = steel_prod_shares.loc[idx, 'BOF']
            pie_values = [h2_eaf_percent, ng_eaf_percent, bof_percent]
            
            # If all values are 0, skip plotting
            if sum(pie_values) == 0:
                continue
            
            # Normalize the size for the pie chart
            size = 100  # Adjust size scaling as necessary
            
            # If any ratio is 1, plot a full circle
            if h2_eaf_percent == 1 or ng_eaf_percent == 1 or bof_percent == 1:
                # Determine the color based on which percentage is 1
                if h2_eaf_percent == 1:
                    full_circle_color = piecolors[0]  
                elif ng_eaf_percent == 1:
                    full_circle_color = piecolors[1]  
                else:
                    full_circle_color = piecolors[2]  
                
                ax.scatter(
                    centroid.x, centroid.y,
                    s=size,
                    color=full_circle_color,
                    edgecolor='black',  # Black border for the full circle
                    linewidth=0.5,
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
                markers.append({'marker': xy, 's': size, 'facecolor': color, 'edgecolor': 'black', 'linewidth': 0.5})
            
            # Scatter each pie slice
            for marker in markers:
                ax.scatter(centroid.x, centroid.y, **marker)
    
    ax.set_facecolor("white")
    
    # Add a title and subtitle (if provided)
    ax.set_title(year, fontsize=14, loc="center")  # Main title


# %%

regions = gpd.read_file(regions_fn).set_index("name")
map_opts = config["plotting"]["map"]

if map_opts["boundaries"] is None:
    map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

config["plotting"]["projection"]["name"] = "EqualEarth"
proj = load_projection(config["plotting"])

years = [2030, 2040, 2050]
scenarios = ["baseline_eu_dem", "climate_policy_eu_dem"]
fig, axes = plt.subplots(len(scenarios), len(years), figsize=(3*len(years), 3*len(scenarios)), subplot_kw={"projection": proj})


for i, year in enumerate(years):
    for j, scenario in enumerate(scenarios):
        fn = root_dir + "results/" + scenario + f"/postnetworks/base_s_39_lvopt___{year}.nc"
        ax = axes[j, i]
        n = pypsa.Network(fn)
        plot_steel_map(n, regions, year, ax=ax)

plt.tight_layout()

plt.savefig("graphs/steel_prod_per_country_pie_chart_new.png", bbox_inches='tight')
plt.show()


