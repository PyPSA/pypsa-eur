# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 08:15:43 2025

@author: Dibella
"""


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
            
def share_green_h2(n):
    # Extract H2 Clean and H2 Dirty production
    timestep = n.snapshot_weightings.iloc[0,0]
    h2_clean = n.links.loc[n.links.index.str.contains('H2 Electrolysis|SMR CC', regex=True, na=False), :].index
    h2_dirty = n.links.loc[n.links.index.str.contains('SMR(?! CC)', regex=True, na=False), :].index
    h2_clean_df = -n.links_t.p1.loc[:, h2_clean].sum() * timestep
    h2_dirty_df = -n.links_t.p1.loc[:, h2_dirty].sum() * timestep
    
    h2_clean_df.index = h2_clean_df.index.str[:5]
    h2_clean_df = h2_clean_df.groupby(h2_clean_df.index).sum()
    
    h2_dirty_df.index = h2_dirty_df.index.str[:5]
    h2_dirty_df = h2_dirty_df.groupby(h2_dirty_df.index).sum()
    
    # Calculate share of green and grey H2
    share_green = round(h2_clean_df / (h2_clean_df + h2_dirty_df), 2)
    
    return share_green

def plot_steel_map(n, regions, year,i, ax=None): 
    
    assign_location(n)
    timestep = n.snapshot_weightings.iloc[0,0]
    
    steel_prod_index = n.links[n.links['bus1'].str.contains('steel', case=False, na=False) &
    ~n.links['bus1'].str.contains('heat', case=False, na=False)].index
    steel_prod = -n.links_t.p1.loc[:,steel_prod_index].sum()#*timestep/1e3 #Mtsteel
    steel_prod.index = steel_prod.index.str.split(' 0 ').str[0] + ' 0'
    steel_prod_df = steel_prod.to_frame(name='steel_prod')
    
    steel_prod_tech = -n.links_t.p1.loc[:,steel_prod_index].sum()#*timestep/1e3 #Mtsteel
    
    steel_prod_eaf = steel_prod_tech[steel_prod_tech.index.str.contains('EAF')]
    steel_prod_eaf.index = steel_prod_eaf.index.str.split(' 0 ').str[0] + ' 0'
    steel_prod_eaf = steel_prod_eaf.groupby(level=0).sum().div(1e3)*timestep
    steel_prod_eaf = steel_prod_eaf.where(steel_prod_eaf > 0.1, 0)
    
    steel_prod_bof = steel_prod_tech[steel_prod_tech.index.str.contains('BOF')]
    steel_prod_bof.index = steel_prod_bof.index.str.split(' 0 ').str[0] + ' 0'
    steel_prod_bof = steel_prod_bof.groupby(level=0).sum().div(1e3)*timestep
    steel_prod_bof = steel_prod_bof.where(steel_prod_bof > 0.1, 0)
    
    # Extract H2 Clean and H2 Dirty production
    share_green = share_green_h2(n)
    #share_green = share_green.reindex(stee.index, fill_value=0)

    # Extract CH4 and H2-based DRI data
    dri_ch4 = -n.links_t.p1.filter(like='CH4 to syn gas DRI', axis=1).sum() * timestep
    dri_h2 = -n.links_t.p1.filter(like='H2 to syn gas DRI', axis=1).sum() * timestep

    # Calculate the share of H2 in DRI production -> at the European level now
    share_h2 = round(dri_h2.sum() / (dri_h2.sum() + dri_ch4.sum()), 2)
    
    # Calculate the adjusted EAF production
    steel_prod_ch4_eaf = steel_prod_eaf * (1 - share_h2)  # EAF with CH4
    steel_prod_grey_h2_eaf = steel_prod_eaf * share_h2 * (1 - share_green)  # EAF with Grey H2
    steel_prod_green_h2_eaf = steel_prod_eaf * share_h2 * share_green  # EAF with Green H2
    
    
    steel_prod_tech = pd.concat([steel_prod_green_h2_eaf,steel_prod_grey_h2_eaf, steel_prod_ch4_eaf, steel_prod_bof], axis=1)

    #steel_prod_tech = pd.concat([steel_prod_eaf, steel_prod_bof], axis=1)
    steel_prod_tech.columns = ['Green H2 EAF','Grey H2 EAF','CH4 EAF','BOF']
    row_sums = steel_prod_tech.sum(axis=1)
    steel_prod_shares = steel_prod_tech.div(row_sums, axis=0).fillna(0)
    
    
    regions["steel"] = (
        steel_prod_df
        .steel_prod.groupby(level=0)
        .sum()
        .div(1e3)*timestep
    )  # TWh
    
    
    regions["steel"] = regions["steel"].where(regions["steel"] > 0.1, 0)
    regions["steel"] = regions["steel"].fillna(0)
    
    
    regions = regions.to_crs(proj.proj4_init)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": proj})
    
    
    #show_legend = i == 0  # Only show legend in column 0
    show_legend = 1
    regions.plot(
        ax=ax,
        column="steel",
        cmap="Blues",
        linewidths=0.3,
        legend=show_legend,
        vmax=50,
        vmin=0,
        edgecolor="black",
        legend_kwds={
            "label": "Steel production [Mt steel/yr]",
            #"fontsize": 10,
            "shrink": 0.7,
            "extend": "max",
        } if show_legend else {},
    )
    """
    regions.plot(
        ax=ax,
        column="steel",
        cmap="Blues",
        linewidths=0.3,
        legend=True,
        vmax=100,
        vmin=0,
        edgecolor="black",  # Add black outline
        legend_kwds={
            "label": "Steel production [Mt steel/yr]",
            #"fontsize": 12,
            "shrink": 0.5,
            "extend": "max",
        },
    )
    """
   # piecolors=[ '#5DCB4C', '#94958C']
    piecolors=['green','gray', '#552C2D', 'black']

    for idx, region in regions.iterrows():
        # Get the centroid of each region
        centroid = region['geometry'].centroid
        
        # Get the EAF and BOF shares for this region
        if idx in steel_prod_shares.index:
            green_eaf_percent = steel_prod_shares.loc[idx, 'Green H2 EAF']
            grey_eaf_percent = steel_prod_shares.loc[idx, 'Grey H2 EAF']
            ch4_eaf_percent = steel_prod_shares.loc[idx, 'CH4 EAF']
            bof_percent = steel_prod_shares.loc[idx, 'BOF']
            pie_values = [green_eaf_percent, grey_eaf_percent, ch4_eaf_percent, bof_percent]
            
            # If both values are 0, skip plotting
            if sum(pie_values) == 0:
                continue
            
            # Normalize the size for the pie chart
            size = 100 #sizes.get(idx, 1) / 15  # Scale size if needed
            
            # If any ratio is 1, plot a full circle
            # Check for full dominance
            if green_eaf_percent == 1:
                full_circle_color = piecolors[0]
            elif grey_eaf_percent == 1:
                full_circle_color = piecolors[1]
            elif ch4_eaf_percent == 1:
                full_circle_color = piecolors[2]
            elif bof_percent == 1:
                full_circle_color = piecolors[3]
            else:
                full_circle_color = None  # fallback, not used here
                
            if full_circle_color:
                
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
            
    if year == 2050 and scenario == 'base_eu_regain':
        legend_elements = [
            Patch(facecolor=piecolors[0], edgecolor='black', label='Green H2 EAF'),
            Patch(facecolor=piecolors[1], edgecolor='black', label='Grey H2 EAF'),
            Patch(facecolor=piecolors[2], edgecolor='black', label='CH4 EAF'),
            Patch(facecolor=piecolors[3], edgecolor='black', label='BOF'),
        ]
        ax.legend(
            handles=legend_elements,
            title="Steel tech share",
            loc='center left',
            bbox_to_anchor=(1.35, 0.5),
            frameon=False
        )
    ax.set_facecolor("white")
    
    # Add a title and subtitle (if provided)
    ax.set_title(year, fontsize=14, loc="center")  # Main title

    
#%%


root_dir = "C:/Users/Dibella/Desktop/CMCC/pypsa-adb-industry/"
scenario = "base_eu_regain/"
res_dir = "results_march/"
regions_fn = root_dir + "resources/" + scenario + "regions_onshore_base_s_39.geojson"

import pypsa
import pandas as pd
import numpy as np
import geopandas as gpd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import yaml
with open(r"C:\Users\Dibella\Desktop\CMCC\pypsa-adb-industry\results_march\base_eu_regain\configs\config.base_s_39___2030.yaml") as config_file: config = yaml.safe_load(config_file)
        #root_dir + res_dir + "base_eu_regain/configs/config.base_s_39_lvopt___2030.yaml") 
    
regions = gpd.read_file(regions_fn).set_index("name")
map_opts = config["plotting"]["map"]

if map_opts["boundaries"] is None:
    map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

config["plotting"]["projection"]["name"] = "EqualEarth"
proj = load_projection(config["plotting"])

years = [2030, 2040, 2050]
scenarios = ["base_eu_regain", "policy_eu_regain"]

fig, axes = plt.subplots(len(scenarios), len(years), figsize=(3*len(years), 3*len(scenarios)), subplot_kw={"projection": proj})


for i, year in enumerate(years):
    for j, scenario in enumerate(scenarios):
        fn = root_dir + "results_march/" + scenario + f"/networks/base_s_39___{year}.nc"
        ax = axes[j, i]
        n = pypsa.Network(fn)
        plot_steel_map(n, regions, year,i, ax=ax)

plt.tight_layout()

plt.savefig("graphs/steel_prod_per_country_pie_chart_eu_dem.png", bbox_inches='tight')
plt.show()
