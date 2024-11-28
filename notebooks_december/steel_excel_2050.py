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
    
    # Hydrogen
    elec_h2prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "Elec_H2Prod", index_col="Bus")
    smrcc_h2prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "SMRCC_H2Prod", index_col="Bus")
    smr_h2prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "SMR_H2Prod", index_col="Bus")
    nh3crack_h2prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "NH3crack_H2Prod", index_col="Bus")
    
    eaf_share = (eaf_prod / (eaf_prod + bof_prod)).fillna(0)
    h2_share = (h2_eaf_prod / (h2_eaf_prod + ng_eaf_prod)).fillna(0)
    green_h2_share = ((elec_h2prod + smrcc_h2prod) / ( elec_h2prod + smrcc_h2prod + smr_h2prod + nh3crack_h2prod)).fillna(0)
    
    if config['sector']['H2_network']:
        # H2 network true -> let's do an average H2 green for Europe
        h2_prod = elec_h2prod + smrcc_h2prod + smr_h2prod + nh3crack_h2prod
        green_h2_prod = green_h2_share * h2_prod
        green_h2_share = round(green_h2_prod.sum() / h2_prod.sum(),1)
    
    green_share = round((eaf_share * h2_share * green_h2_share)*100)
    steel_prod = eaf_prod + bof_prod   
    
    steel_emis = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "Steel_Proc_Emis", index_col="Bus")
    
    return steel_prod, green_share, steel_emis


def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue
            names = ifind.index[ifind == i]
            c.df.loc[names, "location"] = names.str[:i]


def plot_map(n, excel_data, emis_data, regions, year, title,i,j,ncol, ax=None):

    assign_location(n)
    #timestep = n.snapshot_weightings.iloc[0,0]
    
    df = pd.DataFrame(0, index = regions.index, columns = ['data'])
    df["prefix"] = df.index.str.split().str[0].str[:2]
    df["data"] = df["prefix"].map(excel_data[year]) 
    emis_factor = emis_data / (excel_data / 1e3 )
    df["emis"] = df["prefix"].map(emis_factor[year])
  
    
    df.drop(columns=["prefix"], inplace=True)
    
    regions["data"] = df["data"]/1e3 #kt to Mt
    regions["data"] = regions["data"].where(regions["data"] > 0.1, 0)
    regions["data"] = regions["data"].fillna(0)
    regions["emis"] = df["emis"]
    regions["emis"] = regions["emis"].where(regions["emis"] > 0.1, 0)
    regions["emis"] = regions["emis"].fillna(0)
    
    regions = regions.to_crs(proj.proj4_init)
    max_value = excel_data.values.max()/1e3

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": proj})

    regions.plot(
        ax=ax,
        column="data",
        cmap="Blues",
        linewidths=0.5,  # Thickness of the black border
        edgecolor="black",  # Black border for the shapes
        legend=(i == ncol - 1),
        vmax=max_value,
        vmin=0,
        legend_kwds={
            "label": title if i == ncol - 1 else "",
            "shrink": 0.5,
            "extend": "max",
        },
    )
    ax.set_facecolor("white")
    
    # Annotate centroids with emis values
    prev_idx = 0
    max_emis_value = regions["emis"].max()
    for idx, region in regions.iterrows():
        if idx[:2] == prev_idx:
            prev_idx = idx[:2]
            continue
        else:
            centroid = region['geometry'].centroid
            if region["data"] > 1 :#and region["emis"] >= (0.1) * max_emis_value:
                color = "white" if region["data"] > (2/3) * max_value else "black"
                ax.annotate(
                    text=f"{int(region['emis'])}",
                    xy=(centroid.x, centroid.y),  
                    fontsize=8,
                    ha="center",
                    color=color,
                )
            prev_idx = idx[:2]

    #cool_name = cool_names[scenarios.index(scenario)]
    #ax.set_title(cool_name, fontsize=12, loc="center")



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
    
    steel_prod, green_share, steel_emis = steel_preprocessing(excel_dir, scenario, config)

    ax = axes[j, i]
    plot_map(n, steel_prod, steel_emis, regions, year, title, i,j,ncol, ax=ax)
    
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
        
            
fig.text(1, 1, 'Numbers are\nemis. factors [tCO2/ktsteel]', ha='right', va='top', fontsize=10, color='grey')

plt.tight_layout()
fig.suptitle('Year 2050', fontsize=16, y=1.02)
plt.savefig("graphs/steel_eu_dem_2050.png", bbox_inches="tight")
plt.show()

