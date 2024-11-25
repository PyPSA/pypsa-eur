# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:46:12 2024

@author: alice
"""


def load_projection(plotting_params):
    proj_kwargs = plotting_params.get("projection", dict(name="EqualEarth"))
    proj_func = getattr(ccrs, proj_kwargs.pop("name"))
    return proj_func(**proj_kwargs)

def ammonia_preprocessing(excel_dir, scenario, config):
    ammonia_prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "Ammonia_Prod", index_col="Bus")
    
    # Hydrogen
    elec_h2prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "Elec_H2Prod", index_col="Bus")
    smrcc_h2prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "SMRCC_H2Prod", index_col="Bus")
    smr_h2prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "SMR_H2Prod", index_col="Bus")
    nh3crack_h2prod = pd.read_excel(excel_dir + scenario + ".xlsx", sheet_name = "NH3crack_H2Prod", index_col="Bus")
    
    green_h2_share = ((elec_h2prod + smrcc_h2prod) / ( elec_h2prod + smrcc_h2prod + smr_h2prod + nh3crack_h2prod)).fillna(0)
    
    if config['sector']['H2_network']:
        # H2 network true -> let's do an average H2 green for Europe
        h2_prod = elec_h2prod + smrcc_h2prod + smr_h2prod + nh3crack_h2prod
        green_h2_prod = green_h2_share * h2_prod
        green_h2_share = round(green_h2_prod.sum() / h2_prod.sum(),1)
    
    green_share = round((ammonia_prod / ammonia_prod) * green_h2_share*100)
    
    return ammonia_prod, green_share


def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue
            names = ifind.index[ifind == i]
            c.df.loc[names, "location"] = names.str[:i]


def plot_map(n, excel_data, share_data, regions, year, title, ax=None):

    assign_location(n)
    #timestep = n.snapshot_weightings.iloc[0,0]
    
    df = pd.DataFrame(0, index = regions.index, columns = ['data'])
    df["prefix"] = df.index.str.split().str[0].str[:2]
    df["data"] = df["prefix"].map(excel_data[year]) # Assigns same values for one country two nodes for graph
    df["share"] = df["prefix"].map(share_data[year])
  
    
    df.drop(columns=["prefix"], inplace=True)
    
    regions["data"] = df["data"]
    regions["data"] = regions["data"].where(regions["data"] > 0.1, 0)
    regions["data"] = regions["data"].fillna(0)
    regions["share"] = df["share"]
    regions["share"] = regions["share"].where(regions["share"] > 0.1, 0)
    regions["share"] = regions["share"].fillna(0)
    
    regions = regions.to_crs(proj.proj4_init)
    max_value = excel_data.values.max()

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": proj})

    regions.plot(
        ax=ax,
        column="data",
        cmap="Blues",
        linewidths=0.5,  # Thickness of the black border
        edgecolor="black",  # Black border for the shapes
        legend=True,
        vmax=max_value,
        vmin=0,
        legend_kwds={
            "label": title,
            "shrink": 0.5,
            "extend": "max",
        },
    )
    ax.set_facecolor("white")
    
    # Annotate centroids with share values
    prev_idx = 0
    for idx, region in regions.iterrows():
        if idx[:2] == prev_idx:
            prev_idx = idx[:2]
            continue
        else:
            centroid = region['geometry'].centroid
            if region["data"] > 1:
                ax.annotate(
                    text=f"{int(region['share'])}",  # Format share as a percentage
                    xy=(centroid.x, centroid.y),  # Use centroid coordinates
                    fontsize=9,
                    ha="center",
                    color="black",
                )
            prev_idx = idx[:2]
            

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

with open(
    root_dir + res_dir + "baseline_eu_dem/configs/config.base_s_39_lvopt___2030.yaml"
) as config_file:
    config = yaml.safe_load(config_file)

regions = gpd.read_file(regions_fn).set_index("name")
map_opts = config["plotting"]["map"]

if map_opts["boundaries"] is None:
    map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

config["plotting"]["projection"]["name"] = "EqualEarth"
proj = load_projection(config["plotting"])

years = [2030, 2040, 2050]

excel_dir = "excels/"
# Plotting the green share of steel production
title = "Ammonia prod [Mt/yr]"

scenarios = ["baseline_eu_dem", "climate_policy_eu_dem"]

fig, axes = plt.subplots(
    len(scenarios),
    len(years),
    figsize=(3 * len(years), 3 * len(scenarios)),
    subplot_kw={"projection": proj},
)

fn = (root_dir + "results/" + scenarios[0] + "/postnetworks/base_s_39_lvopt___2030.nc")
n = pypsa.Network(fn)

for j, scenario in enumerate(scenarios):
    
    ammonia_prod, green_share = ammonia_preprocessing(excel_dir, scenario, config)

    for i, year in enumerate(years):
        
        #fn = (root_dir + "results/" + scenario + f"/postnetworks/base_s_39_lvopt___{year}.nc")
        ax = axes[j, i]
        #n = pypsa.Network(fn)
        plot_map(n, ammonia_prod, green_share, regions, year, title, ax=ax)

plt.tight_layout()
fig.suptitle('European demand', fontsize=16, y=1.02)
plt.savefig("graphs/ammonia_eu_dem.png", bbox_inches="tight")
plt.show()


scenarios = ["baseline_regional_dem", "climate_policy_regional_dem"]

fig, axes = plt.subplots(
    len(scenarios),
    len(years),
    figsize=(3 * len(years), 3 * len(scenarios)),
    subplot_kw={"projection": proj},
)

fn = (root_dir + "results/" + scenarios[0] + "/postnetworks/base_s_39_lvopt___2030.nc")
n = pypsa.Network(fn)

for j, scenario in enumerate(scenarios):
    
    ammonia_prod, green_share = ammonia_preprocessing(excel_dir, scenario, config)

    for i, year in enumerate(years):
        
        ax = axes[j, i]
        plot_map(n, ammonia_prod, green_share, regions, year, title, ax=ax)

plt.tight_layout()
fig.suptitle('Regional demand', fontsize=16, y=1.02)
plt.savefig("graphs/ammonia_regional_dem.png", bbox_inches="tight")
plt.show()