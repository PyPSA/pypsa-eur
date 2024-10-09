# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:46:12 2024

@author: alice
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

def plot_steel_map(n, regions, foresight, ax=None): 
    #n = pypsa.Network(root_dir + res_dir + scenario + "postnetworks/base_s_39_lvopt___2030.nc")

    
    # if "H2 pipeline" not in n.links.carrier.unique():
    #     return
    
    assign_location(n)
    timestep = n.snapshot_weightings.iloc[0,0]
    
    steel_prod_index = n.links.query("bus1 == 'EU steel'").index
    steel_prod = -n.links_t.p1.loc[:,steel_prod_index].sum()#*timestep/1e3 #Mtsteel
    steel_prod.index = steel_prod.index.str.split(' 0 ').str[0] + ' 0'
    steel_prod_df = steel_prod.to_frame(name='steel_prod')
    
    regions["steel"] = (
        steel_prod_df
        .steel_prod.groupby(level=0)
        .sum()
        .div(1e3)*timestep
    )  # TWh
    
    
    regions["steel"] = regions["steel"].where(regions["steel"] > 0.1, 0)
    regions["steel"] = regions["steel"].fillna(0)
    
    
    # make a fake MultiIndex so that area is correct for legend
    #bus_sizes.rename(index=lambda x: x.replace(" H2", ""), level=0, inplace=True)
    # drop all links which are not H2 pipelines
    n.links.drop(n.links.index[~n.links.index.str.contains("EAF|BOF", case=False)], inplace=True)
    
    regions = regions.to_crs(proj.proj4_init)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": proj})
    
    
    
    regions.plot(
        ax=ax,
        column="steel",
        cmap="Blues",
        linewidths=0,
        legend=True,
        vmax=6,
        vmin=0,
        legend_kwds={
            "label": "Steel production [Mt steel/yr]",
            "shrink": 0.7,
            "extend": "max",
        },
    )
    
    
    ax.set_facecolor("white")

    
#%%


root_dir = "C:/Users/alice/Desktop/CMCC/pypsa-adb-industry/"
scenario = "baseline/"
res_dir = "results/"
regions_fn = root_dir + "resources/" + scenario + "regions_onshore_base_s_39.geojson"

import pypsa
import pandas as pd
import geopandas as gpd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import yaml
with open(root_dir + res_dir + scenario + "configs/config.base_s_39_lvopt___2030.yaml") as config_file: config = yaml.safe_load(config_file)
    
regions = gpd.read_file(regions_fn).set_index("name")
map_opts = config["plotting"]["map"]

if map_opts["boundaries"] is None:
    map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

config["plotting"]["projection"]["name"] = "EqualEarth"
proj = load_projection(config["plotting"])

fig, axes = plt.subplots(1, 3, figsize=(20, 10), subplot_kw={"projection": proj})

for i, year in enumerate([2030, 2040, 2050]):
    for j, scenario in enumerate(["baseline"]): #"policy"]):
        fn = root_dir + "results/" + scenario + f"/postnetworks/base_s_39_lvopt___{year}.nc"
        ax = axes[i]
        n = pypsa.Network(fn)
        plot_steel_map(n, regions, foresight="myopic", ax=ax)

plt.tight_layout()

plt.savefig("graphs/steel_prod_per_country.png", bbox_inches='tight')
plt.show()
