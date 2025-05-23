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


def retrieve_loads(df, suffix):

    filtered_df = df[[col for col in df.columns if col.endswith(suffix)]]
    filtered_df.columns = filtered_df.columns.str.split(" 0").str[0] + " 0"
    return filtered_df


def retrieve_links(df, suffix):

    filtered_df = df.loc[
        :, n.links[(n.links.index.str.contains(suffix, case=False))].index
    ]  # .sum()
    filtered_df.columns = filtered_df.columns.str.split(" 0").str[0] + " 0"
    filtered_df = filtered_df.T.groupby(filtered_df.columns).sum().T
    filtered_df = filtered_df[
        [col for col in filtered_df.columns if not col.startswith("EU")]
    ]
    return filtered_df


def plot_elec_prices_map(n, regions, year, ax=None):

    assign_location(n)

    mprice = n.buses_t.marginal_price
    elec_prices = mprice[[col for col in mprice.columns if col.endswith(" 0")]]

    # Get electricity loads per hours

    # DF from the network
    loads = n.loads_t.p
    links0 = n.links_t.p0
    links1 = n.links_t.p1
    links2 = n.links_t.p2

    elec_loads_resi = retrieve_loads(loads, "0")
    elec_loads_ind = retrieve_loads(
        loads, "industry electricity"
    )  # This will disappear
    elec_loads_agri = retrieve_loads(loads, "agriculture electricity")
    elec_loads_tra = retrieve_loads(loads, "EV")
    elec_links_th = retrieve_links(links0, "thermal|heat")
    elec_links_dac = retrieve_links(links1, "DAC")
    elec_links_meth = retrieve_links(links2, "methanolisation")

    elec_loads = (
        elec_loads_resi
        + elec_loads_ind
        + elec_loads_agri
        + elec_loads_tra
        + elec_links_th
        + elec_links_dac
        + elec_links_meth
    )

    # Get the average of prices per node based on power load
    total_elec_expen = elec_prices * elec_loads
    weighted_average = total_elec_expen.sum() / elec_loads.sum()  # €/MWh

    regions["elec_price"] = weighted_average  # €/MWh

    # drop all links which are not H2 pipelines
    # n.links.drop(n.links.index[~n.links.index.str.contains("EAF|BOF", case=False)], inplace=True)

    regions = regions.to_crs(proj.proj4_init)

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": proj})

    regions.plot(
        ax=ax,
        column="elec_price",
        cmap="Reds",
        linewidths=0,
        legend=True,
        vmax=60,
        vmin=10,
        legend_kwds={
            "label": "Electricity price [€/MWh]",
            # "fontsize": 12,
            "shrink": 0.5,
            "extend": "max",
        },
    )

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

with open(
    root_dir + res_dir + scenario + "configs/config.base_s_39_lvopt___2030.yaml"
) as config_file:
    config = yaml.safe_load(config_file)

regions = gpd.read_file(regions_fn).set_index("name")
map_opts = config["plotting"]["map"]

if map_opts["boundaries"] is None:
    map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

config["plotting"]["projection"]["name"] = "EqualEarth"
proj = load_projection(config["plotting"])

years = [2030, 2040, 2050]
scenarios = ["baseline", "climate_policy"]

fig, axes = plt.subplots(
    len(scenarios),
    len(years),
    figsize=(3 * len(years), 3 * len(scenarios)),
    subplot_kw={"projection": proj},
)


for i, year in enumerate(years):
    for j, scenario in enumerate(scenarios):
        fn = (
            root_dir
            + "results/"
            + scenario
            + f"/postnetworks/base_s_39_lvopt___{year}.nc"
        )
        ax = axes[j, i]
        n = pypsa.Network(fn)
        plot_elec_prices_map(n, regions, year, ax=ax)

plt.tight_layout()

plt.savefig("graphs/elec_price_per_country.png", bbox_inches="tight")
plt.show()
