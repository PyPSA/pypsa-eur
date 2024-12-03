# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 14:55:16 2024

@author: alice
"""

import cartopy.crs as ccrs
import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import yaml

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

def calculate_steel(scenario):
    # Load data from Excel sheets
    excel_dir = "excels_24h/"
    root_dir = "C:/Users/alice/Desktop/CMCC/pypsa-adb-industry/"
    eaf_prod = pd.read_excel(root_dir + "notebooks_december/" + excel_dir + f"{scenario}.xlsx", sheet_name="EAF_Prod", index_col="Bus")
    bof_prod = pd.read_excel(root_dir + "notebooks_december/" +  excel_dir + f"{scenario}.xlsx", sheet_name="BOF_Prod", index_col="Bus")
    
    # Calculate green share and steel production
    steel_prod = eaf_prod + bof_prod
    steel_prod = steel_prod.sum(axis=1) / 1e3  # Convert to Mton for 3 years

    return steel_prod


# %%

root_dir = "C:/Users/alice/Desktop/CMCC/pypsa-adb-industry/"
scenarios = ["baseline_eu_dem", "policy_eu_dem", "baseline_regional_dem", "policy_regional_dem"]
cool_names = ["Baseline European", "Policy European", "Baseline Regional", "Policy Regional"]
years = [2030, 2040, 2050]
columns_to_drop = ["cluster", "Unnamed: 1", "Unnamed: 3"]

base_eu_nodal = pd.read_csv(root_dir + "results_24h/baseline_eu_dem/csvs/nodal_costs.csv")
base_reg_nodal = pd.read_csv(root_dir + "results_24h/baseline_regional_dem/csvs/nodal_costs.csv")
pol_eu_nodal = pd.read_csv(root_dir + "results_24h/policy_eu_dem/csvs/nodal_costs.csv")
pol_reg_nodal = pd.read_csv(root_dir + "results_24h/policy_regional_dem/csvs/nodal_costs.csv")

csvs = [base_eu_nodal, base_reg_nodal, pol_eu_nodal, pol_reg_nodal]
csv_names = ["base_eu", "base_reg", "pol_eu", "pol_reg"]
processed_csvs = []

for df, name in zip(csvs, csv_names):
        
    df = df.drop(columns=columns_to_drop, errors='ignore')  # 'errors=ignore' avoids issues if some columns don't exist
    df = df.dropna(subset=["Unnamed: 2"])
    df = df.set_index("Unnamed: 2")
    df = df[df.index.str.match(r"^[A-Z]")]
    df = df[~df.index.str.startswith("EU")]
    df = df[~df.index.str.startswith("H2")]
    df.columns = years
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.groupby(df.index).sum(numeric_only=True)
    df["scenario"] = name
    df = df.set_index("scenario", append=True)
    processed_csvs.append(df) 

final_combined_df = pd.concat(processed_csvs)
r = 0.04 # social discount rate

for year in years:
    final_combined_df[f'Annualized_{year}'] = final_combined_df[year] / ((1 + r) ** (year - 2020))

annualized_df = final_combined_df['Annualized_2030']  + final_combined_df['Annualized_2040'] + final_combined_df['Annualized_2050']

# %% Map of difference between pol_eu and pol_reg


annualized_df = annualized_df.reset_index()

annualized_df.columns = ['Country','Scenario','Cumulative_costs']
annualized_df.index = annualized_df['Country']
annualized_df = annualized_df.drop('Country', axis=1)
# Weight by population
pop = pd.read_csv(root_dir + "resources/baseline/pop_layout_base_s_39.csv", index_col = 0)
pop = pop['total']
annualized_df = annualized_df.merge(pop, left_index=True, right_index=True, how='left', suffixes=('', '_pop'))
annualized_df['Costs_per_person'] = annualized_df['Cumulative_costs'] / annualized_df['total'] / 1e3 # €/person

#annualized_df.index = annualized_df.index.str[:2]

annualized_df = annualized_df.groupby(["Country", "Scenario"]).sum()

annualized_df = annualized_df.reset_index()
annualized_df.index = annualized_df['Country']
annualized_df = annualized_df.drop(['Country','Cumulative_costs','total'], axis = 1)

annualized_reg = annualized_df[annualized_df['Scenario'] == 'pol_reg']
annualized_eu = annualized_df[annualized_df['Scenario'] == 'pol_eu'].drop('Scenario', axis = 1)
annualized_change = annualized_reg.merge(annualized_eu, left_index=True, right_index=True, how='left', suffixes=('', '_REF'))
annualized_change['Change'] = (annualized_change['Costs_per_person'] - annualized_change['Costs_per_person_REF'])

annualized_change = annualized_change['Change']
annualized_eu = annualized_eu['Costs_per_person']

# %% Get steel production for the graph

# Call the function for different scenarios
policy_reg_steel_prod = calculate_steel("policy_regional_dem")
policy_eu_steel_prod = calculate_steel("policy_eu_dem") 

delta_steel = policy_reg_steel_prod - policy_eu_steel_prod

# %% Plot

fn = (root_dir + "results/" + scenarios[0] + "/postnetworks/base_s_39_lvopt___2050.nc")
n = pypsa.Network(fn)
res_dir = "results_24h/"

with open(
    root_dir + res_dir + "baseline_eu_dem/configs/config.base_s_39_lvopt___2050.yaml"
) as config_file:
    config = yaml.safe_load(config_file)



config["plotting"]["projection"]["name"] = "EqualEarth"
proj = load_projection(config["plotting"])


fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": proj})

assign_location(n)
regions_fn = root_dir + "resources/baseline/regions_onshore_base_s_39.geojson"

regions = gpd.read_file(regions_fn).set_index("name")
map_opts = config["plotting"]["map"]

if map_opts["boundaries"] is None:
    map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

df = pd.DataFrame(0, index = regions.index, columns = ['change','ref'])
df["change"] = df.index.map(annualized_change) 
df["ref"] = df.index.map(annualized_eu)
df["index_prefix"] = df.index.str[:2]
average_changes = df.groupby("index_prefix")["change"].transform("mean")
df["change"] = average_changes
df['steel'] = df['index_prefix'].map(delta_steel)
df = df.drop(columns=["index_prefix"])

regions["change"] = df["change"]
regions["steel"] = df["steel"]


regions = regions.to_crs(proj.proj4_init)
max_value = df["change"].values.max()
min_value = df["change"].values.min()
if ax is None:
    fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": proj})

regions.plot(
    ax=ax,
    column="change",
    cmap="seismic",
    linewidths=0.5,  # Thickness of the black border
    edgecolor="black",  # Black border for the shapes
    legend=True,
    vmax=max_value,
    vmin=min_value,
    legend_kwds={
        "label": "Delta cumulative costs per person\nwrt European scenario [€/person]",
        "shrink": 0.5,
        "extend": "max",
    },
)
ax.set_facecolor("white")


# Annotate centroids with emis values
prev_idx = 0
max_emis_value = regions["steel"].max()
for idx, region in regions.iterrows():
    if idx[:2] == prev_idx:
        prev_idx = idx[:2]
        continue
    else:
        centroid = region['geometry'].centroid
        ax.annotate(
            text=f"{int(region['steel'])}",
            xy=(centroid.x, centroid.y),  
            fontsize=10,
            ha="center",
            color='black',
        )
        prev_idx = idx[:2]


fig.text(0.95, 0.95, 'Numbers are delta cumulative\nMton of steel produced nationally', ha='right', va='top', fontsize=10, color='grey')

plt.tight_layout()
fig.suptitle('', fontsize=16, y=1.02)
plt.savefig("graphs/costs_changes.png", bbox_inches="tight")
plt.show()


#Plotting in color the (cumulative costs / person)regional - (cumulative costs / person) europe and in the label the (steel production)regional - (steel_production)europe)
