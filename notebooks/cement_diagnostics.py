# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 16:36:18 2024

@author: alice
"""

import matplotlib.pyplot as plt
import pandas as pd
import pypsa

res_directory = "../results/climate_policy_eu_dem/postnetworks/"


# %% Script to get production in different countries and demand for a single year

n = pypsa.Network(res_directory + "base_s_39_lvopt___2030.nc")
alinks = n.links
aloads = n.loads
astores = n.stores
aconstr = n.global_constraints
# cement demand
timestep = n.snapshot_weightings.iloc[0,0]
cement_dem = n.loads_t.p.filter(like='cement', axis=1).sum().sum()*timestep #kt cement

# cement production
cement_p0 = n.links_t.p0.filter(like='Cement', axis=1)
cement_p1 = -n.links_t.p1.filter(like='Cement', axis=1)
cement_p2 = n.links_t.p2.filter(like='Cement', axis=1)
cement_p3 = n.links_t.p3.filter(like='Cement', axis=1)
cement_emissions_p0 = n.links_t.p0.filter(like='cement to total proc emis', axis=1)
cement_emis_p0 = n.links_t.p0.filter(like='cement to total proc emis', axis=1)
cement_emis_p1 = -n.links_t.p1.filter(like='cement to total proc emis', axis=1)
cement_emis_p2 = -n.links_t.p2.filter(like='cement to total proc emis', axis=1)
dac_p0 = n.links_t.p0.filter(like='DAC', axis=1)
cement_prod = cement_p1.loc[:, (cement_p1 > 1e-10).any(axis=0)]
cement_prod_per_country = cement_prod.sum() * timestep / 1e3  # Mt
cement_prod_per_country.index = cement_prod_per_country.index.str[:2]
cement_prod_per_country = cement_prod_per_country.groupby(cement_prod_per_country.index).sum()
cement_prod_per_country = pd.DataFrame(cement_prod_per_country, columns=["Values"])

cement_cap = n.links[n.links.index.str.contains('Cement', regex=True)]


# Create a figure and axis
fig, ax = plt.subplots()

# Plot the stacked bar chart for cement production
cement_prod_per_country.T.plot.bar(stacked=True, ax=ax, title="Cement production")


# Overlay a dot on top of the total cement production
ax.plot("Total", cement_dem / 1e3, "ro", label="Cement Dem")

# Adjust the legend to be in 3 columns
ax.legend(title="", ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.15))
# Add labels and title
ax.set_ylabel("Cement Production (Mt/yr)")

# Adjust the layout to fit labels
plt.tight_layout()

# Show the plot
plt.show()

# %% Script to get production in different countries and demand for 3 investment years

years = [2030, 2040, 2050]
n_dict = {
    year: pypsa.Network(res_directory + "base_s_39_lvopt___" + str(year) + ".nc")
    for year in years
}

# Initialize empty dictionaries to store cement demand and cement production per country for each year
cement_dem_dict = {}
cement_prod_dict = {}
cement_tech_dict = {}

for year, n in n_dict.items():
    # cement demand calculation
    timestep = n_dict[year].snapshot_weightings.iloc[0, 0]
    cement_dem = n_dict[year].loads_t.p.filter(like='cement', axis=1).sum().sum() * timestep  # Mt cement
    cement_dem_dict[year] = cement_dem/1e3  # Store cement demand for this year

    # cement production calculation
    cement_prod = -n_dict[year].links_t.p1.filter(regex="Cement", axis=1)
    cement_prod = cement_prod.loc[:, (cement_prod > 1e-6).any(axis=0)]
    cement_prod_per_country = cement_prod.sum() * timestep / 1e3  # Convert to Gt cement
    cement_prod_per_tech = cement_prod_per_country.copy()
    cement_prod_per_country.index = cement_prod_per_country.index.str[
        :2
    ]  # Keep only the first two characters (country code)
    cement_prod_per_country = cement_prod_per_country.groupby(
        cement_prod_per_country.index
    ).sum()  # Group by country code
    cement_prod_dict[year] = (
        cement_prod_per_country  # Store cement production for this year
    )

    cement_prod_per_tech.index = cement_prod_per_tech.index.str.extract(
        "(Cement)", expand=False
    )
    cement_prod_per_tech = cement_prod_per_tech.groupby(
        cement_prod_per_tech.index
    ).sum()  # Group by country code
    cement_tech_dict[year] = cement_prod_per_tech  # Store cement production for this year


#cement_dem_df = pd.DataFrame(cement_dem_dict)
cement_dem_df = pd.DataFrame.from_dict(cement_dem_dict, orient='index', columns=[0]).T

cement_prod_df = pd.DataFrame(cement_prod_dict)
cement_prod_df = cement_prod_df.fillna(0)

cement_tech_df = pd.DataFrame(cement_tech_dict)
cement_tech_df = cement_tech_df.fillna(0)

fig, ax = plt.subplots(figsize=(10, 6))

cement_prod_df.T.plot.bar(stacked=True, ax=ax, width=0.8, color=plt.cm.tab20.colors)
bar_positions = range(len(cement_prod_df.sum()))
ax.scatter(
    bar_positions, cement_dem_df.iloc[0], color="black", zorder=3, label="Cement Demand"
)

ax.set_ylabel("Cement Production (Gt)")
ax.legend(title="Country", loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
plt.xticks(rotation=45)
ax.grid(True, axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()

# %% Plot of BOF vs EAF

fig, ax = plt.subplots(figsize=(10, 6))

cement_tech_df.T.plot.bar(stacked=True, ax=ax, width=0.8, color=plt.cm.tab20.colors)
bar_positions = range(len(cement_tech_df.sum()))
ax.scatter(
    bar_positions, cement_dem_df.iloc[0], color="black", zorder=3, label="cement Demand"
)

ax.set_ylabel("cement Production (Gt)")
ax.legend(title="Country", loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
plt.xticks(rotation=45)
ax.grid(True, axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()

# %% Get the capacities

# Initialize empty dictionaries to store cement demand and cement production per country for each year
cement_cap_country_dict = {}
cement_cap_tech_dict = {}

for year, n in n_dict.items():
    # cement capacities calculation
    timestep = n_dict[year].snapshot_weightings.iloc[0, 0]
    cement_cap = n_dict[year].links.p_nom.filter(regex="EAF|BOF", axis=0)
    cement_cap = cement_cap[cement_cap > 1e-6]
    cement_cap_per_country = (
        cement_cap * timestep / 1e3
    )  # Convert from kt/h to Mt/yr of cement output
    cement_cap_per_tech = cement_cap_per_country.copy()
    cement_cap_per_country.index = cement_cap_per_country.index.str[
        :2
    ]  # Keep only the first two characters (country code)
    cement_cap_per_country = cement_cap_per_country.groupby(
        cement_cap_per_country.index
    ).sum()  # Group by country code
    cement_cap_country_dict[year] = (
        cement_cap_per_country  # Store cement production for this year
    )

    cement_cap_per_tech.index = cement_cap_per_tech.index.str.extract(
        "(EAF|BOF)", expand=False
    )
    cement_cap_per_tech = cement_cap_per_tech.groupby(
        cement_cap_per_tech.index
    ).sum()  # Group by country code
    cement_cap_tech_dict[year] = cement_cap_per_tech


cement_cap_country_df = pd.DataFrame(cement_cap_country_dict)
cement_cap_country_df = cement_cap_country_df.fillna(0)

cement_cap_tech_df = pd.DataFrame(cement_cap_tech_dict)
cement_cap_tech_df = cement_cap_tech_df.fillna(0)

fig, ax = plt.subplots(figsize=(10, 6))

cement_cap_country_df.T.plot.bar(
    stacked=True, ax=ax, width=0.8, color=plt.cm.tab20.colors
)
bar_positions = range(len(cement_cap_country_df.sum()))
# ax.scatter(bar_positions, cement_cap_country_df.iloc[0], color='black', zorder=3, label='cement Demand')

ax.set_ylabel("cement capacities (Mt/yr cement)")
ax.legend(title="Country", loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
plt.xticks(rotation=45)
ax.grid(True, axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()


#%% Get the capacities with hydrogen

# Initialize empty dictionaries to store cement demand and cement production per country for each year
cement_cap_country_dict = {}
cement_cap_tech_dict = {}

cement_dri = n.links.p_nom_opt.filter(regex='DRI', axis = 0)
dri0 = n.links_t.p0[cement_dri.index]
dri1 = n.links_t.p1[cement_dri.index]

filtered_df = n.links[n.links['bus1'].str.contains('H2', case=False, na=False)]



for year, n in n_dict.items():
    # cement capacities calculation
    timestep = n_dict[year].snapshot_weightings.iloc[0, 0]
    cement_cap = n_dict[year].links.p_nom.filter(regex='EAF|BOF', axis=0)
    cement_cap = cement_cap[cement_cap > 1e-6]
    cement_cap_per_country = cement_cap * timestep / 1e3  # Convert from kt/h to Mt/yr of cement output
    cement_cap_per_tech = cement_cap_per_country.copy()
    cement_cap_per_country.index = cement_cap_per_country.index.str[:2]  # Keep only the first two characters (country code)
    cement_cap_per_country = cement_cap_per_country.groupby(cement_cap_per_country.index).sum()  # Group by country code
    cement_cap_country_dict[year] = cement_cap_per_country  # Store cement production for this year
    
    cement_cap_per_tech.index = cement_cap_per_tech.index.str.extract('(EAF|BOF)', expand=False)
    cement_cap_per_tech = cement_cap_per_tech.groupby(cement_cap_per_tech.index).sum()  # Group by country code
    cement_cap_tech_dict[year] = cement_cap_per_tech


cement_cap_country_df = pd.DataFrame(cement_cap_country_dict)
cement_cap_country_df = cement_cap_country_df.fillna(0)

cement_cap_tech_df = pd.DataFrame(cement_cap_tech_dict)
cement_cap_tech_df = cement_cap_tech_df.fillna(0)

fig, ax = plt.subplots(figsize=(10, 6))

cement_cap_country_df.T.plot.bar(stacked=True, ax=ax, width=0.8, color=plt.cm.tab20.colors)
bar_positions = range(len(cement_cap_country_df.sum()))
 #ax.scatter(bar_positions, cement_cap_country_df.iloc[0], color='black', zorder=3, label='cement Demand')

ax.set_ylabel('cement capacities (Mt/yr cement)')
ax.legend(title='Country', loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
plt.xticks(rotation=45)
ax.grid(True, axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()