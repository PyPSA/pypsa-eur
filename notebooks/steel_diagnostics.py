# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 16:36:18 2024

@author: alice
"""

import matplotlib.pyplot as plt
import pandas as pd
import pypsa

res_directory = "../results/baseline/postnetworks/"


# %% Script to get production in different countries and demand for a single year

n = pypsa.Network(res_directory + "base_s_39_lvopt___2030.nc")

# Steel demand
timestep = n.snapshot_weightings.iloc[0, 0]
steel_dem = n.loads_t.p.filter(like="steel", axis=1).sum() * timestep  # Mt steel

# Steel production
steel_prod = -n.links_t.p1.filter(regex="EAF|BOF", axis=1)
steel_prod = steel_prod.loc[:, (steel_prod > 1e-10).any(axis=0)]
steel_prod_per_country = steel_prod.sum() * timestep / 1e3  # Gt
steel_prod_per_country.index = steel_prod_per_country.index.str[:2]
steel_prod_per_country = steel_prod_per_country.groupby(
    steel_prod_per_country.index
).sum()
steel_prod_per_country = pd.DataFrame(steel_prod_per_country, columns=["Values"])

# Create a figure and axis
fig, ax = plt.subplots()

# Plot the stacked bar chart for steel production
steel_prod_per_country.T.plot.bar(stacked=True, ax=ax, title="Steel production")


# Overlay a dot on top of the total steel production
ax.plot("Total", steel_dem / 1e3, "ro", label="Steel Dem")

# Adjust the legend to be in 3 columns
ax.legend(title="", ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.15))
# Add labels and title
ax.set_ylabel("Steel Production (Gt/yr)")

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

# Initialize empty dictionaries to store steel demand and steel production per country for each year
steel_dem_dict = {}
steel_prod_dict = {}
steel_tech_dict = {}

for year, n in n_dict.items():
    # Steel demand calculation
    timestep = n_dict[year].snapshot_weightings.iloc[0, 0]
    steel_dem = (
        n_dict[year].loads_t.p.filter(like="steel", axis=1).sum() * timestep
    )  # Mt steel
    steel_dem_dict[year] = steel_dem / 1e3  # Store steel demand for this year

    # Steel production calculation
    steel_prod = -n_dict[year].links_t.p1.filter(regex="EAF|BOF", axis=1)
    steel_prod = steel_prod.loc[:, (steel_prod > 1e-6).any(axis=0)]
    steel_prod_per_country = steel_prod.sum() * timestep / 1e3  # Convert to Gt steel
    steel_prod_per_tech = steel_prod_per_country.copy()
    steel_prod_per_country.index = steel_prod_per_country.index.str[
        :2
    ]  # Keep only the first two characters (country code)
    steel_prod_per_country = steel_prod_per_country.groupby(
        steel_prod_per_country.index
    ).sum()  # Group by country code
    steel_prod_dict[year] = (
        steel_prod_per_country  # Store steel production for this year
    )

    steel_prod_per_tech.index = steel_prod_per_tech.index.str.extract(
        "(EAF|BOF)", expand=False
    )
    steel_prod_per_tech = steel_prod_per_tech.groupby(
        steel_prod_per_tech.index
    ).sum()  # Group by country code
    steel_tech_dict[year] = steel_prod_per_tech  # Store steel production for this year


steel_dem_df = pd.DataFrame(steel_dem_dict)

steel_prod_df = pd.DataFrame(steel_prod_dict)
steel_prod_df = steel_prod_df.fillna(0)

steel_tech_df = pd.DataFrame(steel_tech_dict)
steel_tech_df = steel_tech_df.fillna(0)

fig, ax = plt.subplots(figsize=(10, 6))

steel_prod_df.T.plot.bar(stacked=True, ax=ax, width=0.8, color=plt.cm.tab20.colors)
bar_positions = range(len(steel_prod_df.sum()))
ax.scatter(
    bar_positions, steel_dem_df.iloc[0], color="black", zorder=3, label="Steel Demand"
)

ax.set_ylabel("Steel Production (Gt)")
ax.legend(title="Country", loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
plt.xticks(rotation=45)
ax.grid(True, axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()

# %% Plot of BOF vs EAF

fig, ax = plt.subplots(figsize=(10, 6))

steel_tech_df.T.plot.bar(stacked=True, ax=ax, width=0.8, color=plt.cm.tab20.colors)
bar_positions = range(len(steel_tech_df.sum()))
ax.scatter(
    bar_positions, steel_dem_df.iloc[0], color="black", zorder=3, label="Steel Demand"
)

ax.set_ylabel("Steel Production (Gt)")
ax.legend(title="Country", loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
plt.xticks(rotation=45)
ax.grid(True, axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()

# %% Get the capacities

# Initialize empty dictionaries to store steel demand and steel production per country for each year
steel_cap_country_dict = {}
steel_cap_tech_dict = {}

for year, n in n_dict.items():
    # Steel capacities calculation
    timestep = n_dict[year].snapshot_weightings.iloc[0, 0]
    steel_cap = n_dict[year].links.p_nom.filter(regex="EAF|BOF", axis=0)
    steel_cap = steel_cap[steel_cap > 1e-6]
    steel_cap_per_country = (
        steel_cap * timestep / 1e3
    )  # Convert from kt/h to Mt/yr of steel output
    steel_cap_per_tech = steel_cap_per_country.copy()
    steel_cap_per_country.index = steel_cap_per_country.index.str[
        :2
    ]  # Keep only the first two characters (country code)
    steel_cap_per_country = steel_cap_per_country.groupby(
        steel_cap_per_country.index
    ).sum()  # Group by country code
    steel_cap_country_dict[year] = (
        steel_cap_per_country  # Store steel production for this year
    )

    steel_cap_per_tech.index = steel_cap_per_tech.index.str.extract(
        "(EAF|BOF)", expand=False
    )
    steel_cap_per_tech = steel_cap_per_tech.groupby(
        steel_cap_per_tech.index
    ).sum()  # Group by country code
    steel_cap_tech_dict[year] = steel_cap_per_tech


steel_cap_country_df = pd.DataFrame(steel_cap_country_dict)
steel_cap_country_df = steel_cap_country_df.fillna(0)

steel_cap_tech_df = pd.DataFrame(steel_cap_tech_dict)
steel_cap_tech_df = steel_cap_tech_df.fillna(0)

fig, ax = plt.subplots(figsize=(10, 6))

steel_cap_country_df.T.plot.bar(
    stacked=True, ax=ax, width=0.8, color=plt.cm.tab20.colors
)
bar_positions = range(len(steel_cap_country_df.sum()))
# ax.scatter(bar_positions, steel_cap_country_df.iloc[0], color='black', zorder=3, label='Steel Demand')

ax.set_ylabel("Steel capacities (Mt/yr steel)")
ax.legend(title="Country", loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
plt.xticks(rotation=45)
ax.grid(True, axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()
