# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 14:34:32 2025

@author: alice
"""

import pypsa
import pandas as pd
import matplotlib.pyplot as plt

# Define colors
base_colors = {"regain": "#4F5050", "maintain": "#85877C", "deindustrial": "#B0B2A1"}
policy_colors = {"regain": "#5D8850", "maintain": "#95BF74", "deindustrial": "#C5DEB1"}


# Replace ISO2 country codes with full names
country_names = {
    "AL": "Albania", "AT": "Austria", "BA": "Bosnia & Herzegovina", "BE": "Belgium", "BG": "Bulgaria",
    "CH": "Switzerland", "CZ": "Czechia", "DE": "Germany", "DK": "Denmark", "EE": "Estonia", "ES": "Spain",
    "FI": "Finland", "FR": "France", "GB": "UK", "GR": "Greece", "HR": "Croatia", "HU": "Hungary",
    "IE": "Ireland", "IT": "Italy", "LT": "Lithuania", "LU": "Luxembourg", "LV": "Latvia", "ME": "Montenegro",
    "MK": "North Macedonia", "NL": "Netherlands", "NO": "Norway", "PL": "Poland", "PT": "Portugal", "RO": "Romania",
    "RS": "Serbia", "SE": "Sweden", "SI": "Slovenia", "SK": "Slovakia", "XK": "Kosovo"
}

scenarios_steel = ["base_eu_regain", "base_eu_maintain", "base_eu_deindustrial", "policy_eu_regain","policy_eu_maintain", "policy_eu_deindustrial"]
lhv_ammonia = 5.166 # MWh / t
lhv_methanol = 5.528 # MWh / t


# %% FUNCTIONS

scenario = 'base_eu_regain'
n = pypsa.Network(f"C:/Users/alice/Desktop/CMCC/pypsa-adb-industry/results/{scenario}/networks/base_s_39___2050.nc")

# Steel marginal price
steel_price = n.buses_t.marginal_price.loc[:,n.buses_t.marginal_price.columns.str.contains('steel')].mean().iloc[0]/1e3 # €/t steel
mean_hist_steel_price = 415 # €/t steel https://tradingeconomics.com/commodity/steel


# Cement marginal price
cement_price = n.buses_t.marginal_price.loc[
    :, 
    n.buses_t.marginal_price.columns.str.contains('cement') & 
    ~n.buses_t.marginal_price.columns.str.contains('process emissions')
].mean().iloc[0]/1e3 # €/t cement
mean_hist_cement_price = 93 # €/t steel https://www.cemnet.com/News/story/175146/cement-prices-in-italy-see-27-rise-in-april.html

# Ammonia marginal price
ammonia_price = n.buses_t.marginal_price.loc[:,n.buses_t.marginal_price.columns.str.contains('NH3')].mean().iloc[0] * lhv_ammonia
mean_hist_ammonia_price = 470 # €/t https://pubs.usgs.gov/periodicals/mcs2024/mcs2024-nitrogen.pdf

# Methanol marginal price
methanol_price = n.buses_t.marginal_price.loc[
    :, 
    n.buses_t.marginal_price.columns.str.contains('methanol') & 
    ~n.buses_t.marginal_price.columns.str.contains('industry') & 
    ~n.buses_t.marginal_price.columns.str.contains('shipping') 
].mean().iloc[0] * lhv_methanol # €/t methanol

mean_hist_methanol_price = 326 # €/t https://tradingeconomics.com/commodity/methanol

# HVC marginal price
hvc_prices = n.buses_t.marginal_price.loc[:,n.buses_t.marginal_price.columns.str.contains('HVC')].mean().iloc[0]/1e3 # €/t



# %%

years = [2030, 2040, 2050]  # Years in the dataset
lhv_ammonia = 5.166  # MWh / t
lhv_methanol = 5.528  # MWh / t

# Invented values for 2020 (for visualization purposes)
hist_2020_prices = {
    "steel": 415,     # €/t
    "cement": 93,     # €/t
    "ammonia": 470,   # €/t
    "methanol": 326,  # €/t
    "HVC": 600        # €/t
}

# Initialize a dictionary to store DataFrames for each commodity
price_data = {commodity: pd.DataFrame(index=scenarios_steel, columns=[2020] + years) for commodity in hist_2020_prices.keys()}

# Fill 2020 values
for commodity, value in hist_2020_prices.items():
    price_data[commodity][2020] = value

max_value = 0
# Iterate through scenarios and years to calculate prices
for scenario in scenarios_steel:
    for year in years:
        # Load network for the given year
        n = pypsa.Network(f"C:/Users/alice/Desktop/CMCC/pypsa-adb-industry/results/{scenario}/networks/base_s_39___{year}.nc")

        # Steel price
        steel_price = n.buses_t.marginal_price.loc[:, n.buses_t.marginal_price.columns.str.contains('steel')].mean().iloc[0] / 1e3
        price_data["steel"].loc[scenario, year] = steel_price

        # Cement price
        cement_price = n.buses_t.marginal_price.loc[
            :, 
            n.buses_t.marginal_price.columns.str.contains('cement') & 
            ~n.buses_t.marginal_price.columns.str.contains('process emissions')
        ].mean().iloc[0] / 1e3
        price_data["cement"].loc[scenario, year] = cement_price

        # Ammonia price
        ammonia_price = n.buses_t.marginal_price.loc[:, n.buses_t.marginal_price.columns.str.contains('NH3')].mean().iloc[0] * lhv_ammonia
        price_data["ammonia"].loc[scenario, year] = ammonia_price

        # Methanol price
        methanol_price = n.buses_t.marginal_price.loc[
            :, 
            n.buses_t.marginal_price.columns.str.contains('methanol') & 
            ~n.buses_t.marginal_price.columns.str.contains('industry') & 
            ~n.buses_t.marginal_price.columns.str.contains('shipping')
        ].mean().iloc[0] * lhv_methanol
        price_data["methanol"].loc[scenario, year] = methanol_price

        # HVC price
        hvc_price = n.buses_t.marginal_price.loc[:, n.buses_t.marginal_price.columns.str.contains('HVC')].mean().iloc[0] / 1e3
        price_data["HVC"].loc[scenario, year] = hvc_price
        max_value = max(max(df.max().max() for df in price_data.values()), max_value)
        


fig, axes = plt.subplots(1, 6, figsize=(15, 5), sharex=True, sharey=True)

commodities = ["steel", "cement", "ammonia", "methanol", "HVC"]
colors = ["#4F5050", "#85877C", "#B0B2A1", "#5D8850", "#95BF74", "#C5DEB1"]

for idx, (commodity, ax) in enumerate(zip(commodities, axes)):
    for i, scenario in enumerate(scenarios_steel):
        ax.plot([2020] + years, price_data[commodity].loc[scenario], marker="o", linestyle="-", color=colors[i], label=scenario if idx == 0 else "")
    
    ax.set_title(f"{commodity.capitalize()} Price")
    ax.set_xticks([2020] + years)
    ax.set_ylim(0, max_value)  
    if idx == 0:
        ax.set_ylabel("Price (€/t)")
    ax.grid(True)

# Extra invisible subplot for the legend
ax_legend = axes[-1]
ax_legend.axis("off")

# Creating legend
handles = []
labels = []
for scenario in scenarios_steel:
    strategy = next(key for key in base_colors.keys() if key in scenario)
    color = policy_colors[strategy] if "policy" in scenario else base_colors[strategy]
    handles.append(plt.Line2D([0], [0], color=color, marker="o", linestyle="-"))
    labels.append(scenario)

ax_legend.legend(handles, labels, loc="center", frameon=False)


plt.tight_layout()
plt.savefig("./graphs/commodity_prices.png")
plt.show()
