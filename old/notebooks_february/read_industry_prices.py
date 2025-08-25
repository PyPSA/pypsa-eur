"""
Created on Thu Mar 13 14:34:32 2025

@author: alice
"""

import os

import matplotlib.pyplot as plt
import pandas as pd
import pypsa

# Define colors
base_colors = {"regain": "#4F5050", "maintain": "#85877C", "deindustrial": "#B0B2A1"}
policy_colors = {"regain": "#5D8850", "maintain": "#95BF74", "deindustrial": "#C5DEB1"}


# Replace ISO2 country codes with full names
country_names = {
    "AL": "Albania",
    "AT": "Austria",
    "BA": "Bosnia & Herzegovina",
    "BE": "Belgium",
    "BG": "Bulgaria",
    "CH": "Switzerland",
    "CZ": "Czechia",
    "DE": "Germany",
    "DK": "Denmark",
    "EE": "Estonia",
    "ES": "Spain",
    "FI": "Finland",
    "FR": "France",
    "GB": "UK",
    "GR": "Greece",
    "HR": "Croatia",
    "HU": "Hungary",
    "IE": "Ireland",
    "IT": "Italy",
    "LT": "Lithuania",
    "LU": "Luxembourg",
    "LV": "Latvia",
    "ME": "Montenegro",
    "MK": "North Macedonia",
    "NL": "Netherlands",
    "NO": "Norway",
    "PL": "Poland",
    "PT": "Portugal",
    "RO": "Romania",
    "RS": "Serbia",
    "SE": "Sweden",
    "SI": "Slovenia",
    "SK": "Slovakia",
    "XK": "Kosovo",
}

scenarios_steel = [
    "base_eu_regain",
    "base_eu_deindustrial",
    "policy_eu_regain",
    "policy_eu_deindustrial",
]
lhv_ammonia = 5.166  # MWh / t
lhv_methanol = 5.528  # MWh / t
scenarios_steel = ["base_eu_regain", "policy_eu_regain"]


# %% FUNCTIONS

scenario = "base_eu_regain"

cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)
file_path = os.path.join(
    parent_dir, "results_8h", scenario, "networks", "base_s_39___2050.nc"
)
n = pypsa.Network(file_path)

# Steel marginal price
steel_price = (
    n.buses_t.marginal_price.loc[
        :, n.buses_t.marginal_price.columns.str.contains("steel")
    ]
    .mean()
    .iloc[0]
    / 1e3
)  # €/t steel
mean_hist_steel_price = 415  # €/t steel https://tradingeconomics.com/commodity/steel


# Cement marginal price
cement_price = (
    n.buses_t.marginal_price.loc[
        :,
        n.buses_t.marginal_price.columns.str.contains("cement")
        & ~n.buses_t.marginal_price.columns.str.contains("process emissions"),
    ]
    .mean()
    .iloc[0]
    / 1e3
)  # €/t cement
mean_hist_cement_price = 93  # €/t steel https://www.cemnet.com/News/story/175146/cement-prices-in-italy-see-27-rise-in-april.html

# Ammonia marginal price
ammonia_price = (
    n.buses_t.marginal_price.loc[
        :, n.buses_t.marginal_price.columns.str.contains("NH3")
    ]
    .mean()
    .iloc[0]
    * lhv_ammonia
)
mean_hist_ammonia_price = (
    470  # €/t https://pubs.usgs.gov/periodicals/mcs2024/mcs2024-nitrogen.pdf
)

# Methanol marginal price
methanol_price = (
    n.buses_t.marginal_price.loc[
        :,
        n.buses_t.marginal_price.columns.str.contains("methanol")
        & ~n.buses_t.marginal_price.columns.str.contains("industry")
        & ~n.buses_t.marginal_price.columns.str.contains("shipping"),
    ]
    .mean()
    .iloc[0]
    * lhv_methanol
)  # €/t methanol

mean_hist_methanol_price = 326  # €/t https://tradingeconomics.com/commodity/methanol

# HVC marginal price
hvc_prices = (
    n.buses_t.marginal_price.loc[
        :, n.buses_t.marginal_price.columns.str.contains("HVC")
    ]
    .mean()
    .iloc[0]
    / 1e3
)  # €/t


# %%

years = [2030, 2040, 2050]  # Years in the dataset
lhv_ammonia = 5.166  # MWh / t
lhv_methanol = 5.528  # MWh / t

# Invented values for 2020 (for visualization purposes)
hist_2020_prices = {
    "steel": 415,  # €/t
    "cement": 93,  # €/t
    "ammonia": 470,  # €/t
    "methanol": 326,  # €/t
    "HVC": 600,  # €/t
}

# Initialize a dictionary to store DataFrames for each commodity
price_data = {
    commodity: pd.DataFrame(index=scenarios_steel, columns=[2020] + years)
    for commodity in hist_2020_prices.keys()
}

# Fill 2020 values
for commodity, value in hist_2020_prices.items():
    price_data[commodity][2020] = value

max_value = 0
cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)

# Iterate through scenarios and years to calculate prices
for scenario in scenarios_steel:
    for year in years:
        # Load network for the given year
        file_path = os.path.join(
            parent_dir, "results_8h", scenario, "networks", f"base_s_39___{year}.nc"
        )
        n = pypsa.Network(file_path)

        # Steel price
        steel_price = (
            n.buses_t.marginal_price.loc[
                :, n.buses_t.marginal_price.columns.str.contains("steel")
            ]
            .mean()
            .iloc[0]
            / 1e3
        )
        price_data["steel"].loc[scenario, year] = steel_price

        # Cement price
        cement_price = (
            n.buses_t.marginal_price.loc[
                :,
                n.buses_t.marginal_price.columns.str.contains("cement")
                & ~n.buses_t.marginal_price.columns.str.contains("process emissions"),
            ]
            .mean()
            .iloc[0]
            / 1e3
        )
        price_data["cement"].loc[scenario, year] = cement_price

        # Ammonia price
        ammonia_price = (
            n.buses_t.marginal_price.loc[
                :, n.buses_t.marginal_price.columns.str.contains("NH3")
            ]
            .mean()
            .iloc[0]
            * lhv_ammonia
        )
        price_data["ammonia"].loc[scenario, year] = ammonia_price

        # Methanol price
        methanol_price = (
            n.buses_t.marginal_price.loc[
                :,
                n.buses_t.marginal_price.columns.str.contains("methanol")
                & ~n.buses_t.marginal_price.columns.str.contains("industry")
                & ~n.buses_t.marginal_price.columns.str.contains("shipping"),
            ]
            .mean()
            .iloc[0]
            * lhv_methanol
        )
        price_data["methanol"].loc[scenario, year] = methanol_price

        # HVC price
        hvc_price = (
            n.buses_t.marginal_price.loc[
                :, n.buses_t.marginal_price.columns.str.contains("HVC")
            ]
            .mean()
            .iloc[0]
            / 1e3
        )
        price_data["HVC"].loc[scenario, year] = hvc_price
        max_value = max(max(df.max().max() for df in price_data.values()), max_value)


fig, axes = plt.subplots(1, 4, figsize=(12, 4), sharex=True, sharey=True)

commodities = ["steel", "cement", "ammonia", "methanol"]
# Define color mapping
scenario_colors = {"base_eu_regain": "grey", "policy_eu_regain": "green"}

for idx, (commodity, ax) in enumerate(zip(commodities, axes)):
    for i, scenario in enumerate(scenarios_steel):
        # Extract label: first part after splitting by '_', capitalized
        label = scenario.split("_")[0].capitalize() if idx == 3 else None
        color = scenario_colors.get(scenario, None)
        if scenario == "base_eu_regain":
            label = "Baseline"
        ax.plot(
            [2020] + years,
            price_data[commodity].loc[scenario],
            marker="o",
            linestyle="-",
            label=label,
            color=color,
        )

    ax.set_title(f"{commodity.capitalize()} Price")
    ax.set_xticks([2020] + years)
    ax.set_ylim(0, max_value)
    if idx == 0:
        ax.set_ylabel("Price (€/t)")
    ax.grid(True, linestyle="--")

# Show legend only in the last subplot
axes[-1].legend(loc="upper left", frameon=False)

plt.tight_layout()
plt.savefig("./graphs/commodity_prices.png")
plt.show()
