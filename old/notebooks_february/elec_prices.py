"""
Created on Fri Mar 21 08:25:57 2025

@author: Dibella
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

scenarios = [
    "base_reg_regain",
    "base_reg_maintain",
    "base_reg_deindustrial",
    "policy_reg_regain",
    "policy_reg_maintain",
    "policy_reg_deindustrial",
]
lhv_ammonia = 5.166  # MWh / t
lhv_methanol = 5.528  # MWh / t


# %% FUNCTIONS

scenario = "policy_eu_regain"

cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)
file_path = os.path.join(
    parent_dir, "results", scenario, "networks", "base_s_39___2050.nc"
)
n = pypsa.Network(file_path)

mprice = n.buses_t.marginal_price
loads = n.loads
mprice_loads = mprice[mprice.columns.intersection(loads.index)]
loads_w_mprice = n.loads_t.p[n.loads_t.p.columns.intersection(mprice_loads.columns)]
total_costs = mprice_loads * loads_w_mprice

# Electricity
elec_total_costs = total_costs[
    total_costs.columns[
        total_costs.columns.str.endswith(" 0")
        | total_costs.columns.str.contains("electricity")
    ]
]
elec_loads_w_mprice = loads_w_mprice[
    loads_w_mprice.columns[
        loads_w_mprice.columns.str.endswith(" 0")
        | loads_w_mprice.columns.str.contains("electricity")
    ]
]

elec_total_costs.columns = elec_total_costs.columns.str[:2]
elec_total_costs = elec_total_costs.T.groupby(level=0).sum().T.sum()

elec_loads_w_mprice.columns = elec_loads_w_mprice.columns.str[:2]
elec_loads_w_mprice = elec_loads_w_mprice.T.groupby(level=0).sum().T.sum()

weigh_aver_elec = elec_total_costs / elec_loads_w_mprice
weigh_aver_elec.index = weigh_aver_elec.index.map(country_names)


# %%


years = [2030, 2040, 2050]  # Years in the dataset
scenarios = [
    "base_reg_regain",
    "base_eu_regain",
    "base_reg_deindustrial",
    "base_eu_deindustrial",
    "policy_reg_regain",
    "policy_eu_regain",
    "policy_reg_deindustrial",
    "policy_eu_deindustrial",
]


# Initialize a dictionary to store DataFrames for each commodity
price_data = {
    country: pd.DataFrame(index=scenarios, columns=years)
    for country in country_names.values()
}


max_value = 0
cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)

# Iterate through scenarios and years to calculate prices
for scenario in scenarios:
    for year in years:
        # Load network for the given year
        file_path = os.path.join(
            parent_dir, "results", scenario, "networks", f"base_s_39___{year}.nc"
        )
        n = pypsa.Network(file_path)

        mprice = n.buses_t.marginal_price
        loads = n.loads
        mprice_loads = mprice[mprice.columns.intersection(loads.index)]
        loads_w_mprice = n.loads_t.p[
            n.loads_t.p.columns.intersection(mprice_loads.columns)
        ]
        total_costs = mprice_loads * loads_w_mprice

        # Electricity
        elec_total_costs = total_costs[
            total_costs.columns[
                total_costs.columns.str.endswith(" 0")
                | total_costs.columns.str.contains("electricity")
            ]
        ]
        elec_loads_w_mprice = loads_w_mprice[
            loads_w_mprice.columns[
                loads_w_mprice.columns.str.endswith(" 0")
                | loads_w_mprice.columns.str.contains("electricity")
            ]
        ]

        elec_total_costs.columns = elec_total_costs.columns.str[:2]
        elec_total_costs = elec_total_costs.T.groupby(level=0).sum().T.sum()

        elec_loads_w_mprice.columns = elec_loads_w_mprice.columns.str[:2]
        elec_loads_w_mprice = elec_loads_w_mprice.T.groupby(level=0).sum().T.sum()

        weigh_aver_elec = elec_total_costs / elec_loads_w_mprice
        weigh_aver_elec.index = weigh_aver_elec.index.map(country_names)

        for country in weigh_aver_elec.index:
            price_data[country].loc[scenario, year] = weigh_aver_elec[country]

        max_value = max(max(df.max().max() for df in price_data.values()), max_value)

# %%

fig, axes = plt.subplots(5, 7, figsize=(15, 5), sharex=True, sharey=True)

row = 0
col = 0
for idx, country in enumerate(country_names.values()):
    if idx == 7:
        row += 1
        col = 0
    elif idx == 14:
        row += 1
        col = 0
    elif idx == 21:
        row += 1
        col = 0
    elif idx == 28:
        row += 1
        col = 0

    ax = axes[row, col]
    col += 1
    for i, scenario in enumerate(scenarios):
        ax.plot(
            years,
            price_data[country].loc[scenario],
            marker="o",
            linestyle="-",
            label=scenario,
        )

    ax.set_title(f"{country}")
    ax.set_xticks(years)
    ax.set_ylim(0, max_value)
    if idx == 0:
        ax.set_ylabel("Price (€/t)")
    ax.grid(True)

# Extra invisible subplot for the legend
ax_legend = axes[-1, -1]
ax_legend.axis("off")

# Creating legend
handles = []
labels = []
for scenario in scenarios:
    strategy = next(key for key in base_colors.keys() if key in scenario)
    color = policy_colors[strategy] if "policy" in scenario else base_colors[strategy]
    handles.append(plt.Line2D([0], [0], color=color, marker="o", linestyle="-"))
    labels.append(scenario)

ax_legend.legend()


plt.tight_layout()
plt.savefig("./graphs/electricity_prices.png")
plt.show()


# %%


years = [2030, 2040, 2050]  # Years in the dataset
scenarios = [
    "base_reg_regain",
    "base_eu_regain",
    "base_reg_deindustrial",
    "base_eu_deindustrial",
    "policy_reg_regain",
    "policy_eu_regain",
    "policy_reg_deindustrial",
    "policy_eu_deindustrial",
]


# Initialize a dictionary to store DataFrames for each commodity
price_data = {
    country: pd.DataFrame(index=scenarios, columns=years)
    for country in country_names.values()
}


max_value = 0
cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)

# Iterate through scenarios and years to calculate prices
for scenario in scenarios:
    for year in years:
        # Load network for the given year
        file_path = os.path.join(
            parent_dir, "results", scenario, "networks", f"base_s_39___{year}.nc"
        )
        n = pypsa.Network(file_path)

        mprice = n.buses_t.marginal_price
        loads = n.loads
        mprice_loads = mprice[mprice.columns.intersection(loads.index)]
        loads_w_mprice = n.loads_t.p[
            n.loads_t.p.columns.intersection(mprice_loads.columns)
        ]
        total_costs = mprice_loads * loads_w_mprice

        # Energy
        ener_total_costs = total_costs[
            total_costs.columns[
                ~total_costs.columns.str.endswith(" 0")
                & ~total_costs.columns.str.contains("electricity")
                & ~total_costs.columns.str.match(r"^[a-z]")
            ]
        ]
        ener_loads_w_mprice = loads_w_mprice[
            loads_w_mprice.columns[
                ~loads_w_mprice.columns.str.endswith(" 0")
                & ~loads_w_mprice.columns.str.contains("electricity")
                & ~total_costs.columns.str.match(r"^[a-z]")
            ]
        ]

        ener_total_costs.columns = ener_total_costs.columns.str[:2]
        ener_total_costs = ener_total_costs.T.groupby(level=0).sum().T.sum()

        ener_loads_w_mprice.columns = ener_loads_w_mprice.columns.str[:2]
        ener_loads_w_mprice = ener_loads_w_mprice.T.groupby(level=0).sum().T.sum()

        weigh_aver_ener = ener_total_costs / ener_loads_w_mprice

        # Remove EU and then distribute the value WEIGHTED AVERAGE
        weigh_aver_ener_eu = weigh_aver_ener[weigh_aver_ener.index == "EU"]
        weigh_aver_ener = weigh_aver_ener[weigh_aver_ener.index != "EU"]

        # Remove EU and then distribute the value LOADS
        ener_loads_w_mprice_eu = ener_loads_w_mprice[ener_loads_w_mprice.index == "EU"]
        ener_loads_w_mprice = ener_loads_w_mprice[ener_loads_w_mprice.index != "EU"]
        ener_loads_shares = ener_loads_w_mprice / ener_loads_w_mprice.sum()
        extra_costs_from_eu_bus = ener_loads_shares * weigh_aver_ener_eu["EU"]

        weigh_aver_ener += extra_costs_from_eu_bus
        weigh_aver_ener.index = weigh_aver_ener.index.map(country_names)

        for country in weigh_aver_elec.index:
            price_data[country].loc[scenario, year] = weigh_aver_elec[country]

        max_value = max(max(df.max().max() for df in price_data.values()), max_value)


# %%

fig, axes = plt.subplots(5, 7, figsize=(15, 5), sharex=True, sharey=True)

row = 0
col = 0
for idx, country in enumerate(country_names.values()):
    if idx == 7:
        row += 1
        col = 0
    elif idx == 14:
        row += 1
        col = 0
    elif idx == 21:
        row += 1
        col = 0
    elif idx == 28:
        row += 1
        col = 0

    ax = axes[row, col]
    col += 1
    for i, scenario in enumerate(scenarios):
        ax.plot(
            years,
            price_data[country].loc[scenario],
            marker="o",
            linestyle="-",
            label=scenario,
        )

    ax.set_title(f"{country}")
    ax.set_xticks(years)
    ax.set_ylim(0, max_value)
    if idx == 0:
        ax.set_ylabel("Price (€/t)")
    ax.grid(True)

# Extra invisible subplot for the legend
ax_legend = axes[-1, -1]
ax_legend.axis("off")

# Creating legend
handles = []
labels = []
for scenario in scenarios:
    strategy = next(key for key in base_colors.keys() if key in scenario)
    color = policy_colors[strategy] if "policy" in scenario else base_colors[strategy]
    handles.append(plt.Line2D([0], [0], color=color, marker="o", linestyle="-"))
    labels.append(scenario)

ax_legend.legend()


plt.tight_layout()
plt.savefig("./graphs/energy_prices.png")
plt.show()

# %%

years = [2030, 2040, 2050]  # Years in the dataset
scenarios = [
    "base_reg_regain",
    "base_eu_regain",
    "base_reg_deindustrial",
    "base_eu_deindustrial",
    "policy_reg_regain",
    "policy_eu_regain",
    "policy_reg_deindustrial",
    "policy_eu_deindustrial",
]

# Initialize a dictionary to store DataFrames for each commodity
tot_costs = pd.DataFrame(index=scenarios, columns=years)


cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)

# Iterate through scenarios and years to calculate prices
for scenario in scenarios:
    for year in years:
        # Load network for the given year
        file_path = os.path.join(
            parent_dir, "results", scenario, "networks", f"base_s_39___{year}.nc"
        )
        n = pypsa.Network(file_path)

        tot_costs.loc[scenario, year] = n.objective / 1e9

# %%
colors = {
    "base_eu_regain": "#4F5050",
    "base_eu_deindustrial": "#85877C",
    "base_reg_deindustrial": "#B0B2A1",
    "base_reg_regain": "grey",
    "policy_eu_regain": "#5D8850",
    "policy_eu_deindustrial": "#95BF74",
    "policy_reg_deindustrial": "#C5DEB1",
    "policy_reg_regain": "green",
}


fig, ax = plt.subplots(1, 1, figsize=(10, 10), sharex=True, sharey=True)

row = 0
col = 0
# ax = axes[row,col]
for scenario in scenarios:
    ax.plot(years, tot_costs.loc[scenario], marker="o", linestyle="-", label=scenario)

    ax.set_xticks(years)
    ax.grid(True)

    ax.legend()


plt.tight_layout()
plt.savefig("./graphs/carbon_and_total_costs.png")
plt.show()
