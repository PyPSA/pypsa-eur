"""
Created on Tue Mar 25 12:07:29 2025

@author: Dibella
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pypsa

# Cement

scenarios = [
    "base_eu_regain",
    "policy_eu_regain",
    "base_eu_deindustrial",
    "policy_eu_deindustrial",
]
# Replace ISO2 country codes with full names
country_names = {
    "AL": "Albania",
    "AT": "Austria",
    "BA": "Bosnia",
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

country_totals = {}  # Track total production per country
all_scenario_data = {}  # Store all relevant data for each scenario

# First pass: Collect data across all scenarios
for scenario in scenarios:
    cwd = os.getcwd()
    parent_dir = os.path.dirname(cwd)
    file_path = os.path.join(
        parent_dir, "results", scenario, "networks", "base_s_39___2050.nc"
    )
    n = pypsa.Network(file_path)
    timestep = n.snapshot_weightings.iloc[0, 0]
    threshold = 1

    # Extract cement production
    prod_cement = (
        -n.links_t.p1.filter(like="Cement Plant", axis=1).sum() * timestep / 1e3
    )
    prod_cement.index = prod_cement.index.str[:2]  # Keep only country code
    prod_cement = prod_cement[prod_cement > threshold]
    summed_prod_cement = prod_cement.groupby(prod_cement.index).sum()

    all_countries = summed_prod_cement.index

    # Extract emissions data
    cement_not_captured = (
        -n.links_t.p1.filter(like="cement process emis to atmosphere", axis=1).sum()
        * timestep
    )
    cement_ccs = -n.links_t.p1.filter(like="cement TGR", axis=1).sum() * timestep

    cement_not_captured.index = cement_not_captured.index.str[:2]
    cement_not_captured = cement_not_captured.groupby(cement_not_captured.index).sum()
    cement_ccs.index = cement_ccs.index.str[:2]
    cement_ccs = cement_ccs.groupby(cement_ccs.index).sum()

    # Calculate CCS share
    share_ccs = round(cement_ccs / (cement_ccs + cement_not_captured), 2)
    share_ccs = share_ccs.reindex(summed_prod_cement.index, fill_value=0)

    # Adjusted cement production
    summed_prod_cement_not_captured = summed_prod_cement * (1 - share_ccs)
    summed_prod_cement_captured = summed_prod_cement * share_ccs
    summed_prod_cement_not_captured = summed_prod_cement_not_captured.dropna()
    summed_prod_cement_captured = summed_prod_cement_captured.dropna()

    if not np.allclose(
        summed_prod_cement_not_captured + summed_prod_cement_captured,
        summed_prod_cement,
        atol=1e-6,
    ):
        print(f"Error in scenario {scenario}: Cement data inconsistency detected.")

    # Map country codes to names
    summed_prod_cement_not_captured.index = [
        country_names.get(code, code) for code in summed_prod_cement_not_captured.index
    ]
    summed_prod_cement_captured.index = [
        country_names.get(code, code) for code in summed_prod_cement_captured.index
    ]

    # Store data
    all_scenario_data[scenario] = {
        "summed_prod_cement_not_captured": summed_prod_cement_not_captured,
        "summed_prod_cement_captured": summed_prod_cement_captured,
    }

    # Accumulate total production
    for country in all_countries:
        full_name = country_names.get(country, country)
        country_totals[full_name] = country_totals.get(
            full_name, 0
        ) + summed_prod_cement.get(country, 0)

# Filter out countries with total production <= 1
relevant_countries = [c for c, total in country_totals.items() if total > 1]

# Second pass: Generate plots
fig, axes = plt.subplots(
    2, len(scenarios) // 2, figsize=(20, 12), sharex=True, sharey=True
)
row, col = 0, 0

for i, scenario in enumerate(scenarios):
    summed_prod_cement_not_captured = all_scenario_data[scenario][
        "summed_prod_cement_not_captured"
    ].reindex(relevant_countries, fill_value=0)
    summed_prod_cement_captured = all_scenario_data[scenario][
        "summed_prod_cement_captured"
    ].reindex(relevant_countries, fill_value=0)

    # Sort countries by total production
    total_production = summed_prod_cement_not_captured + summed_prod_cement_captured
    sorted_countries = total_production.sort_values(ascending=False).index

    # Reorder data
    summed_prod_cement_not_captured = summed_prod_cement_not_captured[sorted_countries]
    summed_prod_cement_captured = summed_prod_cement_captured[sorted_countries]

    # Adjust row/col placement
    if i == 2:
        col = 1
        row = 0
    elif i == 4:
        col = 2
        row = 0

    ax = axes[row, col]

    # Plot data
    ax.bar(
        sorted_countries,
        summed_prod_cement_not_captured.values,
        color="gray",
        label="Cement Not Captured",
    )
    ax.bar(
        sorted_countries,
        summed_prod_cement_captured.values,
        color="orange",
        label="Cement Captured",
        bottom=summed_prod_cement_not_captured.values,
    )

    # Formatting
    ax.tick_params(axis="x", labelrotation=90)
    ax.set_ylim(bottom=0)

    if row == 0:
        scenario_title = scenario.split("_")[
            -1
        ].capitalize()  # Extract last part of the scenario
        ax.set_title(f"{scenario_title} production")
        ax.set_facecolor("#CACACE")
    else:
        ax.set_facecolor("#C1D7AE")

    if col == 0 and row == 0:
        ax.set_ylabel("BASELINE\nMt cement/yr")
    elif col == 0 and row == 1:
        ax.set_ylabel("POLICY\nMt cement/yr")

    if row == 0 and col == 1:
        ax.legend()
    row += 1

plt.tight_layout()
os.makedirs("./graphs", exist_ok=True)
plt.savefig("./graphs/cement_production.png")
plt.show()
