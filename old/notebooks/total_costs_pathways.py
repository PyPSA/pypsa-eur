"""
Created on Fri Dec  6 10:42:50 2024

@author: alice
"""

import matplotlib.pyplot as plt
import pandas as pd
import pypsa

years = [2030, 2040, 2050]
scenarios = [
    "baseline_eu_dem",
    "policy_eu_dem",
    "baseline_regional_dem",
    "policy_regional_dem",
]
costs = pd.DataFrame(0, index=scenarios, columns=years)
costs_stat = pd.DataFrame(0, index=scenarios, columns=years)

for year in years:
    for scenario in scenarios:
        n = pypsa.Network(
            f"../results_24h/{scenario}/postnetworks/base_s_39_lvopt___{year}.nc"
        )
        costs.loc[scenario, year] = n.objective / 1e9  # b€

for year in years:
    for scenario in scenarios:
        n = pypsa.Network(
            f"../results_24h/{scenario}/postnetworks/base_s_39_lvopt___{year}.nc"
        )
        capex = n.statistics.capex().sum()
        opex = n.statistics.opex(aggregate_time="sum").sum()
        costs_stat.loc[scenario, year] = (capex + opex) / 1e9  # b€


# %%

colors = ["red", "blue", "green", "grey"]

plt.figure(figsize=(10, 6))  # Adjust size as needed

# Plot each row as a separate line
for idx, (row_label, row_values) in enumerate(costs.iterrows()):
    print(idx)
    plt.plot(
        costs_stat.columns, row_values, label=row_label, color=colors[idx]
    )  # Assign custom color

# Customize the plot
plt.xlabel("Columns (X-axis)", fontsize=12)
plt.ylabel("Costs (b€)", fontsize=12)
plt.title("Costs Over Columns", fontsize=14)
plt.legend(title="Rows", fontsize=10, bbox_to_anchor=(1.05, 1), loc="upper left")
plt.grid(True)
plt.tight_layout()

plt.savefig("graphs/total_costs_stat_pathways.png", bbox_inches="tight")
plt.show()


# %%

steel_techs
mapping = {
    "Italy": "Europe",
    "France": "Europe",
    "Argentina": "Latin America",
    "Japan": "Asia",
}

# Create a new column based on the mapping

for year in years:
    for scenario in scenarios:
        df = pd.read_csv(f"../results_24h/{scenario}/csvs/costs.csv")
        df = df.dropna(subset=["Unnamed: 2"])
        df["map"] = df["Unnamed: 2"].map(mapping)

        var = data["Unnamed: 2"].unique()
