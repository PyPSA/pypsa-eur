"""
Created on Tue Mar 11 09:13:16 2025

@author: alice
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa

# Define colors
base_colors = {"regain": "#4F5050", "maintain": "#85877C", "deindustrial": "#B0B2A1"}
policy_colors = {"regain": "#5D8850", "maintain": "#95BF74", "deindustrial": "#C5DEB1"}


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

scenarios = [
    "base_reg_regain",
    "policy_reg_regain",
    "base_reg_maintain",
    "policy_reg_maintain",
    "base_reg_deindustrial",
    "policy_reg_deindustrial",
]


# %% FUNCTIONS


# Map scenarios to colors
def get_scenario_color(scenario):
    if "base" in scenario:
        if "regain" in scenario:
            return base_colors["regain"]
        elif "maintain" in scenario:
            return base_colors["maintain"]
        elif "deindustrial" in scenario:
            return base_colors["deindustrial"]
    elif "policy" in scenario:
        if "regain" in scenario:
            return policy_colors["regain"]
        elif "maintain" in scenario:
            return policy_colors["maintain"]
        elif "deindustrial" in scenario:
            return policy_colors["deindustrial"]
    return "white"  # Default background


def share_green_h2(n):
    # Extract H2 Clean and H2 Dirty production
    timestep = n.snapshot_weightings.iloc[0, 0]
    h2_clean = n.links.loc[
        n.links.index.str.contains("H2 Electrolysis|SMR CC", regex=True, na=False), :
    ].index
    h2_dirty = n.links.loc[
        n.links.index.str.contains("SMR(?! CC)", regex=True, na=False), :
    ].index
    h2_clean_df = -n.links_t.p1.loc[:, h2_clean].sum() * timestep
    h2_dirty_df = -n.links_t.p1.loc[:, h2_dirty].sum() * timestep

    h2_clean_df.index = h2_clean_df.index.str[:2]
    h2_clean_df = h2_clean_df.groupby(h2_clean_df.index).sum()

    h2_dirty_df.index = h2_dirty_df.index.str[:2]
    h2_dirty_df = h2_dirty_df.groupby(h2_dirty_df.index).sum()

    # Calculate share of green and grey H2
    share_green = round(h2_clean_df / (h2_clean_df + h2_dirty_df), 2)

    return share_green


# Steel


def plot_steel_scenarios(scenarios):
    country_totals = {}  # Dictionary to track total production per country
    all_scenario_data = {}  # To store all relevant data for each scenario

    # First pass: Collect data across all scenarios
    for scenario in scenarios:
        # Load network once
        cwd = os.getcwd()
        parent_dir = os.path.dirname(cwd)
        file_path = os.path.join(
            parent_dir, "results", scenario, "networks", "base_s_39___2050.nc"
        )
        n = pypsa.Network(file_path)
        timestep = n.snapshot_weightings.iloc[0, 0]

        # Extract EAF production
        p_nom_eaf = (
            -n.links_t.p1.filter(like="EAF", axis=1).sum() * timestep / 1e3
        )  # to Mt
        p_nom_eaf.index = p_nom_eaf.index.str[:2]  # Keep only country code
        p_nom_eaf = p_nom_eaf[p_nom_eaf >= 0]  # Remove values below 0

        # Extract BOF production
        p_nom_bof = (
            -n.links_t.p1.filter(like="BF-BOF", axis=1).sum() * timestep / 1e3
        )  # to Mt
        p_nom_bof.index = p_nom_bof.index.str[:2]  # Keep only country code
        p_nom_bof = p_nom_bof[p_nom_bof >= 0]  # Remove values below 0

        # Sum per country
        summed_p_nom_eaf = p_nom_eaf.groupby(p_nom_eaf.index).sum()
        summed_p_nom_bof = p_nom_bof.groupby(p_nom_bof.index).sum()

        # Get the combined index (union of both indices)
        all_countries = summed_p_nom_eaf.index.union(summed_p_nom_bof.index)

        # Reindex both DataFrames to ensure all countries are present
        summed_p_nom_eaf = summed_p_nom_eaf.reindex(all_countries, fill_value=0)
        summed_p_nom_bof = summed_p_nom_bof.reindex(all_countries, fill_value=0)
        all_countries = [country_names.get(code, code) for code in all_countries]

        # Extract H2 Clean and H2 Dirty production
        share_green = share_green_h2(n)

        # Extract CH4 and H2-based DRI data
        dri_ch4 = (
            -n.links_t.p1.filter(like="CH4 to syn gas DRI", axis=1).sum() * timestep
        )
        dri_h2 = -n.links_t.p1.filter(like="H2 to syn gas DRI", axis=1).sum() * timestep

        # Calculate the share of H2 in DRI production -> at the European level now
        share_h2 = round(dri_h2.sum() / (dri_h2.sum() + dri_ch4.sum()), 2)

        # Calculate the adjusted EAF production
        summed_p_nom_ch4_eaf = summed_p_nom_eaf * (1 - share_h2)  # EAF with CH4
        summed_p_nom_grey_h2_eaf = (
            summed_p_nom_eaf * share_h2 * (1 - share_green)
        )  # EAF with Grey H2
        summed_p_nom_green_h2_eaf = (
            summed_p_nom_eaf * share_h2 * share_green
        )  # EAF with Green H2

        # Check if the sum matches the original summed_eaf
        if not np.allclose(
            summed_p_nom_ch4_eaf + summed_p_nom_grey_h2_eaf + summed_p_nom_green_h2_eaf,
            summed_p_nom_eaf,
            atol=1e-6,
        ):
            print(
                f"Error in scenario {scenario}: The sum of CH4, Grey H2, and Green H2 EAF does not match the original EAF within the threshold."
            )

        # Map country codes to full country names
        summed_p_nom_ch4_eaf.index = [
            country_names.get(code, code) for code in summed_p_nom_ch4_eaf.index
        ]
        summed_p_nom_grey_h2_eaf.index = [
            country_names.get(code, code) for code in summed_p_nom_grey_h2_eaf.index
        ]
        summed_p_nom_green_h2_eaf.index = [
            country_names.get(code, code) for code in summed_p_nom_green_h2_eaf.index
        ]
        summed_p_nom_bof.index = [
            country_names.get(code, code) for code in summed_p_nom_bof.index
        ]

        # Store data for plotting later
        all_scenario_data[scenario] = {
            "summed_p_nom_ch4_eaf": summed_p_nom_ch4_eaf,
            "summed_p_nom_grey_h2_eaf": summed_p_nom_grey_h2_eaf,
            "summed_p_nom_green_h2_eaf": summed_p_nom_green_h2_eaf,
            "summed_p_nom_bof": summed_p_nom_bof,
        }

        # Store total production per country (sum across all scenarios)
        all_countries = [country_names.get(code, code) for code in all_countries]

        for country in all_countries:
            if country not in country_totals:
                country_totals[country] = 0
            country_totals[country] += summed_p_nom_eaf.get(
                country, 0
            ) + summed_p_nom_bof.get(country, 0)

    # Identify countries to keep (remove those with total production = 0 across all scenarios)
    relevant_countries = [c for c, total in country_totals.items() if total > 1]

    # Second pass: Generate plots
    fig, axes = plt.subplots(
        2, len(scenarios) // 2, figsize=(20, 12), sharex=True, sharey=True
    )
    # axes = axes.flatten()  # Convert to 1D array for easy iteration
    row = 0
    col = 0

    for i, scenario in enumerate(scenarios):
        # Get the data for this scenario
        summed_p_nom_ch4_eaf = all_scenario_data[scenario]["summed_p_nom_ch4_eaf"]
        summed_p_nom_grey_h2_eaf = all_scenario_data[scenario][
            "summed_p_nom_grey_h2_eaf"
        ]
        summed_p_nom_green_h2_eaf = all_scenario_data[scenario][
            "summed_p_nom_green_h2_eaf"
        ]
        summed_p_nom_bof = all_scenario_data[scenario]["summed_p_nom_bof"]

        # Reindex using only the relevant countries
        summed_p_nom_ch4_eaf = summed_p_nom_ch4_eaf.reindex(
            relevant_countries, fill_value=0
        )
        summed_p_nom_grey_h2_eaf = summed_p_nom_grey_h2_eaf.reindex(
            relevant_countries, fill_value=0
        )
        summed_p_nom_green_h2_eaf = summed_p_nom_green_h2_eaf.reindex(
            relevant_countries, fill_value=0
        )
        summed_p_nom_bof = summed_p_nom_bof.reindex(relevant_countries, fill_value=0)

        # Sort by total production (EAF + BOF)
        total_production = (
            summed_p_nom_ch4_eaf
            + summed_p_nom_grey_h2_eaf
            + summed_p_nom_green_h2_eaf
            + summed_p_nom_bof
        )
        sorted_countries = total_production.sort_values(ascending=0).index

        # Reorder dataframes
        summed_p_nom_ch4_eaf = summed_p_nom_ch4_eaf[sorted_countries]
        summed_p_nom_grey_h2_eaf = summed_p_nom_grey_h2_eaf[sorted_countries]
        summed_p_nom_green_h2_eaf = summed_p_nom_green_h2_eaf[sorted_countries]
        summed_p_nom_bof = summed_p_nom_bof[sorted_countries]

        if i == 2:
            col = 1
            row = 0
        elif i == 4:
            col = 2
            row = 0

        ax = axes[row, col]

        # Plot in subplot
        ax.bar(sorted_countries, summed_p_nom_bof.values, color="black", label="BF-BOF")
        ax.bar(
            sorted_countries,
            summed_p_nom_ch4_eaf.values,
            color="#552C2D",
            label="CH4 DRI EAF",
            bottom=summed_p_nom_bof.values,
        )
        ax.bar(
            sorted_countries,
            summed_p_nom_grey_h2_eaf.values,
            color="gray",
            label="Grey H2 DRI EAF",
            bottom=summed_p_nom_bof.values + summed_p_nom_ch4_eaf.values,
        )
        ax.bar(
            sorted_countries,
            summed_p_nom_green_h2_eaf.values,
            color="green",
            label="Green H2 DRI EAF",
            bottom=summed_p_nom_bof.values
            + summed_p_nom_ch4_eaf.values
            + summed_p_nom_grey_h2_eaf.values,
        )

        """
        # Set the titles (place them as you requested)
        if col == 0 and row == 0:
            fig.text(0.02, 0.975, 'Baseline', fontsize=15, ha='left')
        elif col == 0 and row == 1:
            fig.text(0.02, 0.513, 'Policy', fontsize=15, ha='left')
        """
        # Show x-tick labels for all subplots
        ax.tick_params(axis="x", labelrotation=90)
        ax.set_ylim(bottom=0)  # Only set ymin, ymax will adjust automatically

        if row == 0:
            scenario_title = scenario.split("_")[
                -1
            ].capitalize()  # Extract last part of the scenario
            ax.set_title(f"{scenario_title} production")
            ax.set_facecolor("#CACACE")
        else:
            ax.set_facecolor("#C1D7AE")

        if col == 0 and row == 0:
            ax.set_ylabel("BASELINE\nMt steel/yr")
        elif col == 0 and row == 1:
            ax.set_ylabel("POLICY\nMt steel/yr")

        ax.grid(axis="x", linestyle="--", alpha=0.7)
        if row == 0 and col == 2:
            ax.legend()

        row += 1

    plt.tight_layout()
    # Check if the folder exists, and create it if it doesn't
    os.makedirs("./graphs", exist_ok=True)
    plt.savefig("./graphs/steel_production_regional.png")
    plt.show()


# Cement


def plot_cement_scenarios(scenarios):
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

        # Extract cement production
        prod_cement = (
            -n.links_t.p1.filter(like="Cement Plant", axis=1).sum() * timestep / 1e3
        )
        prod_cement.index = prod_cement.index.str[:2]  # Keep only country code
        prod_cement = prod_cement[prod_cement > 1e-8]
        summed_prod_cement = prod_cement.groupby(prod_cement.index).sum()

        all_countries = summed_prod_cement.index

        # Extract emissions data
        cement_not_captured = (
            -n.links_t.p1.filter(like="cement process emis to atmosphere", axis=1).sum()
            * timestep
        )
        cement_ccs = -n.links_t.p1.filter(like="cement TGR", axis=1).sum() * timestep

        cement_not_captured.index = cement_not_captured.index.str[:2]
        cement_not_captured = cement_not_captured.groupby(
            cement_not_captured.index
        ).sum()
        cement_ccs.index = cement_ccs.index.str[:2]
        cement_ccs = cement_ccs.groupby(cement_ccs.index).sum()

        # Calculate CCS share
        share_ccs = round(cement_ccs / (cement_ccs + cement_not_captured), 2)

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
            country_names.get(code, code)
            for code in summed_prod_cement_not_captured.index
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
        summed_prod_cement_not_captured = summed_prod_cement_not_captured[
            sorted_countries
        ]
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

        if row == 0 and col == 2:
            ax.legend()
        row += 1

    plt.tight_layout()
    os.makedirs("./graphs", exist_ok=True)
    plt.savefig("./graphs/cement_production_regional.png")
    plt.show()


# Ammonia


def plot_ammonia_scenarios(scenarios):
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

        # Extract ammonia production
        prod_nh3 = -n.links_t.p1.filter(like="Haber-Bosch", axis=1).sum() * timestep
        prod_nh3.index = prod_nh3.index.str[:2]  # Keep only country code
        prod_nh3 = prod_nh3[prod_nh3 > 1e-6]  # Remove negligible values
        summed_prod_nh3 = prod_nh3.groupby(prod_nh3.index).sum()

        all_countries = summed_prod_nh3.index

        # Extract H2 Clean and H2 Dirty production
        share_green = share_green_h2(n)

        # Adjusted ammonia production
        summed_prod_nh3_green = summed_prod_nh3 * share_green
        summed_prod_nh3_grey = summed_prod_nh3 * (1 - share_green)
        summed_prod_nh3_green = summed_prod_nh3_green.dropna()
        summed_prod_nh3_grey = summed_prod_nh3_grey.dropna()

        if not np.allclose(
            summed_prod_nh3_green + summed_prod_nh3_grey, summed_prod_nh3, atol=1e-6
        ):
            print(f"Error in scenario {scenario}: Ammonia data inconsistency detected.")

        # Map country codes to names
        summed_prod_nh3_green.index = [
            country_names.get(code, code) for code in summed_prod_nh3_green.index
        ]
        summed_prod_nh3_grey.index = [
            country_names.get(code, code) for code in summed_prod_nh3_grey.index
        ]

        # Store data
        all_scenario_data[scenario] = {
            "summed_prod_nh3_green": summed_prod_nh3_green,
            "summed_prod_nh3_grey": summed_prod_nh3_grey,
        }

        # Accumulate total production
        for country in all_countries:
            full_name = country_names.get(country, country)
            country_totals[full_name] = country_totals.get(
                full_name, 0
            ) + summed_prod_nh3.get(country, 0)

    # Filter out countries with total production <= 1
    relevant_countries = [c for c, total in country_totals.items() if total > 1]

    # Second pass: Generate plots
    fig, axes = plt.subplots(
        2, len(scenarios) // 2, figsize=(20, 12), sharex=True, sharey=True
    )
    row, col = 0, 0

    for i, scenario in enumerate(scenarios):
        summed_prod_nh3_green = all_scenario_data[scenario][
            "summed_prod_nh3_green"
        ].reindex(relevant_countries, fill_value=0)
        summed_prod_nh3_grey = all_scenario_data[scenario][
            "summed_prod_nh3_grey"
        ].reindex(relevant_countries, fill_value=0)

        # Sort countries by total production
        total_production = summed_prod_nh3_green + summed_prod_nh3_grey
        sorted_countries = total_production.sort_values(ascending=False).index

        # Reorder data
        summed_prod_nh3_green = summed_prod_nh3_green[sorted_countries]
        summed_prod_nh3_grey = summed_prod_nh3_grey[sorted_countries]

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
            summed_prod_nh3_grey.values,
            color="#552C2D",
            label="Grey H2 Ammonia",
        )
        ax.bar(
            sorted_countries,
            summed_prod_nh3_green.values,
            color="green",
            label="Green H2 Ammonia",
            bottom=summed_prod_nh3_grey.values,
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
            ax.set_ylabel("BASELINE\nMt NH3/yr")
        elif col == 0 and row == 1:
            ax.set_ylabel("POLICY\nMt NH3/yr")

        if row == 0 and col == 2:
            ax.legend()
        row += 1

    plt.tight_layout()
    os.makedirs("./graphs", exist_ok=True)
    plt.savefig("./graphs/ammonia_production_regional.png")
    plt.show()


# Methanol
def plot_methanol_scenarios(scenarios):
    country_totals = {}  # Track total production per country
    all_scenario_data = {}  # Store data for each scenario

    for scenario in scenarios:
        cwd = os.getcwd()
        parent_dir = os.path.dirname(cwd)
        file_path = os.path.join(
            parent_dir, "results", scenario, "networks", "base_s_39___2050.nc"
        )
        n = pypsa.Network(file_path)
        timestep = n.snapshot_weightings.iloc[0, 0]

        # Extract methanol production from methanolisation plants
        prod_methanol = (
            -n.links_t.p1.filter(like="methanolisation", axis=1).sum() * timestep
        )
        prod_methanol.index = prod_methanol.index.str[:2]  # Keep only country code
        prod_methanol = prod_methanol[prod_methanol >= 0]  # Remove negative values
        prod_methanol = prod_methanol.groupby(prod_methanol.index).sum()

        all_countries = prod_methanol.index

        # Extract H2 Clean and H2 Dirty production

        share_green = share_green_h2(n)

        # Adjust methanol production based on H2 type
        prod_methanol_green = prod_methanol * share_green
        prod_methanol_grey = prod_methanol * (1 - share_green)

        if not np.allclose(
            prod_methanol_green + prod_methanol_grey, prod_methanol, atol=1e-6
        ):
            print(
                f"Error in scenario {scenario}: The sum of Grey H2 and Green H2 methanol does not match the original methanol production."
            )

        all_scenario_data[scenario] = {
            "prod_methanol_grey": prod_methanol_grey,
            "prod_methanol_green": prod_methanol_green,
        }

        for country in all_countries:
            if country not in country_totals:
                country_totals[country] = 0
            country_totals[country] += prod_methanol_green.get(
                country, 0
            ) + prod_methanol_grey.get(country, 0)

    relevant_countries = [c for c, total in country_totals.items() if total > 1]

    fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, scenario in enumerate(scenarios):
        prod_methanol_green = all_scenario_data[scenario][
            "prod_methanol_green"
        ].reindex(relevant_countries, fill_value=0)
        prod_methanol_grey = all_scenario_data[scenario]["prod_methanol_grey"].reindex(
            relevant_countries, fill_value=0
        )

        total_production = prod_methanol_green + prod_methanol_grey
        sorted_countries = total_production.sort_values().index

        prod_methanol_green = prod_methanol_green[sorted_countries]
        prod_methanol_grey = prod_methanol_grey[sorted_countries]

        ax = axes[i]
        ax.bar(
            sorted_countries,
            prod_methanol_grey.values,
            color="brown",
            label="Methanolisation with grey H2",
        )
        ax.bar(
            sorted_countries,
            prod_methanol_green.values,
            color="green",
            label="Methanolisation with green H2",
            left=prod_methanol_grey.values,
        )

        ax.set_xlabel("MWh Methanol/yr")
        ax.set_ylabel("Country")
        ax.set_title(f"Scenario: {scenario}")
        ax.grid(axis="x", linestyle="--", alpha=0.7)
        ax.legend()

    plt.tight_layout()
    os.makedirs("./graphs", exist_ok=True)
    plt.savefig("./graphs/methanol_production_regional.png")
    plt.show()


def plot_methanol_scenarios(scenarios):
    country_totals = {}
    all_scenario_data = {}

    for scenario in scenarios:
        cwd = os.getcwd()
        parent_dir = os.path.dirname(cwd)
        file_path = os.path.join(
            parent_dir, "results", scenario, "networks", "base_s_39___2050.nc"
        )
        n = pypsa.Network(file_path)
        timestep = n.snapshot_weightings.iloc[0, 0]

        prod_methanol = (
            -n.links_t.p1.filter(like="methanolisation", axis=1).sum() * timestep
        )
        prod_methanol.index = prod_methanol.index.str[:2]
        prod_methanol = prod_methanol[prod_methanol > 1e-10]
        prod_methanol = prod_methanol.groupby(prod_methanol.index).sum()

        all_countries = prod_methanol.index
        all_countries = [country_names.get(code, code) for code in all_countries]

        share_green = share_green_h2(n)
        prod_methanol_green = prod_methanol * share_green
        prod_methanol_grey = prod_methanol * (1 - share_green)

        prod_methanol_green.index = [
            country_names.get(code, code) for code in prod_methanol_green.index
        ]
        prod_methanol_grey.index = [
            country_names.get(code, code) for code in prod_methanol_grey.index
        ]

        all_scenario_data[scenario] = {
            "methanol_green": prod_methanol_green,
            "methanol_grey": prod_methanol_grey,
        }

        for country in all_countries:
            full_name = country_names.get(country, country)
            country_totals[full_name] = country_totals.get(full_name, 0) + (
                prod_methanol_green.get(country, 0) + prod_methanol_grey.get(country, 0)
            )

    relevant_countries = [c for c, total in country_totals.items() if total > 1]

    fig, axes = plt.subplots(
        2, len(scenarios) // 2, figsize=(20, 12), sharex=True, sharey=True
    )
    row, col = 0, 0

    for i, scenario in enumerate(scenarios):
        methanol_green = all_scenario_data[scenario]["methanol_green"].reindex(
            relevant_countries, fill_value=0
        )
        methanol_grey = all_scenario_data[scenario]["methanol_grey"].reindex(
            relevant_countries, fill_value=0
        )

        total_production = methanol_green + methanol_grey
        sorted_countries = total_production.sort_values(ascending=False).index

        if i == 2:
            col = 1
            row = 0
        elif i == 4:
            col = 2
            row = 0

        ax = axes[row, col]

        ax.bar(
            sorted_countries,
            methanol_grey.values,
            color="grey",
            label="Methanol with grey H2",
        )
        ax.bar(
            sorted_countries,
            methanol_green.values,
            color="green",
            label="Methanol with green H2",
            bottom=methanol_grey.values,
        )

        ax.tick_params(axis="x", labelrotation=90)
        ax.set_ylim(bottom=0)
        scenario_title = scenario.split("_")[-1].capitalize()
        ax.set_title(f"{scenario_title} Production")

        if col == 0 and row == 0:
            ax.set_ylabel("BASELINE\nkt Methanol/yr")
        elif col == 0 and row == 1:
            ax.set_ylabel("POLICY\nkt Methanol/yr")

        ax.grid(axis="x", linestyle="--", alpha=0.7)
        if row == 0 and col == 2:
            ax.legend()
        row += 1

    plt.tight_layout()
    os.makedirs("./graphs", exist_ok=True)
    plt.savefig("./graphs/methanol_production_regional.png")
    plt.show()


# Validation


def plot_capacities(csv_file):
    # Load the CSV file
    df = pd.read_csv(csv_file, index_col=0)

    # Extract the first two letters of the index (country code) and group by it
    df.index = df.index.str[:2]
    df_grouped = df.groupby(df.index).sum()

    # Combine EAF and DRI + EAF under a single category 'EAF'
    if "EAF" in df_grouped.columns and "DRI + EAF" in df_grouped.columns:
        df_grouped["EAF"] = df_grouped["EAF"] + df_grouped["DRI + EAF"]
        df_grouped.drop(columns=["DRI + EAF"], inplace=True)

    # Apply the multipliers
    if "Ammonia" in df_grouped.columns:
        df_grouped["Ammonia"] *= 5166
    if "Methanol" in df_grouped.columns:
        df_grouped["Methanol"] *= 5528

    # Remove countries with values below 0.1
    df_grouped = df_grouped.loc[(df_grouped > 0.1).any(axis=1)]

    # Define the columns for the stacked bar plot
    stacked_columns = ["Integrated steelworks", "EAF"]
    df_filtered = df_grouped[stacked_columns].dropna()
    df_filtered = df_filtered.sort_values(by=stacked_columns, ascending=False)

    # Remove EAF and Integrated Steelworks from individual plotting
    remaining_columns = [
        col for col in df_grouped.columns if col not in stacked_columns
    ]

    # Create the figure
    fig, axes = plt.subplots(
        nrows=len(remaining_columns) + 1, figsize=(10, 5 * (len(remaining_columns) + 1))
    )

    # Plot the stacked bar for Integrated Steelworks and EAF
    df_filtered.plot(
        kind="bar",
        stacked=True,
        ax=axes[0],
        width=1,
        color=["steelblue", "darkorange"],
        edgecolor="black",
    )
    axes[0].set_title("Steel production")
    axes[0].set_xlabel("Country Code")
    axes[0].set_ylabel("kt steel")
    axes[0].legend(title="Production Type")
    axes[0].grid(axis="y", linestyle="--", alpha=0.7)

    # Plot the remaining columns separately
    for i, column in enumerate(remaining_columns, start=1):
        df_sorted = df_grouped[column].sort_values(ascending=False)
        df_sorted.plot(
            kind="bar", ax=axes[i], color="skyblue", edgecolor="black", width=1
        )
        axes[i].set_title(column)
        axes[i].set_xlabel("Country Code")

        axes[i].grid(axis="y", linestyle="--", alpha=0.7)

    axes[1].set_ylabel("kt cement")
    axes[2].set_ylabel("MWh NH3")
    axes[3].set_ylabel("kt HVC")
    axes[4].set_ylabel("MWh MeOH")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig("./graphs/existing_capacities.png")
    plt.show()


def plot_hvc_scenarios(scenarios):
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

        # Extract HVC production
        prod_naphtha_cracker = (
            -n.links_t.p1.filter(like="naphtha steam cracker", axis=1).sum() * timestep
        )
        prod_naphtha_cracker.index = prod_naphtha_cracker.index.str[
            :2
        ]  # Keep only country code
        prod_naphtha_cracker = prod_naphtha_cracker[prod_naphtha_cracker > 1e-10]
        prod_naphtha_cracker = prod_naphtha_cracker.groupby(
            prod_naphtha_cracker.index
        ).sum()

        prod_methanol_hvc = (
            -n.links_t.p1.filter(like="methanol-to-hvc", axis=1).sum() * timestep
        )
        prod_methanol_hvc.index = prod_methanol_hvc.index.str[:2]
        prod_methanol_hvc = prod_methanol_hvc[prod_methanol_hvc > 1e-10]
        prod_methanol_hvc = prod_methanol_hvc.groupby(prod_methanol_hvc.index).sum()

        all_countries = prod_naphtha_cracker.index.union(prod_methanol_hvc.index)
        all_countries = [country_names.get(code, code) for code in all_countries]

        # Extract emissions data
        oil_prod_biomass = (
            -n.links_t.p1.filter(like="biomass to liquid", axis=1).sum() * timestep
        )
        oil_prod_biomass.index = oil_prod_biomass.index.str[:2]
        oil_prod_biomass = (
            oil_prod_biomass[oil_prod_biomass > 0].groupby(oil_prod_biomass.index).sum()
        )

        oil_prod = -n.links_t.p1.filter(like="EU oil", axis=1).sum() * timestep
        oil_prod.index = oil_prod.index.str[:2]
        oil_prod = oil_prod[oil_prod > 0].groupby(oil_prod.index).sum()

        share_bio_oil = round(oil_prod_biomass / (oil_prod_biomass + oil_prod), 2)
        share_green = share_green_h2(n)

        prod_naphtha_bio = prod_naphtha_cracker * share_bio_oil.iloc[0]
        prod_naphtha_grey = prod_naphtha_cracker * (1 - share_bio_oil.iloc[0])
        prod_methanol_hvc_green = prod_methanol_hvc * share_green
        prod_methanol_hvc_grey = prod_methanol_hvc * (1 - share_green)
        prod_methanol_hvc_green = prod_methanol_hvc_green.dropna()
        prod_methanol_hvc_grey = prod_methanol_hvc_grey.dropna()

        # Map country codes to full country names
        prod_naphtha_bio.index = [
            country_names.get(code, code) for code in prod_naphtha_bio.index
        ]
        prod_naphtha_grey.index = [
            country_names.get(code, code) for code in prod_naphtha_grey.index
        ]
        prod_methanol_hvc_green.index = [
            country_names.get(code, code) for code in prod_methanol_hvc_green.index
        ]
        prod_methanol_hvc_grey.index = [
            country_names.get(code, code) for code in prod_methanol_hvc_grey.index
        ]

        all_scenario_data[scenario] = {
            "hvc_fossil_naphtha": prod_naphtha_grey,
            "hvc_bio_naphtha": prod_naphtha_bio,
            "hvc_methanol_green": prod_methanol_hvc_green,
            "hvc_methanol_grey": prod_methanol_hvc_grey,
        }

        for country in all_countries:
            full_name = country_names.get(country, country)
            country_totals[full_name] = country_totals.get(full_name, 0) + (
                prod_naphtha_grey.get(country, 0)
                + prod_naphtha_bio.get(country, 0)
                + prod_methanol_hvc_green.get(country, 0)
                + prod_methanol_hvc_grey.get(country, 0)
            )

    relevant_countries = [c for c, total in country_totals.items() if total > 1]

    # Second pass: Generate plots
    fig, axes = plt.subplots(
        2, len(scenarios) // 2, figsize=(20, 12), sharex=True, sharey=True
    )
    row, col = 0, 0

    for i, scenario in enumerate(scenarios):
        hvc_fossil_naphtha = all_scenario_data[scenario]["hvc_fossil_naphtha"].reindex(
            relevant_countries, fill_value=0
        )
        hvc_bio_naphtha = all_scenario_data[scenario]["hvc_bio_naphtha"].reindex(
            relevant_countries, fill_value=0
        )
        hvc_methanol_grey = all_scenario_data[scenario]["hvc_methanol_grey"].reindex(
            relevant_countries, fill_value=0
        )
        hvc_methanol_green = all_scenario_data[scenario]["hvc_methanol_green"].reindex(
            relevant_countries, fill_value=0
        )

        total_production = (
            hvc_fossil_naphtha
            + hvc_bio_naphtha
            + hvc_methanol_grey
            + hvc_methanol_green
        )
        sorted_countries = total_production.sort_values(ascending=False).index

        if i == 2:
            col = 1
            row = 0
        elif i == 4:
            col = 2
            row = 0

        ax = axes[row, col]

        ax.bar(
            sorted_countries,
            hvc_fossil_naphtha.values,
            color="black",
            label="Fossil Naphtha Cracker",
        )
        ax.bar(
            sorted_countries,
            hvc_bio_naphtha.values,
            color="brown",
            label="Bio Naphtha Cracker",
            bottom=hvc_fossil_naphtha.values,
        )
        ax.bar(
            sorted_countries,
            hvc_methanol_grey.values,
            color="grey",
            label="Methanol Cracker with grey H2",
            bottom=hvc_fossil_naphtha.values + hvc_bio_naphtha.values,
        )
        ax.bar(
            sorted_countries,
            hvc_methanol_green.values,
            color="green",
            label="Methanol Cracker with green H2",
            bottom=hvc_fossil_naphtha.values
            + hvc_bio_naphtha.values
            + hvc_methanol_grey.values,
        )

        ax.tick_params(axis="x", labelrotation=90)
        ax.set_ylim(bottom=0)
        scenario_title = scenario.split("_")[-1].capitalize()
        ax.set_title(f"{scenario_title} Production")

        if col == 0 and row == 0:
            ax.set_ylabel("BASELINE\nkt HVC/yr")
        elif col == 0 and row == 1:
            ax.set_ylabel("POLICY\nkt HVC/yr")

        ax.grid(axis="x", linestyle="--", alpha=0.7)
        if row == 0 and col == 2:
            ax.legend()
        row += 1

    plt.tight_layout()
    os.makedirs("./graphs", exist_ok=True)
    plt.savefig("./graphs/hvc_production_regional.png")
    plt.show()


# %% GRAPHS

plot_steel_scenarios(scenarios)
plot_cement_scenarios(scenarios)
plot_ammonia_scenarios(scenarios)

plot_hvc_scenarios(scenarios)
plot_methanol_scenarios(scenarios)
