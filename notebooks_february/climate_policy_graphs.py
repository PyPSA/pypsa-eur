# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 10:42:32 2025

@author: Dibella
"""


import pypsa
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np

def share_green_h2(n):
    # Extract H2 Clean and H2 Dirty production
    timestep = n.snapshot_weightings.iloc[0,0]
    h2_clean = n.links.loc[n.links.index.str.contains('H2 Electrolysis|SMR CC', regex=True, na=False), :].index
    h2_dirty = n.links.loc[n.links.index.str.contains('SMR(?! CC)', regex=True, na=False), :].index
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
        file_path = os.path.join(parent_dir, "results_march", scenario, "networks", "base_s_39___2030.nc")
        n = pypsa.Network(file_path)
        timestep = n.snapshot_weightings.iloc[0, 0]
        threshold = 1 # Remove values below 1 kt/yr

        # Extract EAF production
        p_nom_eaf = -n.links_t.p1.filter(like="EAF", axis=1).sum() * timestep / 1e3 # to Mt
        p_nom_eaf.index = p_nom_eaf.index.str[:2]  # Keep only country code
        p_nom_eaf = p_nom_eaf[p_nom_eaf >= threshold]  # Remove values below 0

        # Extract BOF production
        p_nom_bof = -n.links_t.p1.filter(like="BF-BOF", axis=1).sum() * timestep / 1e3 # to Mt
        p_nom_bof.index = p_nom_bof.index.str[:2]  # Keep only country code
        p_nom_bof = p_nom_bof[p_nom_bof >= threshold]  # Remove values below 0

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
        share_green = share_green.reindex(summed_p_nom_eaf.index, fill_value=0)

        # Extract CH4 and H2-based DRI data
        dri_ch4 = -n.links_t.p1.filter(like='CH4 to syn gas DRI', axis=1).sum() * timestep
        dri_h2 = -n.links_t.p1.filter(like='H2 to syn gas DRI', axis=1).sum() * timestep

        # Calculate the share of H2 in DRI production -> at the European level now
        share_h2 = round(dri_h2.sum() / (dri_h2.sum() + dri_ch4.sum()), 2)

        # Calculate the adjusted EAF production
        summed_p_nom_ch4_eaf = summed_p_nom_eaf * (1 - share_h2)  # EAF with CH4
        summed_p_nom_grey_h2_eaf = summed_p_nom_eaf * share_h2 * (1 - share_green)  # EAF with Grey H2
        summed_p_nom_green_h2_eaf = summed_p_nom_eaf * share_h2 * share_green  # EAF with Green H2
        

        # Check if the sum matches the original summed_eaf
        if not np.allclose(summed_p_nom_ch4_eaf + summed_p_nom_grey_h2_eaf + summed_p_nom_green_h2_eaf, summed_p_nom_eaf, atol=1e-6):
            print(f"Error in scenario {scenario}: The sum of CH4, Grey H2, and Green H2 EAF does not match the original EAF within the threshold.")


        # Map country codes to full country names
        summed_p_nom_ch4_eaf.index = [country_names.get(code, code) for code in summed_p_nom_ch4_eaf.index]
        summed_p_nom_grey_h2_eaf.index = [country_names.get(code, code) for code in summed_p_nom_grey_h2_eaf.index]
        summed_p_nom_green_h2_eaf.index = [country_names.get(code, code) for code in summed_p_nom_green_h2_eaf.index]
        summed_p_nom_bof.index = [country_names.get(code, code) for code in summed_p_nom_bof.index]
 
        # Store data for plotting later
        all_scenario_data[scenario] = {
            "summed_p_nom_ch4_eaf": summed_p_nom_ch4_eaf,
            "summed_p_nom_grey_h2_eaf": summed_p_nom_grey_h2_eaf,
            "summed_p_nom_green_h2_eaf": summed_p_nom_green_h2_eaf,
            "summed_p_nom_bof": summed_p_nom_bof
        }

        for country in all_countries:
            if country not in country_totals:
                country_totals[country] = 0
            country_totals[country] += summed_p_nom_ch4_eaf.get(country, 0) + summed_p_nom_grey_h2_eaf.get(country, 0) + summed_p_nom_green_h2_eaf.get(country, 0) + summed_p_nom_bof.get(country, 0)

    # Identify countries to keep (remove those with total production = 0 across all scenarios)
    relevant_countries = [c for c, total in country_totals.items() if total > 1]

    # Second pass: Generate plots
    fig, axes = plt.subplots(1, 2, figsize=(20, 12), sharex=True, sharey=True)
    #axes = axes.flatten()  # Convert to 1D array for easy iteration
    row = 0
    col = 0


    for i, scenario in enumerate(scenarios):
        print(i)
        # Get the data for this scenario
        summed_p_nom_ch4_eaf = all_scenario_data[scenario]["summed_p_nom_ch4_eaf"]
        summed_p_nom_grey_h2_eaf = all_scenario_data[scenario]["summed_p_nom_grey_h2_eaf"]
        summed_p_nom_green_h2_eaf = all_scenario_data[scenario]["summed_p_nom_green_h2_eaf"]
        summed_p_nom_bof = all_scenario_data[scenario]["summed_p_nom_bof"]


        # Reindex using only the relevant countries
        summed_p_nom_ch4_eaf = summed_p_nom_ch4_eaf.reindex(relevant_countries, fill_value=0)
        summed_p_nom_grey_h2_eaf = summed_p_nom_grey_h2_eaf.reindex(relevant_countries, fill_value=0)
        summed_p_nom_green_h2_eaf = summed_p_nom_green_h2_eaf.reindex(relevant_countries, fill_value=0)
        summed_p_nom_bof = summed_p_nom_bof.reindex(relevant_countries, fill_value=0)

        
        # Sort by total production (EAF + BOF)
        total_production = summed_p_nom_ch4_eaf + summed_p_nom_grey_h2_eaf + summed_p_nom_green_h2_eaf + summed_p_nom_bof
        sorted_countries = total_production.sort_values(ascending = 0).index

        # Reorder dataframes
        summed_p_nom_ch4_eaf = summed_p_nom_ch4_eaf[sorted_countries]
        summed_p_nom_grey_h2_eaf = summed_p_nom_grey_h2_eaf[sorted_countries]
        summed_p_nom_green_h2_eaf = summed_p_nom_green_h2_eaf[sorted_countries]
        summed_p_nom_bof = summed_p_nom_bof[sorted_countries]
        
        ax = axes[col]
        col+=1

        # Plot in subplot
        ax.bar(sorted_countries, summed_p_nom_bof.values, color="black", label="BF-BOF")
        ax.bar(sorted_countries, summed_p_nom_ch4_eaf.values, color="#552C2D", label="CH4 DRI EAF", bottom=summed_p_nom_bof.values)
        ax.bar(sorted_countries, summed_p_nom_grey_h2_eaf.values, color="gray", label="Grey H2 DRI EAF", bottom=summed_p_nom_bof.values + summed_p_nom_ch4_eaf.values)
        ax.bar(sorted_countries, summed_p_nom_green_h2_eaf.values, color="green", label="Green H2 DRI EAF", bottom=summed_p_nom_bof.values + summed_p_nom_ch4_eaf.values + summed_p_nom_grey_h2_eaf.values)

        """
        # Set the titles (place them as you requested)
        if col == 0 and row == 0:
            fig.text(0.02, 0.975, 'Baseline', fontsize=15, ha='left')
        elif col == 0 and row == 1:
            fig.text(0.02, 0.513, 'Policy', fontsize=15, ha='left')
        """
                # Show x-tick labels for all subplots
        ax.tick_params(axis='x', labelrotation=90)
        ax.set_ylim(bottom=0)  # Only set ymin, ymax will adjust automatically
        
        if row == 0:
            scenario_title = scenario.split('_')[-1].capitalize()  # Extract last part of the scenario
            ax.set_title(f"{scenario_title} production", size=14)
            ax.set_facecolor("#CACACE")
        else:
            ax.set_facecolor("#C1D7AE")
            
        if col == 0 and row == 0:
            ax.set_ylabel("BASELINE\nMt steel/yr")
        elif col == 0 and row == 1:
            ax.set_ylabel("POLICY\nMt steel/yr")

        
        ax.grid(axis="x", linestyle="--", alpha=0.7)
        if row == 0 and col == 1:
            ax.legend()
            
        

    plt.tight_layout()
    # Check if the folder exists, and create it if it doesn't
    os.makedirs("./graphs", exist_ok=True)
    plt.savefig("./graphs/steel_production_climate_policy.png")
    plt.show()
    
# %%


def plot_steel_scenarios_pathways(years, scenarios):
    all_scenario_data = {}  # Store scenario and year data
    all_countries = set()
    country_totals = {}
    
    for scenario in scenarios:
        for year in years:
            file_path = os.path.join("..", "results_march", scenario, "networks", f"base_s_39___{year}.nc")
            n = pypsa.Network(file_path)
            timestep = n.snapshot_weightings.iloc[0, 0]
            
            # Extract EAF & BOF production
            p_nom_eaf = -n.links_t.p1.filter(like="EAF").sum() * timestep / 1e3  # Convert to Mt
            p_nom_eaf.index = p_nom_eaf.index.str[:2]
            
            p_nom_bof = -n.links_t.p1.filter(like="BF-BOF").sum() * timestep / 1e3
            p_nom_bof.index = p_nom_bof.index.str[:2]
            
            # Sum per country
            summed_p_nom_eaf = p_nom_eaf.groupby(p_nom_eaf.index).sum()
            summed_p_nom_bof = p_nom_bof.groupby(p_nom_bof.index).sum()
            
            # Extract hydrogen-related data
            share_green = share_green_h2(n).reindex(summed_p_nom_eaf.index, fill_value=0)
            dri_ch4 = -n.links_t.p1.filter(like='CH4 to syn gas DRI').sum() * timestep
            dri_h2 = -n.links_t.p1.filter(like='H2 to syn gas DRI').sum() * timestep
            share_h2 = round(dri_h2.sum() / (dri_h2.sum() + dri_ch4.sum()), 2)
            
            # Split EAF production
            summed_p_nom_ch4_eaf = summed_p_nom_eaf * (1 - share_h2)
            summed_p_nom_grey_h2_eaf = summed_p_nom_eaf * share_h2 * (1 - share_green)
            summed_p_nom_green_h2_eaf = summed_p_nom_eaf * share_h2 * share_green
            
            # Store results
            all_scenario_data[(scenario, year)] = {
                "ch4_eaf": summed_p_nom_ch4_eaf,
                "grey_h2_eaf": summed_p_nom_grey_h2_eaf,
                "green_h2_eaf": summed_p_nom_green_h2_eaf,
                "bof": summed_p_nom_bof
            }
            
            for country in summed_p_nom_eaf.index.union(summed_p_nom_bof.index):
                country_totals[country] = country_totals.get(country, 0) + summed_p_nom_eaf.get(country, 0) + summed_p_nom_bof.get(country, 0)
    
    # Filter out countries with total production = 0
    relevant_countries = {code: country_names.get(code, code) for code, total in country_totals.items() if total > 0}
    
    # Plot
    fig, axes = plt.subplots(2, len(years), figsize=(20, 12), sharex=True, sharey=True)
    
    for row, scenario in enumerate(scenarios):
        for col, year in enumerate(years):
            ax = axes[row, col]
            data = all_scenario_data[(scenario, year)]
            
            # Reindex for consistent order
            for key in data:
                data[key] = data[key].reindex(relevant_countries.keys(), fill_value=0)
            
            # Sorting
            total_production = sum(data.values())
            sorted_countries = [relevant_countries[code] for code in total_production.sort_values(ascending=False).index]
            
            # Plot bars
            bottom = np.zeros(len(sorted_countries))
            for label, color in zip(["bof", "ch4_eaf", "grey_h2_eaf", "green_h2_eaf"], ["black", "#552C2D", "gray", "green"]):
                ax.bar(sorted_countries, data[label].values, bottom=bottom, color=color, label=label.replace("_", " ").upper())
                bottom += data[label].values
            
            ax.set_title(f"{scenario.capitalize()} {year}")
            ax.set_xticklabels(sorted_countries, rotation=90)
            if col == 0:
                ax.set_ylabel("Mt steel/yr")
            if row == 0 and col == len(years) - 1:
                ax.legend()
    
    plt.tight_layout()
    os.makedirs("./graphs", exist_ok=True)
    plt.savefig("./graphs/steel_production_climate_policy.png")
    plt.show()



# %%



# CHECKING THE CLIMATE POLICY

# Replace ISO2 country codes with full names
country_names = {
    "AL": "Albania", "AT": "Austria", "BA": "Bosnia", "BE": "Belgium", "BG": "Bulgaria",
    "CH": "Switzerland", "CZ": "Czechia", "DE": "Germany", "DK": "Denmark", "EE": "Estonia", "ES": "Spain",
    "FI": "Finland", "FR": "France", "GB": "UK", "GR": "Greece", "HR": "Croatia", "HU": "Hungary",
    "IE": "Ireland", "IT": "Italy", "LT": "Lithuania", "LU": "Luxembourg", "LV": "Latvia", "ME": "Montenegro",
    "MK": "North Macedonia", "NL": "Netherlands", "NO": "Norway", "PL": "Poland", "PT": "Portugal", "RO": "Romania",
    "RS": "Serbia", "SE": "Sweden", "SI": "Slovenia", "SK": "Slovakia", "XK": "Kosovo"
}
scenarios = ["base_eu_regain", "policy_eu_regain"]

commodities = ['steel','cement','NH3','HVC','industry methanol']
df = pd.DataFrame(index=commodities,columns = scenarios)

# First pass: Collect data across all scenarios
for scenario in scenarios:
    cwd = os.getcwd()
    parent_dir = os.path.dirname(cwd)
    file_path = os.path.join(parent_dir, "results_march", scenario, "networks", "base_s_39___2050.nc")
    n = pypsa.Network(file_path)
    timestep = n.snapshot_weightings.iloc[0, 0]
    
    loads = n.loads
    
    # Check commodities
    for com in commodities:
        com_load = loads[loads.index.str.contains(com)]
        df.loc[com,scenario] = len(com_load)
    


# %%
years = [2030,2040,2050]
plot_steel_scenarios(scenarios)
plot_steel_scenarios_pathways(years,scenarios)


# %%
# COSTS AND EMISSIONS
years = [2030,2040,2050]
# Initialize the data_dict to store the data for both scenarios and years
data_dict = {
    "Annual system cost [b€/yr]": {"base_eu_regain": {}, "policy_eu_regain": {}},
    "CO2 emissions [MtCO2/yr]": {"base_eu_regain": {}, "policy_eu_regain": {}},
}

# Loop through each year and scenario to fill the data_dict
for year in years:
    for scenario in scenarios:
        cwd = os.getcwd()
        parent_dir = os.path.dirname(cwd)
        file_path = os.path.join(parent_dir, "results_march", scenario, "networks", f"base_s_39___{year}.nc")
        n = pypsa.Network(file_path)
        timestep = n.snapshot_weightings.iloc[0, 0]
        
        # Store the data for each metric
        data_dict["Annual system cost [b€/yr]"][scenario][year] = n.objective / 1e9
        data_dict["CO2 emissions [MtCO2/yr]"][scenario][year] = n.stores.loc['co2 atmosphere','e_nom_opt'] / 1e6

# Mapping of internal scenario names to display labels
scenario_labels = {
    "base_eu_regain": "Baseline",
    "policy_eu_regain": "Policy"
}

# Create the figure with two subplots: one for cost, one for emissions
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6), sharex=True)

# Define line styles and colors
line_styles = {"base_eu_regain": "-", "policy_eu_regain": "-"}  # Solid for baseline, dashed for policy
colors = {"base_eu_regain": "grey", "policy_eu_regain": "green"}

# Plot Annual System Cost on ax1
label_cost = "Annual system cost [b€/yr]"
for scenario in scenarios:
    ax1.plot(
        years,
        [data_dict[label_cost][scenario][year] for year in years],
        linestyle=line_styles[scenario],
        marker='o',
        label=scenario_labels[scenario],
        color=colors[scenario]
    )
ax1.set_ylabel(label_cost)
ax1.grid(True, linestyle = '-')
ax1.legend()
ax1.set_title("Annual System Cost")
ax1.set_xticks(years)

# Plot CO2 Emissions on ax2
label_emissions = "CO2 emissions [MtCO2/yr]"
for scenario in scenarios:
    ax2.plot(
        years,
        [data_dict[label_emissions][scenario][year] for year in years],
        linestyle=line_styles[scenario],
        marker='o',
        label=scenario_labels[scenario],
        color=colors[scenario]
    )
ax2.set_ylabel(label_emissions)
ax2.grid(True, linestyle = '-')
ax2.legend()
ax2.set_title("CO2 Emissions")
ax2.set_xticks(years)

# Improve layout and save
plt.tight_layout()
plt.savefig("./graphs/costs_emissions.png")
plt.show()

# %%


# Create the figure and the first axis
fig, ax1 = plt.subplots(figsize=(8, 5))

# Define line styles and colors
line_styles = {"base_eu_regain": "-", "policy_eu_regain": "--"}  # Solid for baseline, dashed for policy
colors = ["b", "r"]  # Define colors for different variables

# Plot Annual System Cost on the first y-axis (ax1)
for (label, values), color in zip(data_dict.items(), colors):
    if label == "Annual system cost [b€/yr]":
        for scenario in scenarios:
            ax1.plot(years, [values[scenario][year] for year in years], 
                     linestyle=line_styles[scenario], marker='o', label=f"{label} ({scenario})", color=color)

# Set labels and title for the first y-axis
ax1.set_ylabel("Annual system cost [b€/yr]", color="b")
ax1.tick_params(axis="y", labelcolor="b")

# Create the second y-axis (for CO2 emissions)
ax2 = ax1.twinx()

# Plot CO2 emissions on the second y-axis (ax2)
for (label, values), color in zip(data_dict.items(), colors):
    if label == "CO2 emissions [MtCO2/yr]":
        for scenario in scenarios:
            ax2.plot(years, [values[scenario][year] for year in years], 
                     linestyle=line_styles[scenario], marker='o', label=f"{label} ({scenario})", color=color)

# Set labels and title for the second y-axis
ax2.set_ylabel("CO2 emissions [MtCO2/yr]", color="r")
ax2.tick_params(axis="y", labelcolor="r")

# Set the title of the plot
#ax1.set_title("Annual System Cost and CO2 Emissions")

# Display the legend
fig.tight_layout()
fig.legend(loc="upper right",)# bbox_to_anchor=(1.1, 0.9))
plt.grid()
plt.show()
plt.savefig("./graphs/costs_emissions.png")

# %%
scenario = 'policy_eu_regain'
year = 2030
cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)
file_path = os.path.join(parent_dir, "results_march", scenario, "networks", f"base_s_39___{year}.nc")
n = pypsa.Network(file_path)
const = n.global_constraints
co2_links = n.links[n.links.bus0 == 'co2 atmosphere']

co2_links1 = n.links[n.links.bus1 == 'co2 atmosphere']

co2_links2 = n.links[n.links.bus2 == 'co2 atmosphere']

co2_links3 = n.links[n.links.bus3 == 'co2 atmosphere']



alinks1 = -n.links_t.p1[n.links_t.p1.columns.intersection(co2_links1.index)]

alinks2 = -n.links_t.p2[n.links_t.p2.columns.intersection(co2_links2.index)]

alinks3 = -n.links_t.p3[n.links_t.p3.columns.intersection(co2_links3.index)]

alinks1 = alinks1.sum()

alinks2 = alinks2.sum()

alinks3 = alinks3.sum()

tot_co2 = (alinks1.sum() + alinks2.sum() + alinks3.sum())/1e6

    
    