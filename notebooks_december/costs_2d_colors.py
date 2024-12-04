# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 17:19:10 2024

@author: alice
"""

import cartopy.crs as ccrs
import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import yaml
import numpy as np

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

annualized_df.index = annualized_df.index.str[:2]

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

df = pd.DataFrame(0, index = annualized_reg.index, columns = ['change','steel'])
df["change"] = df.index.map(annualized_change) 
df["steel"] = df.index.map(delta_steel)

# Extract data
x = df["steel"]
y = df["change"]

# Identify outliers
num_outliers = 1  # Change this to 1 outlier
outliers_x_positive = x.nlargest(num_outliers).index
outliers_x_negative = x.nsmallest(num_outliers).index
outliers_y_positive = y.nlargest(num_outliers).index
outliers_y_negative = y.nsmallest(num_outliers).index

outlier_indices = (
    outliers_x_positive.union(outliers_x_negative)
    .union(outliers_y_positive)
    .union(outliers_y_negative)
)

# Separate outliers from main data
df_outliers = df.loc[outlier_indices]
df_main = df.drop(outlier_indices)

# Create the figure
fig, ax = plt.subplots(figsize=(10, 8))

# Scatter plot for main data
scatter = ax.scatter(
    df['steel'],#_main["steel"],
    df['change'],#_main["change"],
    c="blue",
    alpha=0.4,
    edgecolor="black"
)


# Annotate all data points (country names)
for idx, row in df.iterrows():#_main.iterrows():
    ax.text(
        row["steel"]*1.05, 
        row["change"]*1.05, 
        f"{idx}",  # Show the country name
        ha="center", 
        va="center", 
        fontsize=13, 
        alpha=1
    )
    
reduction = 0.1

"""
# Annotate outliers on the borders (outside the plot area)
for idx, row in df_outliers.iterrows():
    if idx in outliers_x_positive:
        ax.text(
            row["steel"], 
            ax.get_ylim()[1] * reduction,  # Move the annotation above the plot
            f"{idx}: {row['steel']:.2f}, {row['change']:.2f}", 
            ha="center", 
            va="bottom", 
            color="red", 
            fontsize=10
        )
    elif idx in outliers_x_negative:
        ax.text(
            row["steel"], 
            ax.get_ylim()[0] * reduction,  # Move the annotation below the plot
            f"{idx}: {row['steel']:.2f}, {row['change']:.2f}", 
            ha="center", 
            va="top", 
            color="red", 
            fontsize=10
        )
    elif idx in outliers_y_positive:
        ax.text(
            ax.get_xlim()[1] * reduction,  # Move the annotation to the right of the plot
            row["change"], 
            f"{idx}: {row['steel']:.2f}, {row['change']:.2f}", 
            ha="left", 
            va="center", 
            color="green", 
            fontsize=10
        )
    elif idx in outliers_y_negative:
        ax.text(
            ax.get_xlim()[0] * reduction,  # Move the annotation to the left of the plot
            row["change"], 
            f"{idx}: {row['steel']:.2f}, {row['change']:.2f}", 
            ha="right", 
            va="center", 
            color="green", 
            fontsize=10
        )
"""
# Add axis labels
ax.set_xlabel("Steel Production (Mtons)", fontsize=12)
ax.set_ylabel("Change in Costs (€/person)", fontsize=12)
ax.set_title("Steel Production vs. Change in Costs (Outliers Highlighted)", fontsize=14)

# Add grid for readability
ax.grid(visible=True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

# Save and show the plot
plt.savefig("graphs/steel_vs_change_outliers.png", bbox_inches="tight")
plt.show()


# %%

# Irland check
pol_reg_nodal['Unnamed: 2'] = pol_reg_nodal['Unnamed: 2'].fillna("")
ie_pol_reg_nodal = pol_reg_nodal[pol_reg_nodal['Unnamed: 2'].str.contains('IE')]
#ie_pol_reg_nodal['scenario'] = 'pol_reg'

pol_eu_nodal['Unnamed: 2'] = pol_eu_nodal['Unnamed: 2'].fillna("")
ie_pol_eu_nodal = pol_eu_nodal[pol_eu_nodal['Unnamed: 2'].str.contains('IE')]
#ie_pol_eu_nodal['scenario'] = 'pol_eu'

ie_names = ['pol_eu','pol_reg']

ie_df = [ie_pol_eu_nodal,ie_pol_reg_nodal]

processed_csvs = []
columns_to_drop =  ["cluster", "Unnamed: 1"]

years_tech = ['Tech'] + years

for df, name in zip(ie_df, ie_names):
        
    df = df.drop(columns=columns_to_drop, errors='ignore')  # 'errors=ignore' avoids issues if some columns don't exist
    df = df.dropna(subset=["Unnamed: 2"])
    df = df.set_index("Unnamed: 2")
    #df = df[df.index.str.match(r"^[A-Z]")]
    #df = df[~df.index.str.startswith("EU")]
    #df = df[~df.index.str.startswith("H2")]
    df.columns = years_tech
    #df = df.apply(pd.to_numeric, errors='coerce')
    df["scenario"] = name
    df = df.set_index("scenario", append=True)
    processed_csvs.append(df) 

ie_combined_df = pd.concat(processed_csvs)

r = 0.04 # social discount rate
ie_combined_df[years] = ie_combined_df[years].apply(pd.to_numeric, errors='coerce')

for year in years:
    ie_combined_df[f'Annualized_{year}'] = ie_combined_df[year] / ((1 + r) ** (year - 2020))

ie_combined_df['value'] = ie_combined_df['Annualized_2030']  + ie_combined_df['Annualized_2040'] + ie_combined_df['Annualized_2050']
ie_combined_df = ie_combined_df.reset_index()

ie_combined_df = ie_combined_df[['scenario','Tech','value']]
ie_combined_df = ie_combined_df.groupby(['Tech', 'scenario'], as_index=False)['value'].sum()

ie_combined_df_wide = ie_combined_df.pivot(index='Tech', columns='scenario', values='value').reset_index()
ie_combined_df_wide.index = ie_combined_df_wide['Tech']
ie_combined_df_wide = ie_combined_df_wide.drop('Tech', axis=1)

ie_combined_df_wide['reg-eu'] = ie_combined_df_wide['pol_reg'] - ie_combined_df_wide['pol_eu']

filtered_df = ie_combined_df_wide[ie_combined_df_wide["reg-eu"].abs() >= 1]
