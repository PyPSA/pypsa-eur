#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 06:50:42 2024

@author: lisa
"""
import pandas as pd
import matplotlib.pyplot as plt
from build_energy_totals import car_types
from build_transport_demand import transport_cols
import numpy as np
# heat pumps
wished_scenarios = ["base", "fast", "slow", "low-demand", "high-demand"]
fn = "/home/lisa/Documents/endogenous_transport/data/statistic_id1434431_total-stock-of-heat-pumps-in-europe-2005-2022.xlsx"
historical = pd.read_excel(fn, index_col=[1], sheet_name="Data")
pumps = (heat_df/typical_generation).reindex(columns=wished_scenarios, level=0).iloc[:2].sum().unstack().T

historical = historical.dropna(axis=1, how="all").dropna().iloc[:,0].rename("historical")
historical.index = historical.index.astype(int)
pumps.index = pumps.index.astype(int)

to_plot = pd.concat([historical, pumps])
to_plot.loc[2022] = to_plot.loc[2022].fillna(method="ffill")
fig, ax = plt.subplots()
to_plot.plot(ax=ax, color=["black", "green", "red", "blue", "pink"],
             style = ["-", "--", "--", "--", "--"])
ax.set_ylabel("Number of heat pumps [millions]")
ax.set_xlim([2005, 2050])

fig.savefig(snakemake.output.balances[:-19] + "heat_pump_installation.pdf",
            bbox_inches="tight")
#%%
new_index = pd.Index(range(to_plot.index[0], to_plot.index[-1]+1), name='year')
to_plot_reindexed = to_plot.reindex(new_index)

# Step 2: Interpolate missing values for all scenarios
to_plot_interpolated = to_plot_reindexed.interpolate(method='linear')

# Step 3: Calculate the yearly difference
yearly_diff = to_plot_interpolated.diff()
yearly_diff.loc[historical.reindex(new_index).isna(), "historical"] = np.nan

(yearly_diff.shift(-1)).plot()
plt.ylabel("Annual additional heat pumps \n [million heat pumps]")
plt.savefig(snakemake.output.balances[:-19] + "heat_pump_annual_additional.pdf",
            bbox_inches="tight")
#%%
# renewable capacities
fn = "/home/lisa/Documents/endogenous_transport/data/statistic_id864189_renewable-energy-capacity-in-europe-2010-2023.xlsx"
historical = pd.read_excel(fn, index_col=[1], sheet_name="Data")
historical = historical.dropna(axis=1, how="all").dropna().iloc[:,0].rename("historical")
historical.index = historical.index.astype(int)

caps = capacities.droplevel(0).loc[renewables].droplevel([1,2,3], axis=1).sum().unstack().T.rename(index=lambda x: int(x))[wished_scenarios]

to_plot = pd.concat([historical, caps])
to_plot.loc[2023] = to_plot.loc[2023].fillna(method="ffill")
fig, ax = plt.subplots()
(to_plot/1e3).plot(ax=ax, color=["black", "green", "red", "blue", "pink", "orange"],
             style = ["-", "--", "--", "--", "--", "--"])
ax.set_ylabel("Installed renewable capacity [GW]")
ax.set_xlim([2010, 2050])
fig.savefig(snakemake.output.balances[:-19] + "res_installation.pdf",
            bbox_inches="tight")
#%%
new_index = pd.Index(range(to_plot.index[0], to_plot.index[-1]+1), name='year')
to_plot_reindexed = to_plot.reindex(new_index)

# Step 2: Interpolate missing values for all scenarios
to_plot_interpolated = to_plot_reindexed.interpolate(method='linear')

# Step 3: Calculate the yearly difference
yearly_diff = to_plot_interpolated.diff()
yearly_diff.loc[historical.reindex(new_index).isna(), "historical"] = np.nan

(yearly_diff.shift(-1)/1e3).plot()
plt.ylabel("Annual additional RES capacity [GW]")
plt.savefig(snakemake.output.balances[:-19] + "res_annual_additional.pdf",
            bbox_inches="tight")
#%%
transport_demand = pd.read_csv("/home/lisa/Documents/playground/pypsa-eur/resources/test-after-merge-v2/transport_data.csv",
                            index_col=[0,1])
transport_types = ["light", "heavy"]
fn = "/home/lisa/Documents/endogenous_transport/graphics/demand/"
for transport_type in transport_types:
    fig, ax = plt.subplots()
    cols = [f"mio km-driven {car_type}" for car_type in transport_cols[transport_type]]
    demand = transport_demand.loc[:, cols]
    grouped = demand.groupby(level=1).sum().sum(axis=1)
    
    # convert million to trillion
    grouped /= 1e6
    grouped.plot(ax=ax, 
                 title=f"transport demand {transport_type} duty")
    
    # Get x and y values for the linear fit
    x = np.arange(len(grouped))
    y = grouped.values
    
    # Perform linear regression using polyfit
    slope, intercept = np.polyfit(x, y, 1)
    
    # Create the linear fit line based on the calculated slope and intercept
    linear_fit = slope * x + intercept
    
    # Plot the linear fit line
    ax.plot(grouped.index, linear_fit, color='gray',linestyle='--',
            label=f'Linear fit')
    
    starting_value = grouped.iloc[0]  # Starting value from the first year

    percentage_increase_per_year = (slope / starting_value) * 100
    
    ax.set_title(f"Transport demand {transport_type} duty \n Avgerage increase: {percentage_increase_per_year:.1f}%/year")
    
    ax.set_ylabel("trillion km driven")
    ax.legend()
    
    fig.savefig(fn + f"historic_demand_{transport_type}.pdf",
                bbox_inches="tight")
    
#%% car registration
reg = pd.read_csv("/home/lisa/Documents/playground/pypsa-eur/resources/test-after-merge-v2/car_registration_s_39.csv", index_col=[0,1])
for transport_type in transport_types:
    fig, ax = plt.subplots()
    to_plot = reg.loc[transport_type]
    to_plot = (to_plot.groupby(to_plot.index.str[:2]).mean().iloc[:,0].sort_values()*100)
    to_plot.name = "country"
    to_plot.plot(ax=ax, kind="bar", title=f"Car registration {transport_type} duty \n Average {to_plot.mean():.1f}%")
    ax.set_ylabel("New car registration \n compared to total fleet [%]")
    ax.set_xlabel("Country")
    fig.savefig(fn + f"new_car_registration_{transport_type}.pdf",
                bbox_inches="tight")
#%%
fn = "/home/lisa/Documents/endogenous_transport/graphics/demand/"
demand = pd.read_csv("/home/lisa/Documents/playground/pypsa-eur/resources/test-after-merge-v2/transport_demand_s_39.csv",
                     header=[0,1], index_col=0)


for transport_type in demand.columns.levels[0]:
    fig, ax = plt.subplots()
    to_plot = demand[transport_type].iloc[:7*24,:]/1e6
    to_plot.plot(legend=False, alpha=0.2, ax=ax)
    to_plot.mean(axis=1).plot(color="black", ax=ax, title=transport_type, label="Mean")

    # Calculate and plot the min/max contour area
    ax.fill_between(to_plot.index, to_plot.min(axis=1), to_plot.max(axis=1), color="gray", alpha=0.2)
    
    # Set y-axis label
    ax.set_ylabel("vehicle-km driven [million 100km]")
    
    # Rotate x-labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylim(bottom=0)
    
    # Add legend
    # Get existing handles and labels
    handles, labels = ax.get_legend_handles_labels()
    
    # Filter only the mean line (you could do this by checking labels or handles)
    mean_handle = [h for h, l in zip(handles, labels) if l == "Mean"]
    
    # Reapply the legend with only the mean line handle
    ax.legend(mean_handle, ["Mean"])
   

    fig.savefig(fn + f"example_week_demand_{transport_type}.pdf",
                bbox_inches="tight")
    
to_plot = demand.sum().unstack().T
(to_plot.groupby(to_plot.index.str[:2]).sum()/1e6).plot(kind="bar")
plt.ylabel("vehicle-km driven [million 100 km]")
plt.xlabel("")
plt.savefig(fn + f"total_transport_demand.pdf",
            bbox_inches="tight")
