#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 09:36:58 2024

@author: lisa
"""

# electrification trend
from build_transport_demand import transport_cols
import pandas as pd
idees = pd.read_csv("/home/lisa/Documents/playground/pypsa-eur/resources/test-new-transport/energy_totals.csv",
                    index_col=[0,1])

for transport_type in ["light", "heavy"]:
    electrified = idees.reindex(columns=[f"electricity {light_col}" for light_col in transport_cols[transport_type]]).sum(axis=1)
    total = idees.reindex(columns=[f"total {light_col}" for light_col in transport_cols[transport_type]]).sum(axis=1)
    to_plot = (electrified/total).unstack().T
    to_plot.plot(title=f"{transport_type} transport share of electrified energy")
    plt.legend(bbox_to_anchor=(1,1))

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the logistic function
def logistic_function(t, L, k, t0):
    L = min(1, max_increase_per_year*25)
    return L / (1 + np.exp(-k * (t - t0)))

for transport_type in ["light", "heavy"]:
    
    electrified = idees.reindex(columns=[f"electricity {light_col}" for light_col in transport_cols[transport_type]]).sum(axis=1)
    total = idees.reindex(columns=[f"total {light_col}" for light_col in transport_cols[transport_type]]).sum(axis=1)
    to_plot = (electrified/total).unstack().T
    
    for country in to_plot.columns:
    
        # Get the historical data for the chosen country
        historical_data = to_plot[country].astype(float)
        
        # Convert the index and values to numpy arrays
        years = historical_data.index.to_numpy().astype(float)
        values = historical_data.to_numpy().astype(float)
        
        max_increase_per_year = car_reg.loc[car_reg.index.get_level_values(1).str.contains(country)].loc[transport_type][0]
        
        # Provide initial guesses for L, k, and t0
#         L: Start with a value that represents the expected maximum share.
#           If the current share is around 2.9% and expected to grow to near 100%,
#           a reasonable guess could be between 0.5 and 1.
#         k: A small positive number. A reasonable starting point is between 0.01 and 0.1. Increasing 
#              will make the curve steeper.
# 
#         t0 : A year in the middle of your dataset, around 2010-2020.
        initial_guesses = [0.03, 0.05, 2035]
        
        # Set bounds (lower, upper) for the parameters to constrain the fit
        if transport_type=="light" and country=="DE":
            bounds = ([0.02, 0.001, 2030], [1, 0.59, 2050])
            
        elif transport_type=="light":
            initial_guesses = [0.2, 0.2, 2035]
            bounds = ([0.02, 0.001, 2035], [1, 0.95, 2050])
        else:
            initial_guesses = [0.03, 0.2, 2045]
            bounds = ([0.02, 0.001, 2045], [1, 8, 2050])
        
        # Fit the logistic function to the historical data
        popt, pcov = curve_fit(logistic_function, years, values, p0=initial_guesses, bounds=bounds, maxfev=10000)
        
        # Extract the parameters
        L, k, t0 = popt
        
        # Project the future values up to 2050
        future_years = np.arange(2000, 2051).astype(float)
        projected_values = logistic_function(future_years, L, k, t0)
        
        # Initialize the restricted projected values with the first value from the fitted logistic function
        restricted_values = [projected_values[0]]
        
        # Apply the restriction iteratively
        for i in range(1, len(projected_values)):
            next_value = restricted_values[-1] + min(max_increase_per_year, projected_values[i] - restricted_values[-1])
            restricted_values.append(next_value)
        
        # Convert the restricted_values list to a numpy array
        restricted_values = np.array(restricted_values)
        
        # all new cars electrified
        all_elec = [min(max_increase_per_year*t, 1) for t in range(1, 27)]
        
        # Plot the historical and projected data
        plt.figure(figsize=(10, 6))
        plt.plot(years, values, 'o', label='Historical Data')
        plt.plot(future_years, projected_values, '-', label='Projected S-curve')
        plt.plot(future_years[25:], all_elec, '-', label='All new cars electric')
        if country=="DE" and transport_type=="light":
            plt.scatter(y=0.029, x=2024, color='r', label='Actual 2024 Value (2.9%)')
        if country=="AT" and transport_type=="light":
            plt.scatter(y=0.021, x=2022, color='r', label='Actual 2022 Value (2.1%)')
        if country=="BE" and transport_type=="light":
            plt.scatter(y=0.015, x=2022, color='r', label='Actual 2022 Value (1.5%)')
        if country=="FR" and transport_type=="light":
            plt.scatter(y=0.015, x=2022, color='r', label='Actual 2022 Value (1.5%)')
        if country=="IT" and transport_type=="light":
            plt.scatter(y=0.004, x=2022, color='r', label='Actual 2022 Value 0.4%)')
        if country=="ES" and transport_type=="light":
            plt.scatter(y=0.004, x=2022, color='r', label='Actual 2022 Value 0.4%)')
        if country=="PL" and transport_type=="light":
            plt.scatter(y=0.002, x=2022, color='r', label='Actual 2022 Value 0.2%)')
        plt.plot(future_years, restricted_values, '--', label=f'Restricted S-curve (max {round(max_increase_per_year*100)}% increase/year)')
        plt.xlabel('Year')
        plt.ylabel('Share of Electrified Transport')
        plt.title(f'Share of Electrified Transport for {country} {transport_type} (Historical and Projected)')
        plt.legend()
        plt.grid(True)
        plt.show()
        
        # Print the projected value for 2024
        print(f"-----{country}---")
        print(f"{transport_type}: {pd.Series(projected_values, index=future_years).loc[2022]}")
        print("******")
        #%%
# Polynomial fit
from numpy.polynomial.polynomial import Polynomial

# Choose a country to fit the logistic function
country = 'NL'

# Get the historical data for the chosen country
historical_data = to_plot[country].astype(float)

# Convert the index and values to numpy arrays
years = historical_data.index.astype(float).to_numpy()
values = historical_data.to_numpy()


# Fit a polynomial of degree 3
p = Polynomial.fit(years, values, 3)

# Project the future values up to 2050
future_years = np.arange(2000, 2051)
projected_values = p(future_years)

plt.figure(figsize=(10, 6))
plt.plot(years, values, 'o', label='Historical Data')
plt.plot(future_years, projected_values, '-', label='Polynomial Fit')
plt.xlabel('Year')
plt.ylabel('Share of Electrified Transport')
plt.title(f'Share of Electrified Transport for {country} (Historical and Projected)')
plt.legend()
plt.grid(True)
plt.show()
#%%
from scipy.optimize import curve_fit

def exponential_growth(t, a, b, c):
    return a * np.exp(b * (t - c))

initial_guesses = [0.01, 0.1, 2000]
popt, pcov = curve_fit(exponential_growth, years, values, p0=initial_guesses)

a, b, c = popt
projected_values = exponential_growth(future_years, a, b, c)

plt.figure(figsize=(10, 6))
plt.plot(years, values, 'o', label='Historical Data')
plt.plot(future_years, projected_values, '-', label='Exponential Growth')
plt.xlabel('Year')
plt.ylabel('Share of Electrified Transport')
plt.title(f'Share of Electrified Transport for {country} (Historical and Projected)')
plt.legend()
plt.grid(True)
plt.show()



