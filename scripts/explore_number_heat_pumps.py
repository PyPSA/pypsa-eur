#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 06:50:42 2024

@author: lisa
"""
import pandas as pd
import matplotlib.pyplot as plt

# heat pumps
fn = "/home/lisa/Documents/endogenous_transport/data/statistic_id1434431_total-stock-of-heat-pumps-in-europe-2005-2022.xlsx"
historical = pd.read_excel(fn, index_col=[1], sheet_name="Data")
pumps = (heat_df/typical_generation)[wished_scenarios].iloc[:2].sum().unstack().T

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
(to_plot/1e3).plot(ax=ax, color=["black", "green", "red", "blue", "pink"],
             style = ["-", "--", "--", "--", "--"])
ax.set_ylabel("Installed renewable capacity [GW]")
ax.set_xlim([2010, 2050])
#%%
