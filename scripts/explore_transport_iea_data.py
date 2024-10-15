#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 11:06:06 2024
https://www.iea.org/data-and-statistics/data-product/global-ev-outlook-2024#global-ev-data
@author: lisa
"""
import pandas as pd
import country_converter as coco
cc = coco.CountryConverter()

iea_data = pd.read_csv("/home/lisa/Downloads/IEA Global EV Data 2024.csv", 
                       index_col=[0,1,2,3,4,5])

iea_data.rename(index = lambda x: cc.convert(x, to="ISO2"), level=0, inplace=True)

iea_data = iea_data[iea_data.index.get_level_values(0).isin(countries)]

#%%
category = "Historical"
parameter = "EV sales share"
mode = "Cars"
powertrain = "EV"

for mode in ['Buses', 'Cars', 'Trucks', 'Vans']:
    # Select the rows based on the criteria
    selected_data = iea_data.loc[
        (slice(None), category, parameter, mode, powertrain, slice(None)),
        :
    ]
    to_plot = selected_data.droplevel([1,2,3,4])["value"].drop(["not found"]).unstack().sort_values(by=2023, ascending=False).T  
    to_plot.plot(title=mode)
    plt.ylabel(parameter)
    plt.legend(bbox_to_anchor=(1,1))
    
    
data = (iea_data.xs(category, level="category").xs(parameter, level="parameter")
        .xs("EV", level="powertrain")["value"].drop(["not found"]).unstack().T)
data.to_csv("/home/lisa/Documents/endogenous_transport/data/iea_outlook2024_EVsalesshare_per_country.csv")
#%%
share = (iea_data.xs("EV stock share", level=2).xs(category, level="category")
        .xs("EV", level="powertrain").value)
stock = (iea_data.xs("EV stock", level=2).xs(category, level="category")
         .xs("BEV", level="powertrain").value)
