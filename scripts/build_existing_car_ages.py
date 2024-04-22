# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Read in existing car and truck ages from ACEA report (2024)
"""


import tabula
import numpy as np
import country_converter as coco
import pandas as pd
import logging
logger = logging.getLogger(__name__)
from _helpers import configure_logging, set_scenario_config

def distirbute_cars(row):   
    if row.name=="EU": return row
    i = 0
    for j, number in enumerate(row): 
        if np.isnan(number):
            i += 1
        else:
            distributed = number / (i+1)
            row.iloc[j-i:j+1] = distributed
            i = 0
    return row

def distribute_older_ages(ds, number_of_years=8):
    "Distribute columne <10 years based on averaged data"
    share_last_year = ds.fillna(0).loc['>10 years']
    ds.drop('>10 years', inplace=True)
    mean = ds.fillna(0).mean()
    # convert index to int
    ds.index = ds.index.astype(int)
    years_i = np.arange(2012, 2012-number_of_years, -1)
    av = share_last_year/mean
    distributed = pd.Series(data=0, index=years_i, dtype="float64")
    for i, year in enumerate(years_i): 
        if i+1 < av:
            distributed[year] = mean
        elif (i == int(av)):
            distributed[year] = (av - i)  * mean
        else:
            distributed[year] = 0
        
    if share_last_year - distributed.iloc[:-1].sum()> mean:
        distributed.iloc[-1] = share_last_year - distributed.iloc[:-1].sum()
        
    return pd.concat([ds, distributed])

def read_ages_from_pdf(fn, page, name="CARS BY AGE"):
    tables = tabula.read_pdf(fn, pages=page)
    table = tables[0]
    table.set_index(name, inplace=True)

    table.columns = table.iloc[1,:]
    table = table.iloc[2:,:]

    table = table.applymap(lambda x: np.nan if isinstance(x, str) and '–' in x else x)
    table = table.apply(lambda x: (x.str.replace(",","")).astype(float))
    
    # convert into to iso-code
    cc = coco.CountryConverter()
    table.index = cc.convert(table.index, to="iso2")
    table.drop('not found', inplace=True)
    
    table = table.apply(lambda x: distirbute_cars(x), axis=1)
    
    # get shares instead of total numbers
    table.iloc[:, :-2] = table.iloc[:, :-2].div(table["Total"], axis=0)
    table.rename(columns={'(in years)': "Average age"}, inplace=True)
    df = table.iloc[:,:-2]
    # swedish data is missing, filling with average data
    df = df.fillna(df.mean())
    final_table = df.apply(distribute_older_ages, axis=1)
    
    return pd.concat([final_table, table.iloc[:, -2:]], axis=1)


#%%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake(
            "build_existing_car_ages",
        )
    
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    fn = snakemake.input.ACEA_report
    # car ages
    car_ages = read_ages_from_pdf(fn, page=11, name="CARS BY AGE")
    
    # trucks ages
    truck_ages = read_ages_from_pdf(fn, page=13, name="TRUCKS BY AGE")
    
    
    car_ages.to_csv(snakemake.output.car_ages)
    truck_ages.to_csv(snakemake.output.truck_ages)
