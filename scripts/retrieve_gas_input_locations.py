"""
Retrieve LNG gas infrastructure data from "https://globalenergymonitor.org/wp-content/uploads/2022/09/Europe-Gas-Tracker-August-2022.xlsx"
Note: This is the script that removes non EU grid countries and terminals that are not closely connected to the europe network
"""

import logging
import pandas as pd
import geopandas as gpd

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('retrieve_gas_input_locations')
        rootpath = '..'
    else:
        rootpath = '.'

    # LNG terminals

    lng_url="https://globalenergymonitor.org/wp-content/uploads/2022/09/Europe-Gas-Tracker-August-2022.xlsx"
    df = pd.read_excel(lng_url,storage_options={'User-Agent': 'Mozilla/5.0'}, sheet_name = 'LNG terminals - data')
    df = df.set_index("ComboID")

    remove_status = ['Cancelled']
    remove_country = ['Cyprus','Turkey']
    remove_terminal = ['Puerto de la Luz LNG Terminal','Gran Canaria LNG Terminal']

    for index in df.index:
        if df.Status[index] in remove_status:
            df = df.drop([index])
        elif df.Country[index] in remove_country:
            df = df.drop([index])
        elif df.TerminalName[index] in remove_terminal:
            df = df.drop([index])
        elif df.CapacityInMtpa[index] == "--":
            df = df.drop([index])
        else:
            continue
    
    geometry = gpd.points_from_xy(df['Longitude'], df['Latitude'])
    lng = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

    lng.to_file(snakemake.output[0])