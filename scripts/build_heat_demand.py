"""Build heat demand time series."""

import geopandas as gpd
import atlite
import pandas as pd
import xarray as xr
import numpy as np

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_heat_demands',
            weather_year='',
            simpl='',
            clusters=48,
        )

    cutout_name = snakemake.input.cutout
    year = snakemake.wildcards.weather_year
    drop_leap_day = snakemake.config["atlite"].get("drop_leap_day", False)

    if year:
        snapshots = dict(start=year, end=str(int(year)+1), closed="left")
        cutout_name = cutout_name.format(weather_year=year)
    else:
        snapshots = snakemake.config['snapshots']
    
    time = pd.date_range(freq='m', **snapshots)
    if drop_leap_day:
        time = time[~((time.month == 2) & (time.day == 29))]

    cutout = atlite.Cutout(cutout_name).sel(time=time)

    clustered_regions = gpd.read_file(
        snakemake.input.regions_onshore).set_index('name').buffer(0).squeeze()

    I = cutout.indicatormatrix(clustered_regions)

    for area in ["rural", "urban", "total"]:

        pop_layout = xr.open_dataarray(snakemake.input[f'pop_layout_{area}'])

        stacked_pop = pop_layout.stack(spatial=('y', 'x'))
        M = I.T.dot(np.diag(I.dot(stacked_pop)))

        heat_demand = cutout.heat_demand(
            matrix=M.T, index=clustered_regions.index)

        heat_demand.to_netcdf(snakemake.output[f"heat_demand_{area}"])
