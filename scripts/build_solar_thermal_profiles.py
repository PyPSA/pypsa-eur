"""Build solar thermal collector time series."""

import geopandas as gpd
import atlite
import pandas as pd
import xarray as xr
import numpy as np

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_solar_thermal_profiles',
            weather_year='',
            simpl='',
            clusters=48,
        )

    if 'snakemake' not in globals():
        from vresutils import Dict
        import yaml
        snakemake = Dict()
        with open('config.yaml') as f:
            snakemake.config = yaml.safe_load(f)
        snakemake.input = Dict()
        snakemake.output = Dict()

    config = snakemake.config['solar_thermal']
    
    year = snakemake.wildcards.weather_year
    snapshots = dict(start=year, end=str(int(year)+1), closed="left") if year else snakemake.config['snapshots']
    time = pd.date_range(freq='m', **snapshots)

    cutout_config = snakemake.config['atlite']['cutout']
    if year: cutout_name = cutout_config.format(weather_year=year)
    cutout = atlite.Cutout(cutout_config).sel(time=time)

    clustered_regions = gpd.read_file(
        snakemake.input.regions_onshore).set_index('name').buffer(0).squeeze()

    I = cutout.indicatormatrix(clustered_regions)

    for area in ["total", "rural", "urban"]:

        pop_layout = xr.open_dataarray(snakemake.input[f'pop_layout_{area}'])

        stacked_pop = pop_layout.stack(spatial=('y', 'x'))
        M = I.T.dot(np.diag(I.dot(stacked_pop)))

        nonzero_sum = M.sum(axis=0, keepdims=True)
        nonzero_sum[nonzero_sum == 0.] = 1.
        M_tilde = M / nonzero_sum

        solar_thermal = cutout.solar_thermal(**config, matrix=M_tilde.T,
                                             index=clustered_regions.index)

        solar_thermal.to_netcdf(snakemake.output[f"solar_thermal_{area}"])
