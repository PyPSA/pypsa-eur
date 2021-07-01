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

    time = pd.date_range(freq='h', **snakemake.config['snapshots'])
    cutout_config = snakemake.config['atlite']['cutout']
    cutout = atlite.Cutout(cutout_config).sel(time=time)

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
