
import geopandas as gpd
import atlite
import pandas as pd
import xarray as xr
import scipy as sp
import helper

if 'snakemake' not in globals():
    from vresutils import Dict
    import yaml
    snakemake = Dict()
    with open('config.yaml') as f:
        snakemake.config = yaml.safe_load(f)
    snakemake.input = Dict()
    snakemake.output = Dict()

time = pd.date_range(freq='m', **snakemake.config['snapshots'])
params = dict(years=slice(*time.year[[0, -1]]), months=slice(*time.month[[0, -1]]))



cutout = atlite.Cutout(snakemake.config['atlite']['cutout_name'],
                       cutout_dir=snakemake.config['atlite']['cutout_dir'],
                       **params)

clustered_busregions_as_geopd = gpd.read_file(snakemake.input.regions_onshore).set_index('name', drop=True)

clustered_busregions = pd.Series(clustered_busregions_as_geopd.geometry, index=clustered_busregions_as_geopd.index)

helper.clean_invalid_geometries(clustered_busregions)

I = cutout.indicatormatrix(clustered_busregions)


for item in ["total","rural","urban"]:

    pop_layout = xr.open_dataarray(snakemake.input['pop_layout_'+item])

    M = I.T.dot(sp.diag(I.dot(pop_layout.stack(spatial=('y', 'x')))))
    nonzero_sum = M.sum(axis=0, keepdims=True)
    nonzero_sum[nonzero_sum == 0.] = 1.
    M_tilde = M/nonzero_sum

    solar_thermal_angle = 45.
    #should clearsky_model be "simple" or "enhanced"?
    solar_thermal = cutout.solar_thermal(clearsky_model="simple",
                                         orientation={'slope': solar_thermal_angle, 'azimuth': 180.},
                                         matrix = M_tilde.T,
                                         index=clustered_busregions.index)

    solar_thermal.to_netcdf(snakemake.output["solar_thermal_"+item])
