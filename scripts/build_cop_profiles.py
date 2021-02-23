
import xarray as xr

#quadratic regression based on Staffell et al. (2012)
#https://doi.org/10.1039/C2EE22653G

# COP is function of temp difference source to sink

cop_f = {"air" : lambda d_t: 6.81 -0.121*d_t + 0.000630*d_t**2,
         "soil" : lambda d_t: 8.77 -0.150*d_t + 0.000734*d_t**2}


for area in ["total", "urban", "rural"]:
    for source in ["air", "soil"]:

        source_T = xr.open_dataarray(snakemake.input["temp_{}_{}".format(source,area)])

        delta_T = snakemake.config['sector']['heat_pump_sink_T'] - source_T

        cop = cop_f[source](delta_T)

        cop.to_netcdf(snakemake.output["cop_{}_{}".format(source,area)])
