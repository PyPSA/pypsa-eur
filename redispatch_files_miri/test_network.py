import pandas as pd
import pypsa 
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
import atlite

plt.style.use("bmh")
#matplotlib inline

n = pypsa.Network("results/networks/elec_s_37_ec_lv1.0_Co2L-3H.nc")
                
#/results/test-elec/networks/elec_s_6_ec_lcopt_Co2L-3H.nc")
n.plot();
#"results/networks/elec_s_6_ec_lcopt_Co2L-24H.nc"

#temporal resolution
n.snapshots[:10]
len(n.snapshots)

n.lines.s_max_pu = 0.7
n.lines.loc[["316", "527", "602"], "s_nom"] = 1715

#static component data
n.lines.head()
n.snapshots()
n.generators.head()
n.storage_units.head()
n.loads.head()
n.loads_t.p_set.head()


#carriers
n.generators.groupby("carrier").p_nom.sum().div(1e3).plot.barh()
plt.xlabel('GW')

#plot
fig = plt.figure()
ax = plt.axes(projection=ccrs.EqualEarth())

n.plot(
    ax=ax,
    bus_sizes=load / 2e5,
);

#optimize
n.lopf(solver_name='cbc')