
# Build mapping between grid cells and population (total, urban, rural)

import atlite
import pandas as pd
import xarray as xr

from vresutils import shapes as vshapes

import geopandas as gpd


if 'snakemake' not in globals():
    from vresutils import Dict
    import yaml
    snakemake = Dict()
    with open('config.yaml') as f:
        snakemake.config = yaml.safe_load(f)
    snakemake.input = Dict()
    snakemake.output = Dict()

    snakemake.input["urban_percent"] = "data/urban_percent.csv"

cutout = atlite.Cutout(snakemake.config['atlite']['cutout_name'],
                       cutout_dir=snakemake.config['atlite']['cutout_dir'])

grid_cells = cutout.grid_cells()

#nuts3 has columns country, gdp, pop, geometry
#population is given in dimensions of 1e3=k
nuts3 = gpd.read_file(snakemake.input.nuts3_shapes).set_index('index')


# Indicator matrix NUTS3 -> grid cells
I = atlite.cutout.compute_indicatormatrix(nuts3.geometry, grid_cells)

# Indicator matrix grid_cells -> NUTS3; inprinciple Iinv*I is identity
# but imprecisions mean not perfect
Iinv = cutout.indicatormatrix(nuts3.geometry)

countries = nuts3.country.value_counts().index.sort_values()

urban_fraction = pd.read_csv(snakemake.input.urban_percent,
                             header=None,index_col=0,squeeze=True)/100.

#fill missing Balkans values
missing = ["AL","ME","MK"]
reference = ["RS","BA"]
urban_fraction = urban_fraction.reindex(urban_fraction.index.union(missing))
urban_fraction.loc[missing] = urban_fraction[reference].mean()


#population in each grid cell
pop_cells = pd.Series(I.dot(nuts3['pop']))

#in km^2
cell_areas = pd.Series(cutout.grid_cells()).map(vshapes.area)/1e6

#pop per km^2
density_cells = pop_cells/cell_areas


#rural or urban population in grid cell
pop_rural = pd.Series(0.,density_cells.index)
pop_urban = pd.Series(0.,density_cells.index)

for ct in countries:
    print(ct,urban_fraction[ct])

    indicator_nuts3_ct = pd.Series(0.,nuts3.index)
    indicator_nuts3_ct[nuts3.index[nuts3.country==ct]] = 1.

    indicator_cells_ct = pd.Series(Iinv.T.dot(indicator_nuts3_ct))

    density_cells_ct = indicator_cells_ct*density_cells

    pop_cells_ct = indicator_cells_ct*pop_cells

    #correct for imprecision of Iinv*I
    pop_ct = nuts3['pop'][indicator_nuts3_ct.index[indicator_nuts3_ct == 1.]].sum()
    pop_cells_ct = pop_cells_ct*pop_ct/pop_cells_ct.sum()

    # The first low density grid cells to reach rural fraction are rural
    index_from_low_d_to_high_d = density_cells_ct.sort_values().index
    pop_ct_rural_b = pop_cells_ct[index_from_low_d_to_high_d].cumsum()/pop_cells_ct.sum() < (1-urban_fraction[ct])
    pop_ct_urban_b = ~pop_ct_rural_b

    pop_ct_rural_b[indicator_cells_ct==0.] = False
    pop_ct_urban_b[indicator_cells_ct==0.] = False

    pop_rural += pop_cells_ct.where(pop_ct_rural_b,0.)
    pop_urban += pop_cells_ct.where(pop_ct_urban_b,0.)

pop_cells = {"total" : pop_cells}

pop_cells["rural"] = pop_rural
pop_cells["urban"] = pop_urban

for key in pop_cells.keys():
    layout = xr.DataArray(pop_cells[key].values.reshape(cutout.shape),
                          [('y', cutout.coords['y']), ('x', cutout.coords['x'])])

    layout.to_netcdf(snakemake.output["pop_layout_"+key])
