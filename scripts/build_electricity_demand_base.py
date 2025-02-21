# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Builds the electricity demand for base regions based on population and GDP.
"""

import logging
from itertools import product

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import scipy.sparse as sparse
import xarray as xr
from _helpers import configure_logging, set_scenario_config
from shapely.prepared import prep

logger = logging.getLogger(__name__)


def normed(s: pd.Series) -> pd.Series:
    return s / s.sum()


def shapes_to_shapes(orig: gpd.GeoSeries, dest: gpd.GeoSeries) -> sparse.lil_matrix:
    """
    Adopted from vresutils.transfer.Shapes2Shapes()
    """
    orig_prepped = list(map(prep, orig))
    transfer = sparse.lil_matrix((len(dest), len(orig)), dtype=float)

    for i, j in product(range(len(dest)), range(len(orig))):
        if orig_prepped[j].intersects(dest.iloc[i]):
            area = orig.iloc[j].intersection(dest.iloc[i]).area
            transfer[i, j] = area / dest.iloc[i].area

    return transfer


def upsample_load(
    n: pypsa.Network,
    regions_fn: str,
    load_fn: str,
    nuts3_fn: str,
    distribution_key: dict[str, float],
) -> pd.DataFrame:
    substation_lv_i = n.buses.index[n.buses["substation_lv"]]
    gdf_regions = gpd.read_file(regions_fn).set_index("name").reindex(substation_lv_i)
    load = pd.read_csv(load_fn, index_col=0, parse_dates=True)

    nuts3 = gpd.read_file(nuts3_fn).set_index("index")

    gdp_weight = distribution_key.get("gdp", 0.6)
    pop_weight = distribution_key.get("pop", 0.4)

    data_arrays = []

    for cntry, group in gdf_regions.geometry.groupby(gdf_regions.country):
        load_ct = load[cntry]

        if len(group) == 1:
            factors = pd.Series(1.0, index=group.index)

        else:
            nuts3_cntry = nuts3.loc[nuts3.country == cntry]
            transfer = shapes_to_shapes(group, nuts3_cntry.geometry).T.tocsr()
            gdp_n = pd.Series(
                transfer.dot(nuts3_cntry["gdp"].fillna(1.0).values), index=group.index
            )
            pop_n = pd.Series(
                transfer.dot(nuts3_cntry["pop"].fillna(1.0).values), index=group.index
            )

            factors = normed(gdp_weight * normed(gdp_n) + pop_weight * normed(pop_n))

        data_arrays.append(
            xr.DataArray(
                factors.values * load_ct.values[:, np.newaxis],
                dims=["time", "bus"],
                coords={"time": load_ct.index.values, "bus": factors.index.values},
            )
        )


    return xr.concat(data_arrays, dim="bus")



######################################## PyPSA-Spain
#
# Function to attach loads according to PyPSA-Spain methodology    ### EN PROCESO
#
def upsample_load_vPyPSA_Spain(
    n: pypsa.Network,
    regions_fn: str,
    load_fn: str,
    nuts3_fn: str,
    distribution_key: dict[str, float],
) -> pd.DataFrame:
    substation_lv_i = n.buses.index[n.buses["substation_lv"]]
    gdf_regions = gpd.read_file(regions_fn).set_index("name").reindex(substation_lv_i)
    load = pd.read_csv(load_fn, index_col=0, parse_dates=True)

    nuts3 = gpd.read_file(nuts3_fn).set_index("index")

    gdp_weight = distribution_key.get("gdp", 0.6)
    pop_weight = distribution_key.get("pop", 0.4)

    data_arrays = []

    # In this loop, "cntry" takes country codes 'ES', ..., and "group" is a geoseries with the geometries of the onshore regions (Voronoi cells)
    for cntry, group in gdf_regions.geometry.groupby(gdf_regions.country):

        # This is to do nothing if there is only one region within the country
        if len(group) == 1:
            factors = pd.Series(1.0, index=group.index)

        else:

            ##### Define a list where data_arrays for specific NUTSx regions will be obtained
            data_arrays_list = []


            for rr_ID in load.columns:  ##### rr_ID is the NUTSx ID

                load_local = load[rr_ID]

                ##### Take NUTS3 regions in NUTSx 
                nuts3_rr = nuts3[nuts3.index.str.contains(rr_ID)]

                ##### Get a matrix with area overlap percentages
                #  rows are ~522 substations, columns are NUTS3 in región 
                # .tocsr() is to put the sparse matrix in 'Compressed Sparse Row' format.            
                transfer = shapes_to_shapes(group, nuts3_rr.geometry).T.tocsr() 

                ##### Series with substations, (overlapped area) x (gdp) from corresponding NUTS3 in rr_ID
                gdp_n = pd.Series(
                    transfer.dot(nuts3_rr["gdp"].fillna(1.0).values), index=group.index
                        )

                ##### Series with substations, (overlapped area) x (pop) from corresponding NUTS3 in rr_ID
                pop_n = pd.Series(
                    transfer.dot(nuts3_rr["pop"].fillna(1.0).values), index=group.index
                        )

                factors = normed(gdp_weight * normed(gdp_n) + pop_weight * normed(pop_n))


                ##### Add data_array to data_arrays_list
                data_arrays_list.append(
                    xr.DataArray(
                        factors.values * load_local.values[:, np.newaxis],
                        dims=["time", "bus"],
                        coords={"time": load_local.index.values, "bus": factors.index.values},
                    )
                )

            ##### Sum up all the data_arrays in the list
            data_arrays_total =  xr.concat(data_arrays_list, dim="new_dim").sum(dim="new_dim")


        ##### Finally, append data_arrays_total to data_arrays (the list for countries in the original script)
        data_arrays.append(data_arrays_total)


    return xr.concat(data_arrays, dim="bus")

########################################



if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_electricity_demand_base")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params

    n = pypsa.Network(snakemake.input.base_network)



    ################################################## PyPSA-Spain
    #
    # Attach electricity demand according to PyPSA-Spain customisation
    #

    electricity_demand = snakemake.params.electricity_demand


    if electricity_demand['enable']:

        logger.info(f'##### [PyPSA-Spain] <build_electricity_demand_base>: Using upsample_load_vPyPSA_Spain to add electricity demand in ES.')

        load = upsample_load_vPyPSA_Spain(
            n,
            regions_fn=snakemake.input.regions,
            load_fn=snakemake.input.load,
            nuts3_fn=snakemake.input.nuts3,
            distribution_key=params.distribution_key,
        )

    else:
        load = upsample_load(
            n,
            regions_fn=snakemake.input.regions,
            load_fn=snakemake.input.load,
            nuts3_fn=snakemake.input.nuts3,
            distribution_key=params.distribution_key,
        )

    #
    #        
    ##################################################



    load.name = "electricity demand (MW)"
    comp = dict(zlib=True, complevel=9, least_significant_digit=5)
    load.to_netcdf(snakemake.output[0], encoding={load.name: comp})
