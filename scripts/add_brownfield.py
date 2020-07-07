# coding: utf-8

import logging
logger = logging.getLogger(__name__)
import pandas as pd
idx = pd.IndexSlice

import numpy as np
import scipy as sp
import xarray as xr
import re, os

from six import iteritems, string_types

import pypsa

import yaml

import pytz

from vresutils.costdata import annuity

from add_existing_baseyear import add_power_capacities_installed_before_baseyear 

from add_existing_baseyear import add_heating_capacities_installed_before_baseyear 

from add_existing_baseyear import prepare_costs

#First tell PyPSA that links can have multiple outputs by
#overriding the component_attrs. This can be done for
#as many buses as you need with format busi for i = 2,3,4,5,....
#See https://pypsa.org/doc/components.html#link-with-multiple-outputs-or-inputs

override_component_attrs = pypsa.descriptors.Dict({k : v.copy() for k,v in pypsa.components.component_attrs.items()})
override_component_attrs["Link"].loc["bus2"] = ["string",np.nan,np.nan,"2nd bus","Input (optional)"]
override_component_attrs["Link"].loc["bus3"] = ["string",np.nan,np.nan,"3rd bus","Input (optional)"]
override_component_attrs["Link"].loc["efficiency2"] = ["static or series","per unit",1.,"2nd bus efficiency","Input (optional)"]
override_component_attrs["Link"].loc["efficiency3"] = ["static or series","per unit",1.,"3rd bus efficiency","Input (optional)"]
override_component_attrs["Link"].loc["p2"] = ["series","MW",0.,"2nd bus output","Output"]
override_component_attrs["Link"].loc["p3"] = ["series","MW",0.,"3rd bus output","Output"]


def prepare_costs():

    #set all asset costs and other parameters
    #costs = pd.read_csv(snakemake.input.costs,index_col=list(range(3))).sort_index()
    
    costs = pd.read_csv(snakemake.input.costs,index_col=list(range(2))).sort_index()
    
    #correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"),"value"]*=1e3
    costs.loc[costs.unit.str.contains("USD"),"value"]*=snakemake.config['costs']['USD2013_to_EUR2013']

    #cost_year = snakemake.config['costs']['year']
    #costs = costs.loc[idx[:,cost_year,:],"value"].unstack(level=2).groupby(level="technology").sum(min_count=1)
    costs = costs.loc[:, "value"].unstack(level=1).groupby("technology").sum()
    costs = costs.fillna({"CO2 intensity" : 0,
                          "FOM" : 0,
                          "VOM" : 0,
                          "discount rate" : snakemake.config['costs']['discountrate'],
                          "efficiency" : 1,
                          "fuel" : 0,
                          "investment" : 0,
                          "lifetime" : 25
    })

    costs["fixed"] = [(annuity(v["lifetime"],v["discount rate"])+v["FOM"]/100.)*v["investment"]*Nyears for i,v in costs.iterrows()]
    return costs

    
def add_brownfield(n,n_p):
    print("adding brownfield")
    
    previous_timestep=snakemake.config['scenario']['planning_horizons'][snakemake.config['scenario']['planning_horizons'].index(year)-1]
    previous_timesteps=snakemake.config['scenario']['planning_horizons'][0:snakemake.config['scenario']['planning_horizons'].index(year)]

    
    #add generators from previous step     
    n_p.mremove("Generator", [index for index in n_p.generators.index.to_list() if '<'+snakemake.config['scenario']['planning_horizons'][0] in index])
    n_p.mremove("Generator", [index for index in n_p.generators.index.to_list() if 'ror' in index])
        
    n_p.generators.index=np.where(n_p.generators.index.str[-4:].isin(previous_timesteps)==False,
                                  n_p.generators.index + '-' + previous_timestep, 
                                  n_p.generators.index)

    n.madd("Generator", 
            n_p.generators.index,
            bus=n_p.generators.bus,
            carrier=n_p.generators.carrier,
            p_nom=n_p.generators.p_nom_opt,
            marginal_cost=n_p.generators.marginal_cost,
            capital_cost=n_p.generators.capital_cost,
            efficiency=n_p.generators.efficiency,
            p_max_pu=n_p.generators_t.p_max_pu)
            
    #add stores from previous steps
    n_p.mremove("Store", ['co2 atmosphere', 'co2 stored', 'EU gas Store'] )
    n_p.stores.index=np.where(n_p.stores.index.str[-4:].isin(previous_timesteps)==False,
                                  n_p.stores.index + '-' + previous_timestep, 
                                  n_p.stores.index)
    n.madd("Store", 
            n_p.stores.index,
            bus=n_p.stores.bus,
            carrier=n_p.stores.carrier,
            e_nom=n_p.stores.e_nom_opt,
            e_cyclic=True,
            capital_cost=n_p.stores.capital_cost)
    
    ## add links from previous steps          
    # TODO: add_chp_constraint() in solve_network needs to be adjusted
    n_p.mremove("Link", [index for index in n_p.links.index.to_list() if '<'+snakemake.config['scenario']['planning_horizons'][0] in index])
    n_p.mremove("Link", [index for index in n_p.links.index.to_list() if 'CHP' in index])
    
    n_p.links.index=np.where(n_p.links.index.str[-4:].isin(previous_timesteps)==False,
                                  n_p.links.index + '-' + previous_timestep, 
                                  n_p.links.index)
    n.madd("Link", 
            n_p.links.index,
            bus0=n_p.links.bus0,
            bus1=n_p.links.bus1,
            bus2=n_p.links.bus2,
            carrier=n_p.links.carrier,
            p_nom=n_p.links.p_nom_opt,
            marginal_cost=n_p.links.marginal_cost,
            capital_cost=n_p.links.capital_cost,
            efficiency=n_p.links.efficiency,
            efficiency2=n_p.links.efficiency2)
    
    
if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake
        snakemake = MockSnakemake(
            wildcards=dict(network='elec', simpl='', clusters='37', lv='1.0', 
                           sector_opts='Co2L0-168H-T-H-B-I-solar3-dist1',
                           co2_budget_name='go',
                           planning_horizons='2050'),
            input=dict(network='pypsa-eur-sec/results/test/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}__{sector_opts}_{co2_budget_name}_{planning_horizons}.nc', 
                       network_p='pypsa-eur-sec/results/test/postnetworks/{network}_s{simpl}_{clusters}_lv{lv}__{sector_opts}_{co2_budget_name}_2040.nc',                        
                       costs='pypsa-eur-sec/data/costs/costs_{planning_horizons}.csv',
                       cop_air_total="pypsa-eur-sec/resources/cop_air_total_{network}_s{simpl}_{clusters}.nc",
                       cop_soil_total="pypsa-eur-sec/resources/cop_soil_total_{network}_s{simpl}_{clusters}.nc"),                   
            output=['pypsa-eur-sec/results/test/prenetworks_bf/{network}_s{simpl}_{clusters}_lv{lv}__{sector_opts}_{planning_horizons}.nc']
        )
        import yaml
        with open('config.yaml') as f:
            snakemake.config = yaml.load(f)

    print(snakemake.input.network_p)
    logging.basicConfig(level=snakemake.config['logging_level'])

    options = snakemake.config["sector"]
    opts = snakemake.wildcards.sector_opts.split('-')
    
    year=snakemake.wildcards.planning_horizons

    n = pypsa.Network(snakemake.input.network,
                      override_component_attrs=override_component_attrs)
  
    n_p = pypsa.Network(snakemake.input.network_p,
                      override_component_attrs=override_component_attrs)

    add_brownfield(n,n_p)
    
    Nyears = n.snapshot_weightings.sum()/8760.
       
    costs = prepare_costs()   
    
    baseyear=snakemake.config['scenario']["planning_horizons"][0]
    
    add_power_capacities_installed_before_baseyear(n, year, baseyear, costs) # only the capacities with YearDecomissioning < year are added
   
    if "H" in opts:
        time_dep_hp_cop = options["time_dep_hp_cop"] 
        ashp_cop = xr.open_dataarray(snakemake.input.cop_air_total).T.to_pandas().reindex(index=n.snapshots)
        gshp_cop = xr.open_dataarray(snakemake.input.cop_soil_total).T.to_pandas().reindex(index=n.snapshots)
        add_heating_capacities_installed_before_baseyear(n, year, baseyear, ashp_cop, gshp_cop, time_dep_hp_cop, costs) # only the capacities with YearDecomissioning < year are added
    
    n.export_to_netcdf(snakemake.output[0])

   

