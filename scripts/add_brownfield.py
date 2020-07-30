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

from prepare_sector_network import prepare_costs

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

override_component_attrs["Link"].loc["build_year"] = ["integer","year",np.nan,"build year","Input (optional)"]
override_component_attrs["Link"].loc["lifetime"] = ["float","years",np.nan,"build year","Input (optional)"]
override_component_attrs["Generator"].loc["build_year"] = ["integer","year",np.nan,"build year","Input (optional)"]
override_component_attrs["Generator"].loc["lifetime"] = ["float","years",np.nan,"build year","Input (optional)"]
override_component_attrs["Store"].loc["build_year"] = ["integer","year",np.nan,"build year","Input (optional)"]
override_component_attrs["Store"].loc["lifetime"] = ["float","years",np.nan,"build year","Input (optional)"]

def add_brownfield(n, n_p, year):
    print("adding brownfield")

    #first, remove generators, links and stores that track CO2 or global EU values
    n_p.mremove("Generator", [index for index in n_p.generators.index.to_list() if 'ror' in index])
    n_p.mremove("Generator", ['EU fossil gas', 'fossil oil'] )
    n_p.mremove("Store", ['co2 atmosphere', 'co2 stored', 'EU gas Store'] )
    n_p.mremove("Link", ['co2 vent'] )

    if "H" in opts:
        n_p.mremove("Link", [index for index in n_p.links.index.to_list() if 'water tanks charger' in index])
        n_p.mremove("Link", [index for index in n_p.links.index.to_list() if 'water tanks discharger' in index])
    if "B" in opts:
         n_p.mremove("Store", ['EU biogas', 'EU solid biomass'])
         n_p.mremove("Link", ['biogas to gas'])
    if "I" in opts:
         n_p.mremove("Store", ['Fischer-Tropsch Store'])
         n_p.mremove("Link", ['process emissions' , 'gas for industry', 'solid biomass for industry'])
    if "T" in opts:
        n_p.mremove("Store", [index for index in n_p.stores.index.to_list() if 'battery storage' in index])
        n_p.mremove("Link", [index for index in n_p.links.index.to_list() if 'BEV charger' in index])
        n_p.mremove("Link", [index for index in n_p.links.index.to_list() if 'V2G' in index])

    previous_timestep=snakemake.config['scenario']['planning_horizons'][snakemake.config['scenario']['planning_horizons'].index(year)-1]
    previous_timesteps=snakemake.config['scenario']['planning_horizons'][0:snakemake.config['scenario']['planning_horizons'].index(year)]
    grouping_years=snakemake.config['existing_capacities']['grouping_years']


    ### GENERATORS ###
    # generators whose build_year + lifetime < year are removed
    n_p.mremove("Generator", [index for index in n_p.generators.index.to_list()
                              if  (n_p.generators.loc[index, 'build_year']+n_p.generators.loc[index, 'lifetime'] < int(year))])

    # remove generators if their optimized nominal capacity is lower than a threshold
    n_p.mremove("Generator", [index for index in n_p.generators.index.to_list()
                              if (n_p.generators.loc[index, 'p_nom_opt'] < snakemake.config['existing_capacities']['threshold_capacity'])])


    # generators whose capacity was optimized in the previous year are renamed and build year is added
    n_p.generators.index=np.where(n_p.generators.index.str[-4:].isin(previous_timesteps+grouping_years)==False,
                                  n_p.generators.index + '-' + previous_timestep,
                                  n_p.generators.index)
    n_p.generators.loc[[index for index in n_p.generators.index.to_list()
                        if previous_timestep in index], 'build_year']=int(previous_timestep)

    #add generators from previous step
    n.madd("Generator",
            n_p.generators.index,
            bus=n_p.generators.bus,
            carrier=n_p.generators.carrier,
            p_nom=n_p.generators.p_nom_opt,
            marginal_cost=n_p.generators.marginal_cost,
            capital_cost=n_p.generators.capital_cost,
            efficiency=n_p.generators.efficiency,
            p_max_pu=n_p.generators_t.p_max_pu,
            build_year=n_p.generators.build_year,
            lifetime=n_p.generators.lifetime)

    ### STORES ###
    # stores whose installationYear + lifetime < year are removed
    n_p.mremove("Store", [index for index in n_p.stores.index.to_list()
                          if (n_p.stores.loc[index, 'build_year']+n_p.stores.loc[index, 'lifetime'] < int(year))])

    # remove stores if their optimized nominal capacity is lower than a threshold
    n_p.mremove("Store", [index for index in n_p.stores.index.to_list()
                              if (n_p.stores.loc[index, 'e_nom_opt'] < snakemake.config['existing_capacities']['threshold_capacity'])])

    # stores whose capacity was optimized in the previous year are renamed and the build year is added
    n_p.stores.index=np.where(n_p.stores.index.str[-4:].isin(previous_timesteps+grouping_years)==False,
                                  n_p.stores.index + '-' + previous_timestep,
                                  n_p.stores.index)
    n_p.stores.loc[[index for index in n_p.stores.index.to_list()
                    if previous_timestep in index], 'build_year']=int(previous_timestep)
    #add stores from previous steps
    n.madd("Store",
            n_p.stores.index,
            bus=n_p.stores.bus,
            carrier=n_p.stores.carrier,
            e_nom=n_p.stores.e_nom_opt,
            e_cyclic=True,
            capital_cost=n_p.stores.capital_cost,
            build_year=n_p.stores.build_year,
            lifetime=n_p.stores.lifetime)

    ### LINKS ###
    # TODO: add_chp_constraint() in solve_network needs to be adjusted
    n_p.mremove("Link", [index for index in n_p.links.index.to_list() if 'CHP' in index])

    # links whose installationYear + lifetime < year are removed
    n_p.mremove("Link", [index for index in n_p.links.index.to_list()
                              if (n_p.links.loc[index, 'build_year']+n_p.links.loc[index, 'lifetime'] < int(year))])

    # delete links if their optimized nominal capacity is lower than a threshold
    n_p.mremove("Link", [index for index in n_p.links.index.to_list()
                              if (n_p.links.loc[index, 'p_nom_opt'] < snakemake.config['existing_capacities']['threshold_capacity'])])

    # links whose capacity was optimized in the previous year are renamed and the build year is added
    n_p.links.index=np.where(n_p.links.index.str[-4:].isin(previous_timesteps+grouping_years)==False,
                                  n_p.links.index + '-' + previous_timestep,
                                  n_p.links.index)
    n_p.links.loc[[index for index in n_p.links.index.to_list()
                  if previous_timestep in index], 'build_year']=int(previous_timestep)
    #add links from previous steps
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
            efficiency2=n_p.links.efficiency2,
            build_year=n_p.links.build_year,
            lifetime=n_p.links.lifetime)


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake
        snakemake = MockSnakemake(
            wildcards=dict(network='elec', simpl='', clusters='37', lv='1.0',
                           sector_opts='Co2L0-168H-T-H-B-I-solar3-dist1',
                           co2_budget_name='go',
                           planning_horizons='2030'),
            input=dict(network='pypsa-eur-sec/results/test/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}__{sector_opts}_{co2_budget_name}_{planning_horizons}.nc',
                       network_p='pypsa-eur-sec/results/test/postnetworks/{network}_s{simpl}_{clusters}_lv{lv}__{sector_opts}_{co2_budget_name}_2020.nc',
                       costs='pypsa-eur-sec/data/costs/costs_{planning_horizons}.csv',
                       cop_air_total="pypsa-eur-sec/resources/cop_air_total_{network}_s{simpl}_{clusters}.nc",
                       cop_soil_total="pypsa-eur-sec/resources/cop_soil_total_{network}_s{simpl}_{clusters}.nc"),
            output=['pypsa-eur-sec/results/test/prenetworks_brownfield/{network}_s{simpl}_{clusters}_lv{lv}__{sector_opts}_{planning_horizons}.nc']
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
#%%
    add_brownfield(n, n_p, year)

    Nyears = n.snapshot_weightings.sum()/8760.

    costs = prepare_costs(snakemake.input.costs,
                          snakemake.config['costs']['USD2013_to_EUR2013'],
                          snakemake.config['costs']['discountrate'],
                          Nyears)

    baseyear = snakemake.config['scenario']["planning_horizons"][0]


    n.export_to_netcdf(snakemake.output[0])
