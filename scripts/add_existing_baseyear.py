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


def add_power_capacities_installed_before_baseyear(n, year, baseyear, costs):
    """

    Parameters
    ----------
    n : network
    year : capacity that fulfills YearDecomissioning > year is added

    baseyear : capacity name will be e.g. "solar <baseyear"
               this allows adding the solar capacity that was installed before 
               2020 (baseyear) and is still alive in 2030 (year)
    costs : to read lifetime to estimate YearDecomissioning


    """
    print("adding power capacities installed before baseyear")
    
    
    ### add conventional capacities using 'powerplants.csv'
    df_agg = pd.read_csv('../pypsa-eur/resources/powerplants.csv', index_col=0)

    rename_fuel = {'Hard Coal':'coal',
                   'Lignite':'lignite',
                   'Nuclear':'nuclear',
                   'Oil':'oil',
                   'OCGT':'OCGT',
                   'CCGT':'CCGT',
                   'Natural Gas':'gas',}
    fueltype_to_drop = ['Hydro',
                        'Wind',
                        'Solar',
                        'Geothermal',
                        'Bioenergy',
                        'Waste',
                        'Other',
                        'CCGT, Thermal']
    technology_to_drop = ['Pv',
                          'Storage Technologies']

    df_agg.drop(df_agg.index[df_agg.Fueltype.isin(fueltype_to_drop)],inplace=True)
    df_agg.drop(df_agg.index[df_agg.Technology.isin(technology_to_drop)],inplace=True)
    df_agg.Fueltype = df_agg.Fueltype.map(rename_fuel)

    # add existing solar and wind capacities
    # source: https://www.irena.org/Statistics/Download-Data
    cc = pd.read_csv('data/Country_codes.csv', sep=',', index_col=-1)
    name_to_2code = dict(zip(cc['Country'].tolist(),
                             cc['2 letter code (ISO-3166-2)'].tolist()))
    for tech in ['solar', 'onwind', 'offwind']:
        df = pd.read_csv('data/existing_infrastructure/{}_capacity_IRENA.csv'.format(tech), 
                         sep=',', index_col=0)

        df.rename(index={'Czechia':'Czech Republic',
                         'UK':'United Kingdom',
                         'Bosnia Herzg':'Bosnia Herzegovina',
                         'North Macedonia': 'Macedonia'}, inplace=True)

        df.rename(index=lambda country : name_to_2code[country], inplace=True)

        # calculate yearly differences
        df.insert(loc=0, value=.0, column='1999')
        df = df.diff(axis=1).drop('1999', axis=1)
        df = df.clip(lower=0)
        df.replace(to_replace=0.0, value=np.nan, inplace=True)

        for year in df.columns:
            for country in df.index:
                if df.notnull().loc[country,year]:
                    df_agg = df_agg.append({'Fueltype':tech,
                                            'Country':country,
                                            'Capacity':df.loc[country,year],
                                            'YearCommissioned':int(year),
                                            'YearDecommissioning':int(float(year)+costs.at[tech, 'lifetime'])},
                                             ignore_index=True)
           
        # powerplants already include an estimated decommissining year wich has been
        # calculated based on pm config, check if we want to make this calculation more
        # transparent, e.g., by loading lifetime assumptions from costs.csv into pm config
        # check how pm is dealing with retrofitting
        # pm.get_config()['fuel_to_lifetime']
        # index = df_agg.Fueltype.map(lf)+df_agg.YearCommissioned > 2020

        index=df_agg.YearDecommissioning > int(year)
        df = df_agg[index].pivot_table(index='Country', columns='Fueltype',
                                       values='Capacity', aggfunc='sum')

        nodes=set([node[0:2] for node in n.buses.index[n.buses.carrier == "AC"]])

    for carrier in ['coal', 'oil', 'uranium']:
        n.add("Bus",
              "EU " + carrier,
              carrier=carrier)
    
    for node in nodes:
        #if a country has more than one node, selects the first one
        bus_selected=[bus for bus in n.buses.index[n.buses.carrier == "AC"] if bus[0:2]==node][0]  
        for generator,carrier in [("OCGT","gas"), 
                                  ("CCGT", "gas"), 
                                  ("coal", "coal"), 
                                  ("oil","oil"), 
                                  ("nuclear","uranium")]: 
            if node in df.index and not np.isnan(df.loc[node, generator]):           
                n.add("Link",
                      bus_selected + " " + generator +" <" + baseyear,
                      bus0="EU " + carrier,
                      bus1=bus_selected,
                      bus2="co2 atmosphere",                        
                      marginal_cost=costs.at[generator,'efficiency']*costs.at[generator,'VOM'], #NB: VOM is per MWel
                      capital_cost=costs.at[generator,'efficiency']*costs.at[generator,'fixed'], #NB: fixed cost is per MWel
                      p_nom=df.loc[node, generator],
                      efficiency=costs.at[generator,'efficiency'],
                      efficiency2=costs.at[carrier,'CO2 intensity'])        
        
        for generator in ['solar', 'onwind', 'offwind']:
            try: 
                if not np.isnan(df.loc[node, generator]):
                    if generator =='offwind':
                        p_max_pu=n.generators_t.p_max_pu[bus_selected + ' offwind-ac']
                    else:
                        p_max_pu=n.generators_t.p_max_pu[bus_selected + ' ' + generator]
            
                    n.add("Generator", 
                          bus_selected + ' ' + generator +" <"+ baseyear,
                          bus=bus_selected,
                          carrier=generator,                  
                          p_nom=df.loc[node, generator],
                          marginal_cost=costs.at[generator,'VOM'],
                          capital_cost=costs.at[generator,'fixed'],
                          efficiency=costs.at[generator, 'efficiency'],
                          p_max_pu=p_max_pu)
            except:
                    print("No capacity installed before base year of " + generator + " is alive in node " + node)

    # delete generators if their lifetime is over and p_nom=0
    n.mremove("Generator", [index for index in n.generators.index.to_list() if '<'+baseyear in index and n.generators.p_nom[index]==0])
    n.mremove("Link", [index for index in n.links.index.to_list() if '<'+baseyear in index and n.links.p_nom[index]==0])
    
def add_heating_capacities_installed_before_baseyear(n, year, baseyear, ashp_cop, gshp_cop, time_dep_hp_cop, costs, default_lifetime):

    """

    Parameters
    ----------
    n : network
    year : capacity that fulfills YearDecomissioning > year is added 
    baseyear : capacity name will be e.g. "gas boiler <baseyear"
                this allows adding the gas boiler capacity that was installed 
                before 2020 (baseyear) and is still alive in 2030 (year)
    
    linear decomissioning of heating capacities from 2020 to 2045 is 
    currently assumed
    
    heating capacities split between residential and services proportional
    to heating load in both 
    50% capacities in rural busess 50% in urban buses
    """    
    print("adding heating capacities installed before baseyear")        
                
    # Add existing heating capacities, data comes from the study
    # "Mapping and analyses of the current and future (2020 - 2030) 
    # heating/cooling fuel deployment (fossil/renewables) "
    # https://ec.europa.eu/energy/studies/mapping-and-analyses-current-and-future-2020-2030-heatingcooling-fuel-deployment_en?redir=1
    # file: "WP2_DataAnnex_1_BuildingTechs_ForPublication_201603.xls" -> "existing_heating_raw.csv".
        
    # retrieve existing heating capacities
    techs = ['gas boiler',
             'oil boiler',
             'resistive heater',
             'air heat pump',
             'ground heat pump']
    df = pd.read_csv('data/existing_infrastructure/existing_heating_raw.csv',
                     index_col=0,
                     header=0)
    # data for Albania, Montenegro and Macedonia not included in database 
    df.loc['Albania']=np.nan
    df.loc['Montenegro']=np.nan
    df.loc['Macedonia']=np.nan
    df.fillna(0, inplace=True)
    df *= 1e3 # GW to MW
        
    cc = pd.read_csv('data/Country_codes.csv', sep=',', index_col=-1)
    name_to_2code = dict(zip(cc['Country'].tolist(),
                         cc['2 letter code (ISO-3166-2)'].tolist()))
    df.rename(index=lambda country : name_to_2code[country], inplace=True)

        
    # coal and oil boilers are assimilated to oil boilers        
    df['oil boiler'] =df['oil boiler'] + df['coal boiler'] 
    df.drop(['coal boiler'], axis=1, inplace=True)
        
    # rename countries with network buses names
    nodes_elec=[node for node in n.buses.index[n.buses.carrier == "AC"]]
    name_to_busname={ index : [node for node in nodes_elec if index in node][0] for index in df.index}
    df.rename(index=lambda country : name_to_busname[country], inplace=True)

    # split existing capacities between residential and services 
    # proportional to energy demand           
    ratio_residential=pd.Series([(n.loads_t.p_set.sum()['{} residential rural heat'.format(node)] /
               (n.loads_t.p_set.sum()['{} residential rural heat'.format(node)] +
                n.loads_t.p_set.sum()['{} services rural heat'.format(node)] ))
                   for node in df.index], index=df.index)
       
    for tech in techs: 
        df['residential ' + tech] = df[tech]*ratio_residential
        df['services ' + tech] = df[tech]*(1-ratio_residential)
               
    nodes={}
    p_nom={}        
    for name in ["residential rural", 
                 "services rural",
                 "residential urban decentral",
                 "services urban decentral", 
                 "urban central"]:
            
        name_type = "central" if name == "urban central" else "decentral"            
        nodes[name] = pd.Index([index[0:5] for index in n.buses.index[n.buses.index.str.contains(name) & n.buses.index.str.contains('heat')]])
        heat_pump_type = "air" if "urban" in name else "ground"
        heat_type= "residential" if "residential" in name else "services"
            
        if name == "urban central":
            p_nom[name]=df['air heat pump'][nodes[name]]
        else:
            p_nom[name] =  df['{} {} heat pump'.format(heat_type, heat_pump_type)][nodes[name]]
            
        # Add heat pumps        
        costs_name = "{} {}-sourced heat pump".format("decentral", heat_pump_type)


        cop = {"air" : ashp_cop, "ground" : gshp_cop}
        efficiency = cop[heat_pump_type][nodes[name]] if time_dep_hp_cop else costs.at[costs_name,'efficiency']

        n.madd("Link",
               nodes[name],
               suffix=" {} {} heat pump <{}".format(name,heat_pump_type, baseyear),
               bus0=nodes[name],
               bus1=nodes[name] + " " + name + " heat",
               carrier="{} {} heat pump".format(name,heat_pump_type),
               efficiency=efficiency,
               capital_cost=costs.at[costs_name,'efficiency']*costs.at[costs_name,'fixed'],
               p_nom=p_nom[name]*max(0,(int(baseyear)+default_lifetime-int(year))/default_lifetime)) #decomissioning happens lineary from now to 25 years lter  
           
        # add resistive heater, gas boilers and oil boilers 
        # (50% capacities to rural buses, 50% to urban buses)
        n.madd("Link",
               nodes[name],
               suffix= " " + name + " resistive heater <{}".format(baseyear),
               bus0=nodes[name],
               bus1=nodes[name] + " " + name + " heat",
               carrier=name + " resistive heater",
               efficiency=costs.at[name_type + ' resistive heater','efficiency'],
               capital_cost=costs.at[name_type + ' resistive heater','efficiency']*costs.at[name_type + ' resistive heater','fixed'],
               p_nom=0.5*df['{} resistive heater'.format(heat_type)][nodes[name]]*max(0,(int(baseyear)+default_lifetime-int(year))/default_lifetime))

        n.madd("Link",
               nodes[name],
               suffix= " " + name + " gas boiler <{}".format(baseyear),
               bus0=["EU gas"]*len(nodes[name]),
               bus1=nodes[name] + " " + name + " heat",
               bus2="co2 atmosphere",
               carrier=name + " gas boiler",
               efficiency=costs.at[name_type + ' gas boiler','efficiency'],
               efficiency2=costs.at['gas','CO2 intensity'],
               capital_cost=costs.at[name_type + ' gas boiler','efficiency']*costs.at[name_type + ' gas boiler','fixed'],
               p_nom=0.5*df['{} gas boiler'.format(heat_type)][nodes[name]]*max(0,(int(baseyear)+default_lifetime-int(year))/default_lifetime))

        n.madd("Link",
               nodes[name],
               suffix=" " + name + " oil boiler <{}".format(baseyear),
               bus0=["EU oil"]*len(nodes[name]),
               bus1=nodes[name] + " " + name + " heat",
               bus2="co2 atmosphere",
               carrier=name + " oil boiler",
               efficiency=costs.at['decentral oil boiler','efficiency'],
               efficiency2=costs.at['oil','CO2 intensity'],
               capital_cost=costs.at['decentral oil boiler','efficiency']*costs.at['decentral oil boiler','fixed'],
               p_nom=0.5*df['{} oil boiler'.format(heat_type)][nodes[name]]*max(0,(int(baseyear)+default_lifetime-int(year))/default_lifetime))
            
        # delete links with p_nom=nan corresponding to extra nodes in country
        n.mremove("Link", [index for index in n.links.index.to_list() if baseyear in index and np.isnan(n.links.p_nom[index])])
        
        # delete links if their lifetime is over and p_nom=0
        n.mremove("Link", [index for index in n.links.index.to_list() if '<'+baseyear in index and n.links.p_nom[index]==0])

if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake
        snakemake = MockSnakemake(
            wildcards=dict(network='elec', simpl='', clusters='37', lv='1.0', 
                           sector_opts='Co2L0-168H-T-H-B-I-solar3-dist1',
                           co2_budget_name='go',
                           planning_horizons='2020'),
            input=dict(network='pypsa-eur-sec/results/test/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}__{sector_opts}_{co2_budget_name}_{planning_horizons}.nc', 
                       costs='pypsa-eur-sec/data/costs/costs_{planning_horizons}.csv',
                       cop_air_total="pypsa-eur-sec/resources/cop_air_total_{network}_s{simpl}_{clusters}.nc",
                       cop_soil_total="pypsa-eur-sec/resources/cop_soil_total_{network}_s{simpl}_{clusters}.nc"),                      
            output=['pypsa-eur-sec/results/test/prenetworks_bf/{network}_s{simpl}_{clusters}_lv{lv}__{sector_opts}_{planning_horizons}.nc'],
        )
        import yaml
        with open('config.yaml') as f:
            snakemake.config = yaml.load(f)


    logging.basicConfig(level=snakemake.config['logging_level'])

    options = snakemake.config["sector"]
    opts = snakemake.wildcards.sector_opts.split('-')

    baseyear= snakemake.config['scenario']["planning_horizons"][0]
    
    n = pypsa.Network(snakemake.input.network,
                      override_component_attrs=override_component_attrs)
   
    Nyears = n.snapshot_weightings.sum()/8760.
    costs = prepare_costs(snakemake.input.costs,
                          snakemake.config['costs']['USD2013_to_EUR2013'],
                          snakemake.config['costs']['discountrate'],
                          Nyears)

    
    add_power_capacities_installed_before_baseyear(n, baseyear, baseyear, costs)
        
    if "H" in opts:
        time_dep_hp_cop = options["time_dep_hp_cop"] 
        ashp_cop = xr.open_dataarray(snakemake.input.cop_air_total).T.to_pandas().reindex(index=n.snapshots)
        gshp_cop = xr.open_dataarray(snakemake.input.cop_soil_total).T.to_pandas().reindex(index=n.snapshots)
        default_lifetime = snakemake.config['costs']['lifetime']
        add_heating_capacities_installed_before_baseyear(n, baseyear, baseyear, ashp_cop, gshp_cop, time_dep_hp_cop, costs, default_lifetime)
   
    n.export_to_netcdf(snakemake.output[0])


    


   

