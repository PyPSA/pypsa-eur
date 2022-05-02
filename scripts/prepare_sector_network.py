# coding: utf-8

import pypsa
import re
import os

import pandas as pd
import numpy as np
import xarray as xr
import networkx as nx

from itertools import product
from scipy.stats import beta
from vresutils.costdata import annuity

from build_energy_totals import build_eea_co2, build_eurostat_co2, build_co2_totals
from helper import override_component_attrs, generate_periodic_profiles

from networkx.algorithms.connectivity.edge_augmentation import k_edge_augmentation
from networkx.algorithms import complement
from pypsa.geo import haversine_pts

import logging
logger = logging.getLogger(__name__)

from types import SimpleNamespace
spatial = SimpleNamespace()


def define_spatial(nodes, options):
    """
    Namespace for spatial

    Parameters
    ----------
    nodes : list-like
    """

    global spatial

    spatial.nodes = nodes

    # biomass

    spatial.biomass = SimpleNamespace()

    if options["biomass_transport"]:
        spatial.biomass.nodes = nodes + " solid biomass"
        spatial.biomass.locations = nodes
        spatial.biomass.industry = nodes + " solid biomass for industry"
        spatial.biomass.industry_cc = nodes + " solid biomass for industry CC"
    else:
        spatial.biomass.nodes = ["EU solid biomass"]
        spatial.biomass.locations = ["EU"]
        spatial.biomass.industry = ["solid biomass for industry"]
        spatial.biomass.industry_cc = ["solid biomass for industry CC"]

    spatial.biomass.df = pd.DataFrame(vars(spatial.biomass), index=nodes)

    # co2

    spatial.co2 = SimpleNamespace()

    if options["co2_network"]:
        spatial.co2.nodes = nodes + " co2 stored"
        spatial.co2.locations = nodes
        spatial.co2.vents = nodes + " co2 vent"
    else:
        spatial.co2.nodes = ["co2 stored"]
        spatial.co2.locations = ["EU"]
        spatial.co2.vents = ["co2 vent"]

    spatial.co2.df = pd.DataFrame(vars(spatial.co2), index=nodes)

    # gas

    spatial.gas = SimpleNamespace()

    if options["gas_network"]:
        spatial.gas.nodes = nodes + " gas"
        spatial.gas.locations = nodes
        spatial.gas.biogas = nodes + " biogas"
        spatial.gas.industry = nodes + " gas for industry"
        spatial.gas.industry_cc = nodes + " gas for industry CC"
        spatial.gas.biogas_to_gas = nodes + " biogas to gas"
    else:
        spatial.gas.nodes = ["EU gas"]
        spatial.gas.locations = ["EU"]
        spatial.gas.biogas = ["EU biogas"]
        spatial.gas.industry = ["gas for industry"]
        spatial.gas.industry_cc = ["gas for industry CC"]
        spatial.gas.biogas_to_gas = ["EU biogas to gas"]

    spatial.gas.df = pd.DataFrame(vars(spatial.gas), index=nodes)

    # oil
    spatial.oil = SimpleNamespace()
    spatial.oil.nodes = ["EU oil"]
    spatial.oil.locations = ["EU"]

    # uranium
    spatial.uranium = SimpleNamespace()
    spatial.uranium.nodes = ["EU uranium"]
    spatial.uranium.locations = ["EU"]

    # coal
    spatial.coal = SimpleNamespace()
    spatial.coal.nodes = ["EU coal"]
    spatial.coal.locations = ["EU"]

    # lignite
    spatial.lignite = SimpleNamespace()
    spatial.lignite.nodes = ["EU lignite"]
    spatial.lignite.locations = ["EU"]

    return spatial


from types import SimpleNamespace
spatial = SimpleNamespace()


def emission_sectors_from_opts(opts):

    sectors = ["electricity"]
    if "T" in opts:
        sectors += [
            "rail non-elec",
            "road non-elec"
        ]
    if "H" in opts:
        sectors += [
            "residential non-elec",
            "services non-elec"
        ]
    if "I" in opts:
        sectors += [
            "industrial non-elec",
            "industrial processes",
            "domestic aviation",
            "international aviation",
            "domestic navigation",
            "international navigation"
        ]
    if "A" in opts:
        sectors += [
            "agriculture"
        ]

    return sectors


def get(item, investment_year=None):
    """Check whether item depends on investment year"""
    if isinstance(item, dict):
        return item[investment_year]
    else:
        return item


def co2_emissions_year(countries, opts, year):
    """
    Calculate CO2 emissions in one specific year (e.g. 1990 or 2018).
    """

    eea_co2 = build_eea_co2(year)

    # TODO: read Eurostat data from year > 2014
    # this only affects the estimation of CO2 emissions for BA, RS, AL, ME, MK
    if year > 2014:
        eurostat_co2 = build_eurostat_co2(year=2014)
    else:
        eurostat_co2 = build_eurostat_co2(year)

    co2_totals = build_co2_totals(eea_co2, eurostat_co2)

    sectors = emission_sectors_from_opts(opts)

    co2_emissions = co2_totals.loc[countries, sectors].sum().sum()

    # convert MtCO2 to GtCO2
    co2_emissions *= 0.001

    return co2_emissions


# TODO: move to own rule with sector-opts wildcard?
def build_carbon_budget(o, fn):
    """
    Distribute carbon budget following beta or exponential transition path.
    """
    # opts?

    if "be" in o:
        #beta decay
        carbon_budget = float(o[o.find("cb")+2:o.find("be")])
        be = float(o[o.find("be")+2:])
    if "ex" in o:
        #exponential decay
        carbon_budget = float(o[o.find("cb")+2:o.find("ex")])
        r = float(o[o.find("ex")+2:])

    countries = n.buses.country.dropna().unique()

    e_1990 = co2_emissions_year(countries, opts, year=1990)

    #emissions at the beginning of the path (last year available 2018)
    e_0 = co2_emissions_year(countries, opts, year=2018)

    planning_horizons = snakemake.config['scenario']['planning_horizons']
    t_0 = planning_horizons[0]

    if "be" in o:

        # final year in the path
        t_f = t_0 + (2 * carbon_budget / e_0).round(0)

        def beta_decay(t):
            cdf_term = (t - t_0) / (t_f - t_0)
            return (e_0 / e_1990) * (1 - beta.cdf(cdf_term, be, be))

        #emissions (relative to 1990)
        co2_cap = pd.Series({t: beta_decay(t) for t in planning_horizons}, name=o)

    if "ex" in o:

        T = carbon_budget / e_0
        m = (1 + np.sqrt(1 + r * T)) / T

        def exponential_decay(t):
            return (e_0 / e_1990) * (1 + (m + r) * (t - t_0)) * np.exp(-m * (t - t_0))

        co2_cap = pd.Series({t: exponential_decay(t) for t in planning_horizons}, name=o)

    # TODO log in Snakefile
    if not os.path.exists(fn):
        os.makedirs(fn)
    co2_cap.to_csv(fn, float_format='%.3f')


def add_lifetime_wind_solar(n, costs):
    """Add lifetime for solar and wind generators."""
    for carrier in ['solar', 'onwind', 'offwind']:
        gen_i = n.generators.index.str.contains(carrier)
        n.generators.loc[gen_i, "lifetime"] = costs.at[carrier, 'lifetime']


def haversine(p):
    coord0 = n.buses.loc[p.bus0, ['x', 'y']].values
    coord1 = n.buses.loc[p.bus1, ['x', 'y']].values
    return 1.5 * haversine_pts(coord0, coord1)


def create_network_topology(n, prefix, carriers=["DC"], connector=" -> ", bidirectional=True):
    """
    Create a network topology from transmission lines and link carrier selection.

    Parameters
    ----------
    n : pypsa.Network
    prefix : str
    carriers : list-like
    connector : str
    bidirectional : bool, default True
        True: one link for each connection
        False: one link for each connection and direction (back and forth)

    Returns
    -------
    pd.DataFrame with columns bus0, bus1, length, underwater_fraction
    """

    ln_attrs = ["bus0", "bus1", "length"]
    lk_attrs = ["bus0", "bus1", "length", "underwater_fraction"]
    lk_attrs = n.links.columns.intersection(lk_attrs)

    candidates = pd.concat([
        n.lines[ln_attrs],
        n.links.loc[n.links.carrier.isin(carriers), lk_attrs]
    ]).fillna(0)

    # base network topology purely on location not carrier
    candidates["bus0"] = candidates.bus0.map(n.buses.location)
    candidates["bus1"] = candidates.bus1.map(n.buses.location)

    positive_order = candidates.bus0 < candidates.bus1
    candidates_p = candidates[positive_order]
    swap_buses = {"bus0": "bus1", "bus1": "bus0"}
    candidates_n = candidates[~positive_order].rename(columns=swap_buses)
    candidates = pd.concat([candidates_p, candidates_n])

    def make_index(c):
        return prefix + c.bus0 + connector + c.bus1

    topo = candidates.groupby(["bus0", "bus1"], as_index=False).mean()
    topo.index = topo.apply(make_index, axis=1)

    if not bidirectional:
        topo_reverse = topo.copy()
        topo_reverse.rename(columns=swap_buses, inplace=True)
        topo_reverse.index = topo_reverse.apply(make_index, axis=1)
        topo = pd.concat([topo, topo_reverse])

    return topo


# TODO merge issue with PyPSA-Eur
def update_wind_solar_costs(n, costs):
    """
    Update costs for wind and solar generators added with pypsa-eur to those
    cost in the planning year
    """

    #NB: solar costs are also manipulated for rooftop
    #when distribution grid is inserted
    n.generators.loc[n.generators.carrier=='solar', 'capital_cost'] = costs.at['solar-utility', 'fixed']

    n.generators.loc[n.generators.carrier=='onwind', 'capital_cost'] = costs.at['onwind', 'fixed']

    #for offshore wind, need to calculated connection costs

    #assign clustered bus
    #map initial network -> simplified network
    busmap_s = pd.read_csv(snakemake.input.busmap_s, index_col=0).squeeze()
    busmap_s.index = busmap_s.index.astype(str)
    busmap_s = busmap_s.astype(str)
    #map simplified network -> clustered network
    busmap = pd.read_csv(snakemake.input.busmap, index_col=0).squeeze()
    busmap.index = busmap.index.astype(str)
    busmap = busmap.astype(str)
    #map initial network -> clustered network
    clustermaps = busmap_s.map(busmap)

    #code adapted from pypsa-eur/scripts/add_electricity.py
    for connection in ['dc', 'ac']:
        tech = "offwind-" + connection
        profile = snakemake.input['profile_offwind_' + connection]
        with xr.open_dataset(profile) as ds:
            underwater_fraction = ds['underwater_fraction'].to_pandas()
            connection_cost = (snakemake.config['costs']['lines']['length_factor'] *
                               ds['average_distance'].to_pandas() *
                               (underwater_fraction *
                                costs.at[tech + '-connection-submarine', 'fixed'] +
                                (1. - underwater_fraction) *
                                costs.at[tech + '-connection-underground', 'fixed']))

            #convert to aggregated clusters with weighting
            weight = ds['weight'].to_pandas()

            #e.g. clusters == 37m means that VRE generators are left
            #at clustering of simplified network, but that they are
            #connected to 37-node network
            if snakemake.wildcards.clusters[-1:] == "m":
                genmap = busmap_s
            else:
                genmap = clustermaps

            connection_cost = (connection_cost*weight).groupby(genmap).sum()/weight.groupby(genmap).sum()

            capital_cost = (costs.at['offwind', 'fixed'] +
                            costs.at[tech + '-station', 'fixed'] +
                            connection_cost)

            logger.info("Added connection cost of {:0.0f}-{:0.0f} Eur/MW/a to {}"
                        .format(connection_cost[0].min(), connection_cost[0].max(), tech))

            n.generators.loc[n.generators.carrier==tech, 'capital_cost'] = capital_cost.rename(index=lambda node: node + ' ' + tech)


def add_carrier_buses(n, carrier, nodes=None):
    """
    Add buses to connect e.g. coal, nuclear and oil plants
    """

    if nodes is None:
        nodes = vars(spatial)[carrier].nodes
    location = vars(spatial)[carrier].locations

    # skip if carrier already exists
    if carrier in n.carriers.index:
        return

    if not isinstance(nodes, pd.Index):
        nodes = pd.Index(nodes)

    n.add("Carrier", carrier)

    n.madd("Bus",
        nodes,
        location=location,
        carrier=carrier
    )

    #capital cost could be corrected to e.g. 0.2 EUR/kWh * annuity and O&M
    n.madd("Store",
        nodes + " Store",
        bus=nodes,
        e_nom_extendable=True,
        e_cyclic=True,
        carrier=carrier,
    )

    n.madd("Generator",
        nodes,
        bus=nodes,
        p_nom_extendable=True,
        carrier=carrier,
        marginal_cost=costs.at[carrier, 'fuel']
    )


# TODO: PyPSA-Eur merge issue
def remove_elec_base_techs(n):
    """remove conventional generators (e.g. OCGT) and storage units (e.g. batteries and H2)
    from base electricity-only network, since they're added here differently using links
    """

    for c in n.iterate_components(snakemake.config["pypsa_eur"]):
        to_keep = snakemake.config["pypsa_eur"][c.name]
        to_remove = pd.Index(c.df.carrier.unique()).symmetric_difference(to_keep)
        print("Removing", c.list_name, "with carrier", to_remove)
        names = c.df.index[c.df.carrier.isin(to_remove)]
        n.mremove(c.name, names)
        n.carriers.drop(to_remove, inplace=True, errors="ignore")


# TODO: PyPSA-Eur merge issue
def remove_non_electric_buses(n):
    """
    remove buses from pypsa-eur with carriers which are not AC buses
    """
    print("drop buses from PyPSA-Eur with carrier: ", n.buses[~n.buses.carrier.isin(["AC", "DC"])].carrier.unique())
    n.buses = n.buses[n.buses.carrier.isin(["AC", "DC"])]


def patch_electricity_network(n):
    remove_elec_base_techs(n)
    remove_non_electric_buses(n)
    update_wind_solar_costs(n, costs)
    n.loads["carrier"] = "electricity"
    n.buses["location"] = n.buses.index
    # remove trailing white space of load index until new PyPSA version after v0.18.
    n.loads.rename(lambda x: x.strip(), inplace=True)
    n.loads_t.p_set.rename(lambda x: x.strip(), axis=1, inplace=True)


def add_co2_tracking(n, options):

    # minus sign because opposite to how fossil fuels used:
    # CH4 burning puts CH4 down, atmosphere up
    n.add("Carrier", "co2",
          co2_emissions=-1.)

    # this tracks CO2 in the atmosphere
    n.add("Bus",
        "co2 atmosphere",
        location="EU",
        carrier="co2"
    )

    # can also be negative
    n.add("Store",
        "co2 atmosphere",
        e_nom_extendable=True,
        e_min_pu=-1,
        carrier="co2",
        bus="co2 atmosphere"
    )

    # this tracks CO2 stored, e.g. underground
    n.madd("Bus",
        spatial.co2.nodes,
        location=spatial.co2.locations,
        carrier="co2 stored"
    )

    n.madd("Store",
        spatial.co2.nodes,
        e_nom_extendable=True,
        e_nom_max=np.inf,
        capital_cost=options['co2_sequestration_cost'],
        carrier="co2 stored",
        bus=spatial.co2.nodes
    )

    if options['co2_vent']:

        n.madd("Link",
            spatial.co2.vents,
            bus0=spatial.co2.nodes,
            bus1="co2 atmosphere",
            carrier="co2 vent",
            efficiency=1.,
            p_nom_extendable=True
        )


def add_co2_network(n, costs):

    logger.info("Adding CO2 network.")
    co2_links = create_network_topology(n, "CO2 pipeline ")

    cost_onshore = (1 - co2_links.underwater_fraction) * costs.at['CO2 pipeline', 'fixed'] * co2_links.length
    cost_submarine = co2_links.underwater_fraction * costs.at['CO2 submarine pipeline', 'fixed'] * co2_links.length
    capital_cost = cost_onshore + cost_submarine

    n.madd("Link",
        co2_links.index,
        bus0=co2_links.bus0.values + " co2 stored",
        bus1=co2_links.bus1.values + " co2 stored",
        p_min_pu=-1,
        p_nom_extendable=True,
        length=co2_links.length.values,
        capital_cost=capital_cost.values,
        carrier="CO2 pipeline",
        lifetime=costs.at['CO2 pipeline', 'lifetime']
    )


def add_dac(n, costs):

    heat_carriers = ["urban central heat", "services urban decentral heat"]
    heat_buses = n.buses.index[n.buses.carrier.isin(heat_carriers)]
    locations = n.buses.location[heat_buses]

    efficiency2 = -(costs.at['direct air capture', 'electricity-input'] + costs.at['direct air capture', 'compression-electricity-input'])
    efficiency3 = -(costs.at['direct air capture', 'heat-input'] - costs.at['direct air capture', 'compression-heat-output'])

    n.madd("Link",
        heat_buses.str.replace(" heat", " DAC"),
        bus0="co2 atmosphere",
        bus1=spatial.co2.df.loc[locations, "nodes"].values,
        bus2=locations.values,
        bus3=heat_buses,
        carrier="DAC",
        capital_cost=costs.at['direct air capture', 'fixed'],
        efficiency=1.,
        efficiency2=efficiency2,
        efficiency3=efficiency3,
        p_nom_extendable=True,
        lifetime=costs.at['direct air capture', 'lifetime']
    )


def add_co2limit(n, Nyears=1., limit=0.):

    logger.info(f"Adding CO2 budget limit as per unit of 1990 levels of {limit}")

    countries = n.buses.country.dropna().unique()

    sectors = emission_sectors_from_opts(opts)

    # convert Mt to tCO2
    co2_totals = 1e6 * pd.read_csv(snakemake.input.co2_totals_name, index_col=0)

    co2_limit = co2_totals.loc[countries, sectors].sum().sum()

    co2_limit *= limit * Nyears

    n.add("GlobalConstraint",
        "CO2Limit",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=co2_limit
    )

# TODO PyPSA-Eur merge issue
def average_every_nhours(n, offset):
    logger.info(f'Resampling the network to {offset}')
    m = n.copy(with_time=False)

    snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name+"_t")
        for k, df in c.pnl.items():
            if not df.empty:
                if c.list_name == "stores" and k == "e_max_pu":
                    pnl[k] = df.resample(offset).min()
                elif c.list_name == "stores" and k == "e_min_pu":
                    pnl[k] = df.resample(offset).max()
                else:
                    pnl[k] = df.resample(offset).mean()

    return m


def cycling_shift(df, steps=1):
    """Cyclic shift on index of pd.Series|pd.DataFrame by number of steps"""
    df = df.copy()
    new_index = np.roll(df.index, steps)
    df.values[:] = df.reindex(index=new_index).values
    return df


# TODO checkout PyPSA-Eur script
def prepare_costs(cost_file, USD_to_EUR, discount_rate, Nyears, lifetime):

    #set all asset costs and other parameters
    costs = pd.read_csv(cost_file, index_col=[0,1]).sort_index()

    #correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD"), "value"] *= USD_to_EUR

    #min_count=1 is important to generate NaNs which are then filled by fillna
    costs = costs.loc[:, "value"].unstack(level=1).groupby("technology").sum(min_count=1)
    costs = costs.fillna({"CO2 intensity" : 0,
                          "FOM" : 0,
                          "VOM" : 0,
                          "discount rate" : discount_rate,
                          "efficiency" : 1,
                          "fuel" : 0,
                          "investment" : 0,
                          "lifetime" : lifetime
    })

    annuity_factor = lambda v: annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100
    costs["fixed"] = [annuity_factor(v) * v["investment"] * Nyears for i, v in costs.iterrows()]

    return costs


def add_generation(n, costs):

    logger.info("adding electricity generation")

    nodes = pop_layout.index

    fallback = {"OCGT": "gas"}
    conventionals = options.get("conventional_generation", fallback)

    for generator, carrier in conventionals.items():


        carrier_nodes = vars(spatial)[carrier].nodes

        add_carrier_buses(n, carrier, carrier_nodes)

        n.madd("Link",
            nodes + " " + generator,
            bus0=carrier_nodes,
            bus1=nodes,
            bus2="co2 atmosphere",
            marginal_cost=costs.at[generator, 'efficiency'] * costs.at[generator, 'VOM'], #NB: VOM is per MWel
            capital_cost=costs.at[generator, 'efficiency'] * costs.at[generator, 'fixed'], #NB: fixed cost is per MWel
            p_nom_extendable=True,
            carrier=generator,
            efficiency=costs.at[generator, 'efficiency'],
            efficiency2=costs.at[carrier, 'CO2 intensity'],
            lifetime=costs.at[generator, 'lifetime']
        )


def add_wave(n, wave_cost_factor):

    # TODO: handle in Snakefile
    wave_fn = "data/WindWaveWEC_GLTB.xlsx"

    #in kW
    capacity = pd.Series({"Attenuator": 750,
                          "F2HB": 1000,
                          "MultiPA": 600})

    #in EUR/MW
    annuity_factor = annuity(25,0.07) + 0.03
    costs = 1e6 * wave_cost_factor * annuity_factor * pd.Series({"Attenuator": 2.5,
                                                                 "F2HB": 2,
                                                                 "MultiPA": 1.5})

    sheets = pd.read_excel(wave_fn, sheet_name=["FirthForth", "Hebrides"],
                           usecols=["Attenuator", "F2HB", "MultiPA"],
                           index_col=0, skiprows=[0], parse_dates=True)

    wave = pd.concat([sheets[l].divide(capacity, axis=1) for l in locations],
                     keys=locations,
                     axis=1)

    for wave_type in costs.index:
        n.add("Generator",
            "Hebrides " + wave_type,
            bus="GB4 0", # TODO this location is hardcoded
            p_nom_extendable=True,
            carrier="wave",
            capital_cost=costs[wave_type],
            p_max_pu=wave["Hebrides", wave_type]
        )


def insert_electricity_distribution_grid(n, costs):
    # TODO pop_layout?
    # TODO options?

    print("Inserting electricity distribution grid with investment cost factor of",
          options['electricity_distribution_grid_cost_factor'])

    nodes = pop_layout.index

    cost_factor = options['electricity_distribution_grid_cost_factor']

    n.madd("Bus",
        nodes + " low voltage",
        location=nodes,
        carrier="low voltage"
    )

    n.madd("Link",
        nodes + " electricity distribution grid",
        bus0=nodes,
        bus1=nodes + " low voltage",
        p_nom_extendable=True,
        p_min_pu=-1,
        carrier="electricity distribution grid",
        efficiency=1,
        lifetime=costs.at['electricity distribution grid', 'lifetime'],
        capital_cost=costs.at['electricity distribution grid', 'fixed'] * cost_factor
    )

    # this catches regular electricity load and "industry electricity" and
    # "agriculture machinery electric" and "agriculture electricity"
    loads = n.loads.index[n.loads.carrier.str.contains("electric")]
    n.loads.loc[loads, "bus"] += " low voltage"

    bevs = n.links.index[n.links.carrier == "BEV charger"]
    n.links.loc[bevs, "bus0"] += " low voltage"

    v2gs = n.links.index[n.links.carrier == "V2G"]
    n.links.loc[v2gs, "bus1"] += " low voltage"

    hps = n.links.index[n.links.carrier.str.contains("heat pump")]
    n.links.loc[hps, "bus0"] += " low voltage"

    rh = n.links.index[n.links.carrier.str.contains("resistive heater")]
    n.links.loc[rh, "bus0"] += " low voltage"

    mchp = n.links.index[n.links.carrier.str.contains("micro gas")]
    n.links.loc[mchp, "bus1"] += " low voltage"

    # set existing solar to cost of utility cost rather the 50-50 rooftop-utility
    solar = n.generators.index[n.generators.carrier == "solar"]
    n.generators.loc[solar, "capital_cost"] = costs.at['solar-utility', 'fixed']
    if snakemake.wildcards.clusters[-1:] == "m":
        simplified_pop_layout = pd.read_csv(snakemake.input.simplified_pop_layout, index_col=0)
        pop_solar = simplified_pop_layout.total.rename(index = lambda x: x + " solar")
    else:
        pop_solar = pop_layout.total.rename(index = lambda x: x + " solar")

    # add max solar rooftop potential assuming 0.1 kW/m2 and 10 m2/person,
    # i.e. 1 kW/person (population data is in thousands of people) so we get MW
    potential = 0.1 * 10 * pop_solar

    n.madd("Generator",
        solar,
        suffix=" rooftop",
        bus=n.generators.loc[solar, "bus"] + " low voltage",
        carrier="solar rooftop",
        p_nom_extendable=True,
        p_nom_max=potential,
        marginal_cost=n.generators.loc[solar, 'marginal_cost'],
        capital_cost=costs.at['solar-rooftop', 'fixed'],
        efficiency=n.generators.loc[solar, 'efficiency'],
        p_max_pu=n.generators_t.p_max_pu[solar],
        lifetime=costs.at['solar-rooftop', 'lifetime']
    )

    n.add("Carrier", "home battery")

    n.madd("Bus",
        nodes + " home battery",
        location=nodes,
        carrier="home battery"
    )

    n.madd("Store",
        nodes + " home battery",
        bus=nodes + " home battery",
        e_cyclic=True,
        e_nom_extendable=True,
        carrier="home battery",
        capital_cost=costs.at['home battery storage', 'fixed'],
        lifetime=costs.at['battery storage', 'lifetime']
    )

    n.madd("Link",
        nodes + " home battery charger",
        bus0=nodes + " low voltage",
        bus1=nodes + " home battery",
        carrier="home battery charger",
        efficiency=costs.at['battery inverter', 'efficiency']**0.5,
        capital_cost=costs.at['home battery inverter', 'fixed'],
        p_nom_extendable=True,
        lifetime=costs.at['battery inverter', 'lifetime']
    )

    n.madd("Link",
        nodes + " home battery discharger",
        bus0=nodes + " home battery",
        bus1=nodes + " low voltage",
        carrier="home battery discharger",
        efficiency=costs.at['battery inverter', 'efficiency']**0.5,
        marginal_cost=options['marginal_cost_storage'],
        p_nom_extendable=True,
        lifetime=costs.at['battery inverter', 'lifetime']
    )


def insert_gas_distribution_costs(n, costs):
    # TODO options?

    f_costs = options['gas_distribution_grid_cost_factor']

    print("Inserting gas distribution grid with investment cost factor of", f_costs)

    capital_cost = costs.loc['electricity distribution grid']["fixed"] * f_costs

    # gas boilers
    gas_b = n.links.index[n.links.carrier.str.contains("gas boiler") &
                          (~n.links.carrier.str.contains("urban central"))]
    n.links.loc[gas_b, "capital_cost"] += capital_cost

    # micro CHPs
    mchp = n.links.index[n.links.carrier.str.contains("micro gas")]
    n.links.loc[mchp,  "capital_cost"] += capital_cost


def add_electricity_grid_connection(n, costs):

    carriers = ["onwind", "solar"]

    gens = n.generators.index[n.generators.carrier.isin(carriers)]

    n.generators.loc[gens, "capital_cost"] += costs.at['electricity grid connection', 'fixed']


def add_storage_and_grids(n, costs):

    logger.info("Add hydrogen storage")

    nodes = pop_layout.index

    n.add("Carrier", "H2")

    n.madd("Bus",
        nodes + " H2",
        location=nodes,
        carrier="H2"
    )

    n.madd("Link",
        nodes + " H2 Electrolysis",
        bus1=nodes + " H2",
        bus0=nodes,
        p_nom_extendable=True,
        carrier="H2 Electrolysis",
        efficiency=costs.at["electrolysis", "efficiency"],
        capital_cost=costs.at["electrolysis", "fixed"],
        lifetime=costs.at['electrolysis', 'lifetime']
    )

    n.madd("Link",
        nodes + " H2 Fuel Cell",
        bus0=nodes + " H2",
        bus1=nodes,
        p_nom_extendable=True,
        carrier ="H2 Fuel Cell",
        efficiency=costs.at["fuel cell", "efficiency"],
        capital_cost=costs.at["fuel cell", "fixed"] * costs.at["fuel cell", "efficiency"], #NB: fixed cost is per MWel
        lifetime=costs.at['fuel cell', 'lifetime']
    )

    cavern_types = snakemake.config["sector"]["hydrogen_underground_storage_locations"]
    h2_caverns = pd.read_csv(snakemake.input.h2_cavern, index_col=0)

    if not h2_caverns.empty and options['hydrogen_underground_storage']:

        h2_caverns = h2_caverns[cavern_types].sum(axis=1)

        # only use sites with at least 2 TWh potential
        h2_caverns = h2_caverns[h2_caverns > 2]

        # convert TWh to MWh
        h2_caverns = h2_caverns * 1e6

        # clip at 1000 TWh for one location
        h2_caverns.clip(upper=1e9, inplace=True)

        logger.info("Add hydrogen underground storage")

        h2_capital_cost = costs.at["hydrogen storage underground", "fixed"]

        n.madd("Store",
            h2_caverns.index + " H2 Store",
            bus=h2_caverns.index + " H2",
            e_nom_extendable=True,
            e_nom_max=h2_caverns.values,
            e_cyclic=True,
            carrier="H2 Store",
            capital_cost=h2_capital_cost,
            lifetime=costs.at["hydrogen storage underground", "lifetime"]
        )

    # hydrogen stored overground (where not already underground)
    h2_capital_cost = costs.at["hydrogen storage tank incl. compressor", "fixed"]
    nodes_overground = h2_caverns.index.symmetric_difference(nodes)

    n.madd("Store",
        nodes_overground + " H2 Store",
        bus=nodes_overground + " H2",
        e_nom_extendable=True,
        e_cyclic=True,
        carrier="H2 Store",
        capital_cost=h2_capital_cost
    )

    if options["gas_network"] or options["H2_retrofit"]:

        fn = snakemake.input.clustered_gas_network
        gas_pipes = pd.read_csv(fn, index_col=0)

    if options["gas_network"]:

        logger.info("Add natural gas infrastructure, incl. LNG terminals, production and entry-points.")

        if options["H2_retrofit"]:
            gas_pipes["p_nom_max"] = gas_pipes.p_nom
            gas_pipes["p_nom_min"] = 0.
            # 0.1 EUR/MWkm/a to prefer decommissioning to address degeneracy
            gas_pipes["capital_cost"] = 0.1 * gas_pipes.length
        else:
            gas_pipes["p_nom_max"] = np.inf
            gas_pipes["p_nom_min"] = gas_pipes.p_nom
            gas_pipes["capital_cost"] = gas_pipes.length * costs.at['CH4 (g) pipeline', 'fixed']

        n.madd("Link",
            gas_pipes.index,
            bus0=gas_pipes.bus0 + " gas",
            bus1=gas_pipes.bus1 + " gas",
            p_min_pu=gas_pipes.p_min_pu,
            p_nom=gas_pipes.p_nom,
            p_nom_extendable=True,
            p_nom_max=gas_pipes.p_nom_max,
            p_nom_min=gas_pipes.p_nom_min,
            length=gas_pipes.length,
            capital_cost=gas_pipes.capital_cost,
            tags=gas_pipes.name,
            carrier="gas pipeline",
            lifetime=costs.at['CH4 (g) pipeline', 'lifetime']
        )

        # remove fossil generators where there is neither
        # production, LNG terminal, nor entry-point beyond system scope

        fn = snakemake.input.gas_input_nodes_simplified
        gas_input_nodes = pd.read_csv(fn, index_col=0)

        unique = gas_input_nodes.index.unique()
        gas_i = n.generators.carrier == 'gas'
        internal_i = ~n.generators.bus.map(n.buses.location).isin(unique)

        remove_i = n.generators[gas_i & internal_i].index
        n.generators.drop(remove_i, inplace=True)

        p_nom = gas_input_nodes.sum(axis=1).rename(lambda x: x + " gas")
        n.generators.loc[gas_i, "p_nom_extendable"] = False
        n.generators.loc[gas_i, "p_nom"] = p_nom

        # add candidates for new gas pipelines to achieve full connectivity

        G = nx.Graph()

        gas_buses = n.buses.loc[n.buses.carrier=='gas', 'location']
        G.add_nodes_from(np.unique(gas_buses.values))

        sel = gas_pipes.p_nom > 1500
        attrs = ["bus0", "bus1", "length"]
        G.add_weighted_edges_from(gas_pipes.loc[sel, attrs].values)

        # find all complement edges
        complement_edges = pd.DataFrame(complement(G).edges, columns=["bus0", "bus1"])
        complement_edges["length"] = complement_edges.apply(haversine, axis=1)

        # apply k_edge_augmentation weighted by length of complement edges
        k_edge = options.get("gas_network_connectivity_upgrade", 3)
        augmentation = list(k_edge_augmentation(G, k_edge, avail=complement_edges.values))

        if augmentation:

            new_gas_pipes = pd.DataFrame(augmentation, columns=["bus0", "bus1"])
            new_gas_pipes["length"] = new_gas_pipes.apply(haversine, axis=1)

            new_gas_pipes.index = new_gas_pipes.apply(
                lambda x: f"gas pipeline new {x.bus0} <-> {x.bus1}", axis=1)

            n.madd("Link",
                new_gas_pipes.index,
                bus0=new_gas_pipes.bus0 + " gas",
                bus1=new_gas_pipes.bus1 + " gas",
                p_min_pu=-1, # new gas pipes are bidirectional
                p_nom_extendable=True,
                length=new_gas_pipes.length,
                capital_cost=new_gas_pipes.length * costs.at['CH4 (g) pipeline', 'fixed'],
                carrier="gas pipeline new",
                lifetime=costs.at['CH4 (g) pipeline', 'lifetime']
            )

    if options["H2_retrofit"]:

        logger.info("Add retrofitting options of existing CH4 pipes to H2 pipes.")

        fr = "gas pipeline"
        to = "H2 pipeline retrofitted"
        h2_pipes = gas_pipes.rename(index=lambda x: x.replace(fr, to))

        n.madd("Link",
            h2_pipes.index,
            bus0=h2_pipes.bus0 + " H2",
            bus1=h2_pipes.bus1 + " H2",
            p_min_pu=-1.,  # allow that all H2 retrofit pipelines can be used in both directions
            p_nom_max=h2_pipes.p_nom * options["H2_retrofit_capacity_per_CH4"],
            p_nom_extendable=True,
            length=h2_pipes.length,
            capital_cost=costs.at['H2 (g) pipeline repurposed', 'fixed'] * h2_pipes.length,
            tags=h2_pipes.name,
            carrier="H2 pipeline retrofitted",
            lifetime=costs.at['H2 (g) pipeline repurposed', 'lifetime']
        )

    if options.get("H2_network", True):

        logger.info("Add options for new hydrogen pipelines.")

        h2_pipes = create_network_topology(n, "H2 pipeline ", carriers=["DC", "gas pipeline"])

        # TODO Add efficiency losses
        n.madd("Link",
            h2_pipes.index,
            bus0=h2_pipes.bus0.values + " H2",
            bus1=h2_pipes.bus1.values + " H2",
            p_min_pu=-1,
            p_nom_extendable=True,
            length=h2_pipes.length.values,
            capital_cost=costs.at['H2 (g) pipeline', 'fixed'] * h2_pipes.length.values,
            carrier="H2 pipeline",
            lifetime=costs.at['H2 (g) pipeline', 'lifetime']
        )

    n.add("Carrier", "battery")

    n.madd("Bus",
        nodes + " battery",
        location=nodes,
        carrier="battery"
    )

    n.madd("Store",
        nodes + " battery",
        bus=nodes + " battery",
        e_cyclic=True,
        e_nom_extendable=True,
        carrier="battery",
        capital_cost=costs.at['battery storage', 'fixed'],
        lifetime=costs.at['battery storage', 'lifetime']
    )

    n.madd("Link",
        nodes + " battery charger",
        bus0=nodes,
        bus1=nodes + " battery",
        carrier="battery charger",
        efficiency=costs.at['battery inverter', 'efficiency']**0.5,
        capital_cost=costs.at['battery inverter', 'fixed'],
        p_nom_extendable=True,
        lifetime=costs.at['battery inverter', 'lifetime']
    )

    n.madd("Link",
        nodes + " battery discharger",
        bus0=nodes + " battery",
        bus1=nodes,
        carrier="battery discharger",
        efficiency=costs.at['battery inverter', 'efficiency']**0.5,
        marginal_cost=options['marginal_cost_storage'],
        p_nom_extendable=True,
        lifetime=costs.at['battery inverter', 'lifetime']
    )

    if options['methanation']:

        n.madd("Link",
            spatial.nodes,
            suffix=" Sabatier",
            bus0=nodes + " H2",
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            p_nom_extendable=True,
            carrier="Sabatier",
            efficiency=costs.at["methanation", "efficiency"],
            efficiency2=-costs.at["methanation", "efficiency"] * costs.at['gas', 'CO2 intensity'],
            capital_cost=costs.at["methanation", "fixed"] * costs.at["methanation", "efficiency"],  # costs given per kW_gas
            lifetime=costs.at['methanation', 'lifetime']
        )

    if options['helmeth']:

        n.madd("Link",
            spatial.nodes,
            suffix=" helmeth",
            bus0=nodes,
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            carrier="helmeth",
            p_nom_extendable=True,
            efficiency=costs.at["helmeth", "efficiency"],
            efficiency2=-costs.at["helmeth", "efficiency"] * costs.at['gas', 'CO2 intensity'],
            capital_cost=costs.at["helmeth", "fixed"],
            lifetime=costs.at['helmeth', 'lifetime']
        )


    if options['SMR']:

        n.madd("Link",
            spatial.nodes,
            suffix=" SMR CC",
            bus0=spatial.gas.nodes,
            bus1=nodes + " H2",
            bus2="co2 atmosphere",
            bus3=spatial.co2.nodes,
            p_nom_extendable=True,
            carrier="SMR CC",
            efficiency=costs.at["SMR CC", "efficiency"],
            efficiency2=costs.at['gas', 'CO2 intensity'] * (1 - options["cc_fraction"]),
            efficiency3=costs.at['gas', 'CO2 intensity'] * options["cc_fraction"],
            capital_cost=costs.at["SMR CC", "fixed"],
            lifetime=costs.at['SMR CC', 'lifetime']
        )

        n.madd("Link",
            nodes + " SMR",
            bus0=spatial.gas.nodes,
            bus1=nodes + " H2",
            bus2="co2 atmosphere",
            p_nom_extendable=True,
            carrier="SMR",
            efficiency=costs.at["SMR", "efficiency"],
            efficiency2=costs.at['gas', 'CO2 intensity'],
            capital_cost=costs.at["SMR", "fixed"],
            lifetime=costs.at['SMR', 'lifetime']
        )


def add_land_transport(n, costs):
    # TODO options?

    logger.info("Add land transport")

    transport = pd.read_csv(snakemake.input.transport_demand, index_col=0, parse_dates=True)
    number_cars = pd.read_csv(snakemake.input.transport_data, index_col=0)["number cars"]
    avail_profile = pd.read_csv(snakemake.input.avail_profile, index_col=0, parse_dates=True)
    dsm_profile = pd.read_csv(snakemake.input.dsm_profile, index_col=0, parse_dates=True)

    fuel_cell_share = get(options["land_transport_fuel_cell_share"], investment_year)
    electric_share = get(options["land_transport_electric_share"], investment_year)
    ice_share = 1 - fuel_cell_share - electric_share

    print("FCEV share", fuel_cell_share)
    print("EV share", electric_share)
    print("ICEV share", ice_share)

    assert ice_share >= 0, "Error, more FCEV and EV share than 1."

    nodes = pop_layout.index

    if electric_share > 0:

        n.add("Carrier", "Li ion")

        n.madd("Bus",
            nodes,
            location=nodes,
            suffix=" EV battery",
            carrier="Li ion"
        )

        p_set = electric_share * (transport[nodes] + cycling_shift(transport[nodes], 1) + cycling_shift(transport[nodes], 2)) / 3

        n.madd("Load",
            nodes,
            suffix=" land transport EV",
            bus=nodes + " EV battery",
            carrier="land transport EV",
            p_set=p_set
        )

        p_nom = number_cars * options.get("bev_charge_rate", 0.011) * electric_share

        n.madd("Link",
            nodes,
            suffix= " BEV charger",
            bus0=nodes,
            bus1=nodes + " EV battery",
            p_nom=p_nom,
            carrier="BEV charger",
            p_max_pu=avail_profile[nodes],
            efficiency=options.get("bev_charge_efficiency", 0.9),
            #These were set non-zero to find LU infeasibility when availability = 0.25
            #p_nom_extendable=True,
            #p_nom_min=p_nom,
            #capital_cost=1e6,  #i.e. so high it only gets built where necessary
        )

    if electric_share > 0 and options["v2g"]:

        n.madd("Link",
            nodes,
            suffix=" V2G",
            bus1=nodes,
            bus0=nodes + " EV battery",
            p_nom=p_nom,
            carrier="V2G",
            p_max_pu=avail_profile[nodes],
            efficiency=options.get("bev_charge_efficiency", 0.9),
        )

    if electric_share > 0 and options["bev_dsm"]:

        e_nom = number_cars * options.get("bev_energy", 0.05) * options["bev_availability"] * electric_share

        n.madd("Store",
            nodes,
            suffix=" battery storage",
            bus=nodes + " EV battery",
            carrier="battery storage",
            e_cyclic=True,
            e_nom=e_nom,
            e_max_pu=1,
            e_min_pu=dsm_profile[nodes]
        )

    if fuel_cell_share > 0:

        n.madd("Load",
            nodes,
            suffix=" land transport fuel cell",
            bus=nodes + " H2",
            carrier="land transport fuel cell",
            p_set=fuel_cell_share / options['transport_fuel_cell_efficiency'] * transport[nodes]
        )

    if ice_share > 0:

        if "oil" not in n.buses.carrier.unique():
            n.madd("Bus",
                spatial.oil.nodes,
                location=spatial.oil.locations,
                carrier="oil"
            )

        ice_efficiency = options['transport_internal_combustion_efficiency']

        n.madd("Load",
            nodes,
            suffix=" land transport oil",
            bus=spatial.oil.nodes,
            carrier="land transport oil",
            p_set=ice_share / ice_efficiency * transport[nodes]
        )

        co2 = ice_share / ice_efficiency * transport[nodes].sum().sum() / 8760 * costs.at["oil", 'CO2 intensity']

        n.add("Load",
            "land transport oil emissions",
            bus="co2 atmosphere",
            carrier="land transport oil emissions",
            p_set=-co2
        )


def build_heat_demand(n):

    # copy forward the daily average heat demand into each hour, so it can be multipled by the intraday profile
    daily_space_heat_demand = xr.open_dataarray(snakemake.input.heat_demand_total).to_pandas().reindex(index=n.snapshots, method="ffill")

    intraday_profiles = pd.read_csv(snakemake.input.heat_profile, index_col=0)

    sectors = ["residential", "services"]
    uses = ["water", "space"]

    heat_demand = {}
    electric_heat_supply = {}
    for sector, use in product(sectors, uses):
        weekday = list(intraday_profiles[f"{sector} {use} weekday"])
        weekend = list(intraday_profiles[f"{sector} {use} weekend"])
        weekly_profile = weekday * 5 + weekend * 2
        intraday_year_profile = generate_periodic_profiles(
            daily_space_heat_demand.index.tz_localize("UTC"),
            nodes=daily_space_heat_demand.columns,
            weekly_profile=weekly_profile
        )

        if use == "space":
            heat_demand_shape = daily_space_heat_demand * intraday_year_profile
        else:
            heat_demand_shape = intraday_year_profile

        heat_demand[f"{sector} {use}"] = (heat_demand_shape/heat_demand_shape.sum()).multiply(pop_weighted_energy_totals[f"total {sector} {use}"]) * 1e6
        electric_heat_supply[f"{sector} {use}"] = (heat_demand_shape/heat_demand_shape.sum()).multiply(pop_weighted_energy_totals[f"electricity {sector} {use}"]) * 1e6

    heat_demand = pd.concat(heat_demand, axis=1)
    electric_heat_supply = pd.concat(electric_heat_supply, axis=1)

    # subtract from electricity load since heat demand already in heat_demand
    electric_nodes = n.loads.index[n.loads.carrier == "electricity"]
    n.loads_t.p_set[electric_nodes] = n.loads_t.p_set[electric_nodes] - electric_heat_supply.groupby(level=1, axis=1).sum()[electric_nodes]

    return heat_demand


def add_heat(n, costs):

    logger.info("Add heat sector")

    sectors = ["residential", "services"]

    heat_demand = build_heat_demand(n)

    nodes, dist_fraction, urban_fraction = create_nodes_for_heat_sector()

    #NB: must add costs of central heating afterwards (EUR 400 / kWpeak, 50a, 1% FOM from Fraunhofer ISE)

    # exogenously reduce space heat demand
    if options["reduce_space_heat_exogenously"]:
        dE = get(options["reduce_space_heat_exogenously_factor"], investment_year)
        print(f"assumed space heat reduction of {dE*100} %")
        for sector in sectors:
            heat_demand[sector + " space"] = (1 - dE) * heat_demand[sector + " space"]

    heat_systems = [
        "residential rural",
        "services rural",
        "residential urban decentral",
        "services urban decentral",
        "urban central"
    ]

    cop = {
        "air": xr.open_dataarray(snakemake.input.cop_air_total).to_pandas().reindex(index=n.snapshots),
        "ground": xr.open_dataarray(snakemake.input.cop_soil_total).to_pandas().reindex(index=n.snapshots)
    }

    solar_thermal = xr.open_dataarray(snakemake.input.solar_thermal_total).to_pandas().reindex(index=n.snapshots)
    # 1e3 converts from W/m^2 to MW/(1000m^2) = kW/m^2
    solar_thermal = options['solar_cf_correction'] * solar_thermal / 1e3

    for name in heat_systems:

        name_type = "central" if name == "urban central" else "decentral"

        n.add("Carrier", name + " heat")

        n.madd("Bus",
            nodes[name] + f" {name} heat",
            location=nodes[name],
            carrier=name + " heat"
        )

        ## Add heat load

        for sector in sectors:
            # heat demand weighting
            if "rural" in name:
                factor = 1 - urban_fraction[nodes[name]]
            elif "urban central" in name:
                factor = dist_fraction[nodes[name]]
            elif "urban decentral" in name:
                factor = urban_fraction[nodes[name]] - \
                    dist_fraction[nodes[name]]
            else:
                raise NotImplementedError(f" {name} not in " f"heat systems: {heat_systems}")

            if sector in name:
                heat_load = heat_demand[[sector + " water",sector + " space"]].groupby(level=1,axis=1).sum()[nodes[name]].multiply(factor)

        if name == "urban central":
            heat_load = heat_demand.groupby(level=1,axis=1).sum()[nodes[name]].multiply(factor * (1 + options['district_heating']['district_heating_loss']))

        n.madd("Load",
            nodes[name],
            suffix=f" {name} heat",
            bus=nodes[name] + f" {name} heat",
            carrier=name + " heat",
            p_set=heat_load
        )

        ## Add heat pumps

        heat_pump_type = "air" if "urban" in name else "ground"

        costs_name = f"{name_type} {heat_pump_type}-sourced heat pump"
        efficiency = cop[heat_pump_type][nodes[name]] if options["time_dep_hp_cop"] else costs.at[costs_name, 'efficiency']

        n.madd("Link",
            nodes[name],
            suffix=f" {name} {heat_pump_type} heat pump",
            bus0=nodes[name],
            bus1=nodes[name] + f" {name} heat",
            carrier=f"{name} {heat_pump_type} heat pump",
            efficiency=efficiency,
            capital_cost=costs.at[costs_name, 'efficiency'] * costs.at[costs_name, 'fixed'],
            p_nom_extendable=True,
            lifetime=costs.at[costs_name, 'lifetime']
        )

        if options["tes"]:

            n.add("Carrier", name + " water tanks")

            n.madd("Bus",
                nodes[name] + f" {name} water tanks",
                location=nodes[name],
                carrier=name + " water tanks"
            )

            n.madd("Link",
                nodes[name] + f" {name} water tanks charger",
                bus0=nodes[name] + f" {name} heat",
                bus1=nodes[name] + f" {name} water tanks",
                efficiency=costs.at['water tank charger', 'efficiency'],
                carrier=name + " water tanks charger",
                p_nom_extendable=True
            )

            n.madd("Link",
                nodes[name] + f" {name} water tanks discharger",
                bus0=nodes[name] + f" {name} water tanks",
                bus1=nodes[name] + f" {name} heat",
                carrier=name + " water tanks discharger",
                efficiency=costs.at['water tank discharger', 'efficiency'],
                p_nom_extendable=True
            )


            if isinstance(options["tes_tau"], dict):
                tes_time_constant_days = options["tes_tau"][name_type]
            else:
                logger.warning("Deprecated: a future version will require you to specify 'tes_tau' ",
                               "for 'decentral' and 'central' separately.")
                tes_time_constant_days = options["tes_tau"] if name_type == "decentral" else 180.

            # conversion from EUR/m^3 to EUR/MWh for 40 K diff and 1.17 kWh/m^3/K
            capital_cost = costs.at[name_type + ' water tank storage', 'fixed'] / 0.00117 / 40

            n.madd("Store",
                nodes[name] + f" {name} water tanks",
                bus=nodes[name] + f" {name} water tanks",
                e_cyclic=True,
                e_nom_extendable=True,
                carrier=name + " water tanks",
                standing_loss=1 - np.exp(- 1 / 24 / tes_time_constant_days),
                capital_cost=capital_cost,
                lifetime=costs.at[name_type + ' water tank storage', 'lifetime']
            )

        if options["boilers"]:

            key = f"{name_type} resistive heater"

            n.madd("Link",
                nodes[name] + f" {name} resistive heater",
                bus0=nodes[name],
                bus1=nodes[name] + f" {name} heat",
                carrier=name + " resistive heater",
                efficiency=costs.at[key, 'efficiency'],
                capital_cost=costs.at[key, 'efficiency'] * costs.at[key, 'fixed'],
                p_nom_extendable=True,
                lifetime=costs.at[key, 'lifetime']
            )

            key = f"{name_type} gas boiler"

            n.madd("Link",
                nodes[name] + f" {name} gas boiler",
                p_nom_extendable=True,
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name] + f" {name} heat",
                bus2="co2 atmosphere",
                carrier=name + " gas boiler",
                efficiency=costs.at[key, 'efficiency'],
                efficiency2=costs.at['gas', 'CO2 intensity'],
                capital_cost=costs.at[key, 'efficiency'] * costs.at[key, 'fixed'],
                lifetime=costs.at[key, 'lifetime']
            )

        if options["solar_thermal"]:

            n.add("Carrier", name + " solar thermal")

            n.madd("Generator",
                nodes[name],
                suffix=f" {name} solar thermal collector",
                bus=nodes[name] + f" {name} heat",
                carrier=name + " solar thermal",
                p_nom_extendable=True,
                capital_cost=costs.at[name_type + ' solar thermal', 'fixed'],
                p_max_pu=solar_thermal[nodes[name]],
                lifetime=costs.at[name_type + ' solar thermal', 'lifetime']
            )

        if options["chp"] and name == "urban central":

            # add gas CHP; biomass CHP is added in biomass section
            n.madd("Link",
                nodes[name] + " urban central gas CHP",
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name],
                bus2=nodes[name] + " urban central heat",
                bus3="co2 atmosphere",
                carrier="urban central gas CHP",
                p_nom_extendable=True,
                capital_cost=costs.at['central gas CHP', 'fixed'] * costs.at['central gas CHP', 'efficiency'],
                marginal_cost=costs.at['central gas CHP', 'VOM'],
                efficiency=costs.at['central gas CHP', 'efficiency'],
                efficiency2=costs.at['central gas CHP', 'efficiency'] / costs.at['central gas CHP', 'c_b'],
                efficiency3=costs.at['gas', 'CO2 intensity'],
                lifetime=costs.at['central gas CHP', 'lifetime']
            )

            n.madd("Link",
                nodes[name] + " urban central gas CHP CC",
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name],
                bus2=nodes[name] + " urban central heat",
                bus3="co2 atmosphere",
                bus4=spatial.co2.df.loc[nodes[name], "nodes"].values,
                carrier="urban central gas CHP CC",
                p_nom_extendable=True,
                capital_cost=costs.at['central gas CHP', 'fixed']*costs.at['central gas CHP', 'efficiency'] + costs.at['biomass CHP capture', 'fixed']*costs.at['gas', 'CO2 intensity'],
                marginal_cost=costs.at['central gas CHP', 'VOM'],
                efficiency=costs.at['central gas CHP', 'efficiency'] - costs.at['gas', 'CO2 intensity'] * (costs.at['biomass CHP capture', 'electricity-input'] + costs.at['biomass CHP capture', 'compression-electricity-input']),
                efficiency2=costs.at['central gas CHP', 'efficiency'] / costs.at['central gas CHP', 'c_b'] + costs.at['gas', 'CO2 intensity'] * (costs.at['biomass CHP capture', 'heat-output'] + costs.at['biomass CHP capture', 'compression-heat-output'] - costs.at['biomass CHP capture', 'heat-input']),
                efficiency3=costs.at['gas', 'CO2 intensity'] * (1-costs.at['biomass CHP capture', 'capture_rate']),
                efficiency4=costs.at['gas', 'CO2 intensity'] * costs.at['biomass CHP capture', 'capture_rate'],
                lifetime=costs.at['central gas CHP', 'lifetime']
            )

        if options["chp"] and options["micro_chp"] and name != "urban central":

            n.madd("Link",
                nodes[name] + f" {name} micro gas CHP",
                p_nom_extendable=True,
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name],
                bus2=nodes[name] + f" {name} heat",
                bus3="co2 atmosphere",
                carrier=name + " micro gas CHP",
                efficiency=costs.at['micro CHP', 'efficiency'],
                efficiency2=costs.at['micro CHP', 'efficiency-heat'],
                efficiency3=costs.at['gas', 'CO2 intensity'],
                capital_cost=costs.at['micro CHP', 'fixed'],
                lifetime=costs.at['micro CHP', 'lifetime']
            )


    if options['retrofitting']['retro_endogen']:

        logger.info("Add retrofitting endogenously")

        # resample heat demand temporal 'heat_demand_r' depending on in config
        # specified temporal resolution, to not overestimate retrofitting
        hours = list(filter(re.compile(r'^\d+h$', re.IGNORECASE).search, opts))
        if len(hours)==0:
            hours = [n.snapshots[1] - n.snapshots[0]]
        heat_demand_r =  heat_demand.resample(hours[0]).mean()

        # retrofitting data 'retro_data' with 'costs' [EUR/m^2] and heat
        # demand 'dE' [per unit of original heat demand] for each country and
        # different retrofitting strengths [additional insulation thickness in m]
        retro_data = pd.read_csv(snakemake.input.retro_cost_energy,
                                 index_col=[0, 1], skipinitialspace=True,
                                 header=[0, 1])
        # heated floor area [10^6 * m^2] per country
        floor_area = pd.read_csv(snakemake.input.floor_area, index_col=[0, 1])

        n.add("Carrier", "retrofitting")

        # share of space heat demand 'w_space' of total heat demand
        w_space = {}
        for sector in sectors:
            w_space[sector] = heat_demand_r[sector + " space"] / \
                (heat_demand_r[sector + " space"] + heat_demand_r[sector + " water"])
        w_space["tot"] = ((heat_demand_r["services space"] +
                           heat_demand_r["residential space"]) /
                           heat_demand_r.groupby(level=[1], axis=1).sum())


        for name in n.loads[n.loads.carrier.isin([x + " heat" for x in heat_systems])].index:

            node = n.buses.loc[name, "location"]
            ct = pop_layout.loc[node, "ct"]

            # weighting 'f' depending on the size of the population at the node
            f = urban_fraction[node] if "urban" in name else (1-urban_fraction[node])
            if f == 0:
                continue
            # get sector name ("residential"/"services"/or both "tot" for urban central)
            sec = [x if x in name else "tot" for x in sectors][0]

            # get floor aread at node and region (urban/rural) in m^2
            floor_area_node = ((pop_layout.loc[node].fraction
                                  * floor_area.loc[ct, "value"] * 10**6).loc[sec] * f)
            # total heat demand at node [MWh]
            demand = (n.loads_t.p_set[name].resample(hours[0])
                      .mean())

            # space heat demand at node [MWh]
            space_heat_demand = demand * w_space[sec][node]
            # normed time profile of space heat demand 'space_pu' (values between 0-1),
            # p_max_pu/p_min_pu of retrofitting generators
            space_pu = (space_heat_demand / space_heat_demand.max()).to_frame(name=node)

            # minimum heat demand 'dE' after retrofitting in units of original heat demand (values between 0-1)
            dE = retro_data.loc[(ct, sec), ("dE")]
            # get addtional energy savings 'dE_diff' between the different retrofitting strengths/generators at one node
            dE_diff = abs(dE.diff()).fillna(1-dE.iloc[0])
            # convert costs Euro/m^2 -> Euro/MWh
            capital_cost =  retro_data.loc[(ct, sec), ("cost")] * floor_area_node / \
                            ((1 - dE) * space_heat_demand.max())
            # number of possible retrofitting measures 'strengths' (set in list at config.yaml 'l_strength')
            # given in additional insulation thickness [m]
            # for each measure, a retrofitting generator is added at the node
            strengths = retro_data.columns.levels[1]

            # check that ambitious retrofitting has higher costs per MWh than moderate retrofitting
            if (capital_cost.diff() < 0).sum():
                logger.warning(f"Costs are not linear for {ct} {sec}")
                s = capital_cost[(capital_cost.diff() < 0)].index
                strengths = strengths.drop(s)

            # reindex normed time profile of space heat demand back to hourly resolution
            space_pu = space_pu.reindex(index=heat_demand.index).fillna(method="ffill")

            # add for each retrofitting strength a generator with heat generation profile following the profile of the heat demand
            for strength in strengths:
                n.madd('Generator',
                    [node],
                    suffix=' retrofitting ' + strength + " " + name[6::],
                    bus=name,
                    carrier="retrofitting",
                    p_nom_extendable=True,
                    p_nom_max=dE_diff[strength] * space_heat_demand.max(), # maximum energy savings for this renovation strength
                    p_max_pu=space_pu,
                    p_min_pu=space_pu,
                    country=ct,
                    capital_cost=capital_cost[strength] * options['retrofitting']['cost_factor']
                )


def create_nodes_for_heat_sector():
    # TODO pop_layout

    # rural are areas with low heating density and individual heating
    # urban are areas with high heating density
    # urban can be split into district heating (central) and individual heating (decentral)

    ct_urban = pop_layout.urban.groupby(pop_layout.ct).sum()
    # distribution of urban population within a country
    pop_layout["urban_ct_fraction"] = pop_layout.urban / pop_layout.ct.map(ct_urban.get)

    sectors = ["residential", "services"]

    nodes = {}
    urban_fraction = pop_layout.urban / pop_layout[["rural", "urban"]].sum(axis=1)

    for sector in sectors:
        nodes[sector + " rural"] = pop_layout.index
        nodes[sector + " urban decentral"] = pop_layout.index

    district_heat_share = pop_weighted_energy_totals["district heat share"]

    # maximum potential of urban demand covered by district heating
    central_fraction = options['district_heating']["potential"]
    # district heating share at each node
    dist_fraction_node = district_heat_share * pop_layout["urban_ct_fraction"] / pop_layout["fraction"]
    nodes["urban central"] = dist_fraction_node.index
    # if district heating share larger than urban fraction -> set urban
    # fraction to district heating share
    urban_fraction = pd.concat([urban_fraction, dist_fraction_node],
                               axis=1).max(axis=1)
    # difference of max potential and today's share of district heating
    diff = (urban_fraction * central_fraction) - dist_fraction_node
    progress = get(options["district_heating"]["progress"], investment_year)
    dist_fraction_node += diff * progress
    print(
        "The current district heating share compared to the maximum",
        f"possible is increased by a progress factor of\n{progress}",
        f"resulting in a district heating share of\n{dist_fraction_node}"
    )

    return nodes, dist_fraction_node, urban_fraction


def add_biomass(n, costs):

    logger.info("Add biomass")

    biomass_potentials = pd.read_csv(snakemake.input.biomass_potentials, index_col=0)

    # need to aggregate potentials if gas not nodally resolved
    if options["gas_network"]:
        biogas_potentials_spatial = biomass_potentials["biogas"].rename(index=lambda x: x + " biogas")
    else:
        biogas_potentials_spatial = biomass_potentials["biogas"].sum()

    if options["biomass_transport"]:
        solid_biomass_potentials_spatial = biomass_potentials["solid biomass"].rename(index=lambda x: x + " solid biomass")
    else:
        solid_biomass_potentials_spatial = biomass_potentials["solid biomass"].sum()


    n.add("Carrier", "biogas")
    n.add("Carrier", "solid biomass")

    n.madd("Bus",
        spatial.gas.biogas,
        location=spatial.gas.locations,
        carrier="biogas"
    )

    n.madd("Bus",
        spatial.biomass.nodes,
        location=spatial.biomass.locations,
        carrier="solid biomass"
    )

    n.madd("Store",
        spatial.gas.biogas,
        bus=spatial.gas.biogas,
        carrier="biogas",
        e_nom=biogas_potentials_spatial,
        marginal_cost=costs.at['biogas', 'fuel'],
        e_initial=biogas_potentials_spatial
    )

    n.madd("Store",
        spatial.biomass.nodes,
        bus=spatial.biomass.nodes,
        carrier="solid biomass",
        e_nom=solid_biomass_potentials_spatial,
        marginal_cost=costs.at['solid biomass', 'fuel'],
        e_initial=solid_biomass_potentials_spatial
    )

    n.madd("Link",
        spatial.gas.biogas_to_gas,
        bus0=spatial.gas.biogas,
        bus1=spatial.gas.nodes,
        bus2="co2 atmosphere",
        carrier="biogas to gas",
        capital_cost=costs.loc["biogas upgrading", "fixed"],
        marginal_cost=costs.loc["biogas upgrading", "VOM"],
        efficiency2=-costs.at['gas', 'CO2 intensity'],
        p_nom_extendable=True
    )

    if options["biomass_transport"]:

        transport_costs = pd.read_csv(
            snakemake.input.biomass_transport_costs,
            index_col=0,
        ).squeeze()

        # add biomass transport
        biomass_transport = create_network_topology(n, "biomass transport ", bidirectional=False)

        # costs
        bus0_costs = biomass_transport.bus0.apply(lambda x: transport_costs[x[:2]])
        bus1_costs = biomass_transport.bus1.apply(lambda x: transport_costs[x[:2]])
        biomass_transport["costs"] = pd.concat([bus0_costs, bus1_costs], axis=1).mean(axis=1)

        n.madd("Link",
            biomass_transport.index,
            bus0=biomass_transport.bus0 + " solid biomass",
            bus1=biomass_transport.bus1 + " solid biomass",
            p_nom_extendable=True,
            length=biomass_transport.length.values,
            marginal_cost=biomass_transport.costs * biomass_transport.length.values,
            capital_cost=1,
            carrier="solid biomass transport"
        )

    #AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if not urban_central.empty and options["chp"]:
        urban_central = urban_central.str[:-len(" urban central heat")]

        key = 'central solid biomass CHP'

        n.madd("Link",
            urban_central + " urban central solid biomass CHP",
            bus0=spatial.biomass.df.loc[urban_central, "nodes"].values,
            bus1=urban_central,
            bus2=urban_central + " urban central heat",
            carrier="urban central solid biomass CHP",
            p_nom_extendable=True,
            capital_cost=costs.at[key, 'fixed'] * costs.at[key, 'efficiency'],
            marginal_cost=costs.at[key, 'VOM'],
            efficiency=costs.at[key, 'efficiency'],
            efficiency2=costs.at[key, 'efficiency-heat'],
            lifetime=costs.at[key, 'lifetime']
        )

        n.madd("Link",
            urban_central + " urban central solid biomass CHP CC",
            bus0=spatial.biomass.df.loc[urban_central, "nodes"].values,
            bus1=urban_central,
            bus2=urban_central + " urban central heat",
            bus3="co2 atmosphere",
            bus4=spatial.co2.df.loc[urban_central, "nodes"].values,
            carrier="urban central solid biomass CHP CC",
            p_nom_extendable=True,
            capital_cost=costs.at[key, 'fixed'] * costs.at[key, 'efficiency'] + costs.at['biomass CHP capture', 'fixed'] * costs.at['solid biomass', 'CO2 intensity'],
            marginal_cost=costs.at[key, 'VOM'],
            efficiency=costs.at[key, 'efficiency'] - costs.at['solid biomass', 'CO2 intensity'] * (costs.at['biomass CHP capture', 'electricity-input'] + costs.at['biomass CHP capture', 'compression-electricity-input']),
            efficiency2=costs.at[key, 'efficiency-heat'] + costs.at['solid biomass', 'CO2 intensity'] * (costs.at['biomass CHP capture', 'heat-output'] + costs.at['biomass CHP capture', 'compression-heat-output'] - costs.at['biomass CHP capture', 'heat-input']),
            efficiency3=-costs.at['solid biomass', 'CO2 intensity'] * costs.at['biomass CHP capture', 'capture_rate'],
            efficiency4=costs.at['solid biomass', 'CO2 intensity'] * costs.at['biomass CHP capture', 'capture_rate'],
            lifetime=costs.at[key, 'lifetime']
        )


def add_industry(n, costs):

    logger.info("Add industrial demand")

    nodes = pop_layout.index

    # 1e6 to convert TWh to MWh
    industrial_demand = pd.read_csv(snakemake.input.industrial_demand, index_col=0) * 1e6

    n.madd("Bus",
        spatial.biomass.industry,
        location=spatial.biomass.locations,
        carrier="solid biomass for industry"
    )

    if options["biomass_transport"]:
        p_set = industrial_demand.loc[spatial.biomass.locations, "solid biomass"].rename(index=lambda x: x + " solid biomass for industry") / 8760
    else:
        p_set = industrial_demand["solid biomass"].sum() / 8760

    n.madd("Load",
        spatial.biomass.industry,
        bus=spatial.biomass.industry,
        carrier="solid biomass for industry",
        p_set=p_set
    )

    n.madd("Link",
        spatial.biomass.industry,
        bus0=spatial.biomass.nodes,
        bus1=spatial.biomass.industry,
        carrier="solid biomass for industry",
        p_nom_extendable=True,
        efficiency=1.
    )

    n.madd("Link",
        spatial.biomass.industry_cc,
        bus0=spatial.biomass.nodes,
        bus1=spatial.biomass.industry,
        bus2="co2 atmosphere",
        bus3=spatial.co2.nodes,
        carrier="solid biomass for industry CC",
        p_nom_extendable=True,
        capital_cost=costs.at["cement capture", "fixed"] * costs.at['solid biomass', 'CO2 intensity'],
        efficiency=0.9, # TODO: make config option
        efficiency2=-costs.at['solid biomass', 'CO2 intensity'] * costs.at["cement capture", "capture_rate"],
        efficiency3=costs.at['solid biomass', 'CO2 intensity'] * costs.at["cement capture", "capture_rate"],
        lifetime=costs.at['cement capture', 'lifetime']
    )

    n.madd("Bus",
        spatial.gas.industry,
        location=spatial.gas.locations,
        carrier="gas for industry")

    gas_demand = industrial_demand.loc[nodes, "methane"] / 8760.

    if options["gas_network"]:
        spatial_gas_demand = gas_demand.rename(index=lambda x: x + " gas for industry")
    else:
        spatial_gas_demand = gas_demand.sum()

    n.madd("Load",
        spatial.gas.industry,
        bus=spatial.gas.industry,
        carrier="gas for industry",
        p_set=spatial_gas_demand
    )

    n.madd("Link",
        spatial.gas.industry,
        bus0=spatial.gas.nodes,
        bus1=spatial.gas.industry,
        bus2="co2 atmosphere",
        carrier="gas for industry",
        p_nom_extendable=True,
        efficiency=1.,
        efficiency2=costs.at['gas', 'CO2 intensity']
    )

    n.madd("Link",
        spatial.gas.industry_cc,
        bus0=spatial.gas.nodes,
        bus1=spatial.gas.industry,
        bus2="co2 atmosphere",
        bus3=spatial.co2.nodes,
        carrier="gas for industry CC",
        p_nom_extendable=True,
        capital_cost=costs.at["cement capture", "fixed"] * costs.at['gas', 'CO2 intensity'],
        efficiency=0.9,
        efficiency2=costs.at['gas', 'CO2 intensity'] * (1 - costs.at["cement capture", "capture_rate"]),
        efficiency3=costs.at['gas', 'CO2 intensity'] * costs.at["cement capture", "capture_rate"],
        lifetime=costs.at['cement capture', 'lifetime']
    )

    n.madd("Load",
        nodes,
        suffix=" H2 for industry",
        bus=nodes + " H2",
        carrier="H2 for industry",
        p_set=industrial_demand.loc[nodes, "hydrogen"] / 8760
    )

    if options["shipping_hydrogen_liquefaction"]:

        n.madd("Bus",
            nodes,
            suffix=" H2 liquid",
            carrier="H2 liquid",
            location=nodes
        )

        n.madd("Link",
            nodes + " H2 liquefaction",
            bus0=nodes + " H2",
            bus1=nodes + " H2 liquid",
            carrier="H2 liquefaction",
            efficiency=costs.at["H2 liquefaction", 'efficiency'],
            capital_cost=costs.at["H2 liquefaction", 'fixed'],
            p_nom_extendable=True,
            lifetime=costs.at['H2 liquefaction', 'lifetime']
        )

        shipping_bus = nodes + " H2 liquid"
    else:
        shipping_bus = nodes + " H2"

    all_navigation = ["total international navigation", "total domestic navigation"]
    efficiency = options['shipping_average_efficiency'] / costs.at["fuel cell", "efficiency"]
    shipping_hydrogen_share = get(options['shipping_hydrogen_share'], investment_year)
    p_set = shipping_hydrogen_share * pop_weighted_energy_totals.loc[nodes, all_navigation].sum(axis=1) * 1e6 * efficiency / 8760

    n.madd("Load",
        nodes,
        suffix=" H2 for shipping",
        bus=shipping_bus,
        carrier="H2 for shipping",
        p_set=p_set
    )

    if shipping_hydrogen_share < 1:

        shipping_oil_share = 1 - shipping_hydrogen_share

        p_set = shipping_oil_share * pop_weighted_energy_totals.loc[nodes, all_navigation].sum(axis=1) * 1e6 / 8760.

        n.madd("Load",
            nodes,
            suffix=" shipping oil",
            bus=spatial.oil.nodes,
            carrier="shipping oil",
            p_set=p_set
        )

        co2 = shipping_oil_share * pop_weighted_energy_totals.loc[nodes, all_navigation].sum().sum() * 1e6 / 8760 * costs.at["oil", "CO2 intensity"]

        n.add("Load",
            "shipping oil emissions",
            bus="co2 atmosphere",
            carrier="shipping oil emissions",
            p_set=-co2
        )

    if "oil" not in n.buses.carrier.unique():
        n.madd("Bus",
            spatial.oil.nodes,
            location=spatial.oil.locations,
            carrier="oil"
        )

    if "oil" not in n.stores.carrier.unique():

        #could correct to e.g. 0.001 EUR/kWh * annuity and O&M
        n.madd("Store",
            [oil_bus + " Store" for oil_bus in spatial.oil.nodes],
            bus=spatial.oil.nodes,
            e_nom_extendable=True,
            e_cyclic=True,
            carrier="oil",
        )

    if "oil" not in n.generators.carrier.unique():

        n.madd("Generator",
            spatial.oil.nodes,
            bus=spatial.oil.nodes,
            p_nom_extendable=True,
            carrier="oil",
            marginal_cost=costs.at["oil", 'fuel']
        )

    if options["oil_boilers"]:

        nodes_heat = create_nodes_for_heat_sector()[0]

        for name in ["residential rural", "services rural", "residential urban decentral", "services urban decentral"]:

            n.madd("Link",
                nodes_heat[name] + f" {name} oil boiler",
                p_nom_extendable=True,
                bus0=spatial.oil.nodes,
                bus1=nodes_heat[name] + f" {name}  heat",
                bus2="co2 atmosphere",
                carrier=f"{name} oil boiler",
                efficiency=costs.at['decentral oil boiler', 'efficiency'],
                efficiency2=costs.at['oil', 'CO2 intensity'],
                capital_cost=costs.at['decentral oil boiler', 'efficiency'] * costs.at['decentral oil boiler', 'fixed'],
                lifetime=costs.at['decentral oil boiler', 'lifetime']
            )

    n.madd("Link",
        nodes + " Fischer-Tropsch",
        bus0=nodes + " H2",
        bus1=spatial.oil.nodes,
        bus2=spatial.co2.nodes,
        carrier="Fischer-Tropsch",
        efficiency=costs.at["Fischer-Tropsch", 'efficiency'],
        capital_cost=costs.at["Fischer-Tropsch", 'fixed'],
        efficiency2=-costs.at["oil", 'CO2 intensity'] * costs.at["Fischer-Tropsch", 'efficiency'],
        p_nom_extendable=True,
        lifetime=costs.at['Fischer-Tropsch', 'lifetime']
    )

    n.madd("Load",
        ["naphtha for industry"],
        bus=spatial.oil.nodes,
        carrier="naphtha for industry",
        p_set=industrial_demand.loc[nodes, "naphtha"].sum() / 8760
    )

    all_aviation = ["total international aviation", "total domestic aviation"]
    p_set = pop_weighted_energy_totals.loc[nodes, all_aviation].sum(axis=1).sum() * 1e6 / 8760

    n.madd("Load",
        ["kerosene for aviation"],
        bus=spatial.oil.nodes,
        carrier="kerosene for aviation",
        p_set=p_set
    )

    #NB: CO2 gets released again to atmosphere when plastics decay or kerosene is burned
    #except for the process emissions when naphtha is used for petrochemicals, which can be captured with other industry process emissions
    #tco2 per hour
    co2_release = ["naphtha for industry", "kerosene for aviation"]
    co2 = n.loads.loc[co2_release, "p_set"].sum() * costs.at["oil", 'CO2 intensity'] - industrial_demand.loc[nodes, "process emission from feedstock"].sum() / 8760

    n.add("Load",
        "oil emissions",
        bus="co2 atmosphere",
        carrier="oil emissions",
        p_set=-co2
    )

    # TODO simplify bus expression
    n.madd("Load",
        nodes,
        suffix=" low-temperature heat for industry",
        bus=[node + " urban central heat" if node + " urban central heat" in n.buses.index else node + " services urban decentral heat" for node in nodes],
        carrier="low-temperature heat for industry",
        p_set=industrial_demand.loc[nodes, "low-temperature heat"] / 8760
    )

    # remove today's industrial electricity demand by scaling down total electricity demand
    for ct in n.buses.country.dropna().unique():
        # TODO map onto n.bus.country
        loads_i = n.loads.index[(n.loads.index.str[:2] == ct) & (n.loads.carrier == "electricity")]
        if n.loads_t.p_set[loads_i].empty: continue
        factor = 1 - industrial_demand.loc[loads_i, "current electricity"].sum() / n.loads_t.p_set[loads_i].sum().sum()
        n.loads_t.p_set[loads_i] *= factor

    n.madd("Load",
        nodes,
        suffix=" industry electricity",
        bus=nodes,
        carrier="industry electricity",
        p_set=industrial_demand.loc[nodes, "electricity"] / 8760
    )

    n.add("Bus",
        "process emissions",
        location="EU",
        carrier="process emissions"
    )

    # this should be process emissions fossil+feedstock
    # then need load on atmosphere for feedstock emissions that are currently going to atmosphere via Link Fischer-Tropsch demand
    n.add("Load",
        "process emissions",
        bus="process emissions",
        carrier="process emissions",
        p_set=-industrial_demand.loc[nodes,["process emission", "process emission from feedstock"]].sum(axis=1).sum() / 8760
    )

    n.add("Link",
        "process emissions",
        bus0="process emissions",
        bus1="co2 atmosphere",
        carrier="process emissions",
        p_nom_extendable=True,
        efficiency=1.
    )

    #assume enough local waste heat for CC
    n.madd("Link",
        spatial.co2.locations,
        suffix=" process emissions CC",
        bus0="process emissions",
        bus1="co2 atmosphere",
        bus2=spatial.co2.nodes,
        carrier="process emissions CC",
        p_nom_extendable=True,
        capital_cost=costs.at["cement capture", "fixed"],
        efficiency=1 - costs.at["cement capture", "capture_rate"],
        efficiency2=costs.at["cement capture", "capture_rate"],
        lifetime=costs.at['cement capture', 'lifetime']
    )


def add_waste_heat(n):
    # TODO options?

    logger.info("Add possibility to use industrial waste heat in district heating")

    #AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if not urban_central.empty:
        urban_central = urban_central.str[:-len(" urban central heat")]

        # TODO what is the 0.95 and should it be a config option?
        if options['use_fischer_tropsch_waste_heat']:
            n.links.loc[urban_central + " Fischer-Tropsch", "bus3"] = urban_central + " urban central heat"
            n.links.loc[urban_central + " Fischer-Tropsch", "efficiency3"] = 0.95 - n.links.loc[urban_central + " Fischer-Tropsch", "efficiency"]

        if options['use_fuel_cell_waste_heat']:
            n.links.loc[urban_central + " H2 Fuel Cell", "bus2"] = urban_central + " urban central heat"
            n.links.loc[urban_central + " H2 Fuel Cell", "efficiency2"] = 0.95 - n.links.loc[urban_central + " H2 Fuel Cell", "efficiency"]


def add_agriculture(n, costs):

    logger.info('Add agriculture, forestry and fishing sector.')

    nodes = pop_layout.index

    # electricity

    n.madd("Load",
        nodes,
        suffix=" agriculture electricity",
        bus=nodes,
        carrier='agriculture electricity',
        p_set=pop_weighted_energy_totals.loc[nodes, "total agriculture electricity"] * 1e6 / 8760
    )

    # heat

    n.madd("Load",
        nodes,
        suffix=" agriculture heat",
        bus=nodes + " services rural heat",
        carrier="agriculture heat",
        p_set=pop_weighted_energy_totals.loc[nodes, "total agriculture heat"] * 1e6 / 8760
    )

    # machinery

    electric_share = get(options["agriculture_machinery_electric_share"], investment_year)
    assert electric_share <= 1.
    ice_share = 1 - electric_share

    machinery_nodal_energy = pop_weighted_energy_totals.loc[nodes, "total agriculture machinery"]

    if electric_share > 0:

        efficiency_gain = options["agriculture_machinery_fuel_efficiency"] / options["agriculture_machinery_electric_efficiency"]

        n.madd("Load",
            nodes,
            suffix=" agriculture machinery electric",
            bus=nodes,
            carrier="agriculture machinery electric",
            p_set=electric_share / efficiency_gain * machinery_nodal_energy * 1e6 / 8760,
        )

    if ice_share > 0:

        n.add("Load",
            "agriculture machinery oil",
            bus=spatial.oil.nodes,
            carrier="agriculture machinery oil",
            p_set=ice_share * machinery_nodal_energy.sum() * 1e6 / 8760
        )

        co2 = ice_share * machinery_nodal_energy.sum() * 1e6 / 8760 * costs.at["oil", 'CO2 intensity']

        n.add("Load",
            "agriculture machinery oil emissions",
            bus="co2 atmosphere",
            carrier="agriculture machinery oil emissions",
            p_set=-co2
        )


def decentral(n):
    """Removes the electricity transmission system."""
    n.lines.drop(n.lines.index, inplace=True)
    n.links.drop(n.links.index[n.links.carrier.isin(["DC", "B2B"])], inplace=True)


def remove_h2_network(n):

    n.links.drop(n.links.index[n.links.carrier.str.contains("H2 pipeline")], inplace=True)

    if "EU H2 Store" in n.stores.index:
        n.stores.drop("EU H2 Store", inplace=True)


def maybe_adjust_costs_and_potentials(n, opts):

    for o in opts:
        if "+" not in o: continue
        oo = o.split("+")
        carrier_list = np.hstack((n.generators.carrier.unique(), n.links.carrier.unique(),
                                n.stores.carrier.unique(), n.storage_units.carrier.unique()))
        suptechs = map(lambda c: c.split("-", 2)[0], carrier_list)
        if oo[0].startswith(tuple(suptechs)):
            carrier = oo[0]
            attr_lookup = {"p": "p_nom_max", "e": "e_nom_max", "c": "capital_cost"}
            attr = attr_lookup[oo[1][0]]
            factor = float(oo[1][1:])
            #beware if factor is 0 and p_nom_max is np.inf, 0*np.inf is nan
            if carrier == "AC":  # lines do not have carrier
                n.lines[attr] *= factor
            else:
                if attr == 'p_nom_max':
                    comps = {"Generator", "Link", "StorageUnit"}
                elif attr == 'e_nom_max':
                    comps = {"Store"}
                else:
                    comps = {"Generator", "Link", "StorageUnit", "Store"}
                for c in n.iterate_components(comps):
                    if carrier=='solar':
                        sel = c.df.carrier.str.contains(carrier) & ~c.df.carrier.str.contains("solar rooftop")
                    else:
                        sel = c.df.carrier.str.contains(carrier)
                    c.df.loc[sel,attr] *= factor
            print("changing", attr , "for", carrier, "by factor", factor)


# TODO this should rather be a config no wildcard
def limit_individual_line_extension(n, maxext):
    logger.info(f"limiting new HVAC and HVDC extensions to {maxext} MW")
    n.lines['s_nom_max'] = n.lines['s_nom'] + maxext
    hvdc = n.links.index[n.links.carrier == 'DC']
    n.links.loc[hvdc, 'p_nom_max'] = n.links.loc[hvdc, 'p_nom'] + maxext

#%%
if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'prepare_sector_network',
            simpl='',
            opts="",
            clusters="37",
            lv=1.0,
            sector_opts='Co2L0-168H-T-H-B-I-solar3-dist1',
            planning_horizons="2020",
        )

    logging.basicConfig(level=snakemake.config['logging_level'])

    options = snakemake.config["sector"]

    opts = snakemake.wildcards.sector_opts.split('-')

    investment_year = int(snakemake.wildcards.planning_horizons[-4:])

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    Nyears = n.snapshot_weightings.generators.sum() / 8760

    costs = prepare_costs(snakemake.input.costs,
                          snakemake.config['costs']['USD2013_to_EUR2013'],
                          snakemake.config['costs']['discountrate'],
                          Nyears,
                          snakemake.config['costs']['lifetime'])

    pop_weighted_energy_totals = pd.read_csv(snakemake.input.pop_weighted_energy_totals, index_col=0)

    patch_electricity_network(n)

    spatial = define_spatial(pop_layout.index, options)

    if snakemake.config["foresight"] == 'myopic':

        add_lifetime_wind_solar(n, costs)

        conventional = snakemake.config['existing_capacities']['conventional_carriers']
        for carrier in conventional:
            add_carrier_buses(n, carrier)

    add_co2_tracking(n, options)

    add_generation(n, costs)

    add_storage_and_grids(n, costs)

    # TODO merge with opts cost adjustment below
    for o in opts:
        if o[:4] == "wave":
            wave_cost_factor = float(o[4:].replace("p", ".").replace("m", "-"))
            print("Including wave generators with cost factor of", wave_cost_factor)
            add_wave(n, wave_cost_factor)
        if o[:4] == "dist":
            options['electricity_distribution_grid'] = True
            options['electricity_distribution_grid_cost_factor'] = float(o[4:].replace("p", ".").replace("m", "-"))
        if o == "biomasstransport":
            options["biomass_transport"] = True

    if "nodistrict" in opts:
        options["district_heating"]["progress"] = 0.0

    if "T" in opts:
        add_land_transport(n, costs)

    if "H" in opts:
        add_heat(n, costs)

    if "B" in opts:
        add_biomass(n, costs)

    if "I" in opts:
        add_industry(n, costs)

    if "I" in opts and "H" in opts:
        add_waste_heat(n)

    if "A" in opts:  # requires H and I
        add_agriculture(n, costs)

    if options['dac']:
        add_dac(n, costs)

    if "decentral" in opts:
        decentral(n)

    if "noH2network" in opts:
        remove_h2_network(n)

    if options["co2_network"]:
        add_co2_network(n, costs)

    for o in opts:
        m = re.match(r'^\d+h$', o, re.IGNORECASE)
        if m is not None:
            n = average_every_nhours(n, m.group(0))
            break

    limit_type = "config"
    limit = get(snakemake.config["co2_budget"], investment_year)
    for o in opts:
        if not "cb" in o: continue
        limit_type = "carbon budget"
        fn = snakemake.config['results_dir'] + snakemake.config['run'] + '/csvs/carbon_budget_distribution.csv'
        if not os.path.exists(fn):
            build_carbon_budget(o, fn)
        co2_cap = pd.read_csv(fn, index_col=0).squeeze()
        limit = co2_cap[investment_year]
        break
    for o in opts:
        if not "Co2L" in o: continue
        limit_type = "wildcard"
        limit = o[o.find("Co2L")+4:]
        limit = float(limit.replace("p", ".").replace("m", "-"))
        break
    print("Add CO2 limit from", limit_type)
    add_co2limit(n, Nyears, limit)

    for o in opts:
        if not o[:10] == 'linemaxext': continue
        maxext = float(o[10:]) * 1e3
        limit_individual_line_extension(n, maxext)
        break

    if options['electricity_distribution_grid']:
        insert_electricity_distribution_grid(n, costs)

    maybe_adjust_costs_and_potentials(n, opts)

    if options['gas_distribution_grid']:
        insert_gas_distribution_costs(n, costs)

    if options['electricity_grid_connection']:
        add_electricity_grid_connection(n, costs)

    n.export_to_netcdf(snakemake.output[0])
