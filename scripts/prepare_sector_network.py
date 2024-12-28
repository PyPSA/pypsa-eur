# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Adds all sector-coupling components to the network, including demand and supply
technologies for the buildings, transport and industry sectors.
"""

import logging
import os
from itertools import product
from types import SimpleNamespace

import networkx as nx
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import (
    configure_logging,
    get,
    set_scenario_config,
    update_config_from_wildcards,
)
from add_electricity import calculate_annuity, sanitize_carriers, sanitize_locations
from build_energy_totals import (
    build_co2_totals,
    build_eea_co2,
    build_eurostat,
    build_eurostat_co2,
)
from build_transport_demand import transport_degree_factor
from definitions.heat_sector import HeatSector
from definitions.heat_system import HeatSystem
from definitions.heat_system_type import HeatSystemType
from networkx.algorithms import complement
from networkx.algorithms.connectivity.edge_augmentation import k_edge_augmentation
from prepare_network import maybe_adjust_costs_and_potentials
from pypsa.geo import haversine_pts
from scipy.stats import beta

spatial = SimpleNamespace()
logger = logging.getLogger(__name__)


def define_spatial(nodes, options):
    """
    Namespace for spatial.

    Parameters
    ----------
    nodes : list-like
    """
    global spatial

    spatial.nodes = nodes

    # biomass

    spatial.biomass = SimpleNamespace()
    spatial.msw = SimpleNamespace()

    if options.get("biomass_spatial", options["biomass_transport"]):
        spatial.biomass.nodes = nodes + " solid biomass"
        spatial.biomass.nodes_unsustainable = nodes + " unsustainable solid biomass"
        spatial.biomass.bioliquids = nodes + " unsustainable bioliquids"
        spatial.biomass.locations = nodes
        spatial.biomass.industry = nodes + " solid biomass for industry"
        spatial.biomass.industry_cc = nodes + " solid biomass for industry CC"
        spatial.msw.nodes = nodes + " municipal solid waste"
        spatial.msw.locations = nodes
    else:
        spatial.biomass.nodes = ["EU solid biomass"]
        spatial.biomass.nodes_unsustainable = ["EU unsustainable solid biomass"]
        spatial.biomass.bioliquids = ["EU unsustainable bioliquids"]
        spatial.biomass.locations = ["EU"]
        spatial.biomass.industry = ["solid biomass for industry"]
        spatial.biomass.industry_cc = ["solid biomass for industry CC"]
        spatial.msw.nodes = ["EU municipal solid waste"]
        spatial.msw.locations = ["EU"]

    spatial.biomass.df = pd.DataFrame(vars(spatial.biomass), index=nodes)
    spatial.msw.df = pd.DataFrame(vars(spatial.msw), index=nodes)

    # co2

    spatial.co2 = SimpleNamespace()

    if options["co2_spatial"]:
        spatial.co2.nodes = nodes + " co2 stored"
        spatial.co2.locations = nodes
        spatial.co2.vents = nodes + " co2 vent"
        spatial.co2.process_emissions = nodes + " process emissions"
    else:
        spatial.co2.nodes = ["co2 stored"]
        spatial.co2.locations = ["EU"]
        spatial.co2.vents = ["co2 vent"]
        spatial.co2.process_emissions = ["process emissions"]

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
        spatial.gas.biogas_to_gas_cc = nodes + " biogas to gas CC"
    else:
        spatial.gas.nodes = ["EU gas"]
        spatial.gas.locations = ["EU"]
        spatial.gas.biogas = ["EU biogas"]
        spatial.gas.industry = ["gas for industry"]
        spatial.gas.biogas_to_gas = ["EU biogas to gas"]
        if options.get("biomass_spatial", options["biomass_transport"]):
            spatial.gas.biogas_to_gas_cc = nodes + " biogas to gas CC"
        else:
            spatial.gas.biogas_to_gas_cc = ["EU biogas to gas CC"]
        if options.get("co2_spatial", options["co2_network"]):
            spatial.gas.industry_cc = nodes + " gas for industry CC"
        else:
            spatial.gas.industry_cc = ["gas for industry CC"]

    spatial.gas.df = pd.DataFrame(vars(spatial.gas), index=nodes)

    # ammonia

    if options["ammonia"]:
        spatial.ammonia = SimpleNamespace()
        if options["ammonia"] == "regional":
            spatial.ammonia.nodes = nodes + " NH3"
            spatial.ammonia.locations = nodes
        else:
            spatial.ammonia.nodes = ["EU NH3"]
            spatial.ammonia.locations = ["EU"]

        spatial.ammonia.df = pd.DataFrame(vars(spatial.ammonia), index=nodes)

    # hydrogen
    spatial.h2 = SimpleNamespace()
    spatial.h2.nodes = nodes + " H2"
    spatial.h2.locations = nodes

    # methanol

    # beware: unlike other carriers, uses locations rather than locations+carriername
    # this allows to avoid separation between nodes and locations

    spatial.methanol = SimpleNamespace()

    spatial.methanol.nodes = ["EU methanol"]
    spatial.methanol.locations = ["EU"]

    if options["methanol"]["regional_methanol_demand"]:
        spatial.methanol.demand_locations = nodes
        spatial.methanol.industry = nodes + " industry methanol"
        spatial.methanol.shipping = nodes + " shipping methanol"
    else:
        spatial.methanol.demand_locations = ["EU"]
        spatial.methanol.shipping = ["EU shipping methanol"]
        spatial.methanol.industry = ["EU industry methanol"]

    # oil
    spatial.oil = SimpleNamespace()

    spatial.oil.nodes = ["EU oil"]
    spatial.oil.locations = ["EU"]

    if options["regional_oil_demand"]:
        spatial.oil.demand_locations = nodes
        spatial.oil.naphtha = nodes + " naphtha for industry"
        spatial.oil.kerosene = nodes + " kerosene for aviation"
        spatial.oil.shipping = nodes + " shipping oil"
        spatial.oil.agriculture_machinery = nodes + " agriculture machinery oil"
        spatial.oil.land_transport = nodes + " land transport oil"
    else:
        spatial.oil.demand_locations = ["EU"]
        spatial.oil.naphtha = ["EU naphtha for industry"]
        spatial.oil.kerosene = ["EU kerosene for aviation"]
        spatial.oil.shipping = ["EU shipping oil"]
        spatial.oil.agriculture_machinery = ["EU agriculture machinery oil"]
        spatial.oil.land_transport = ["EU land transport oil"]

    # uranium
    spatial.uranium = SimpleNamespace()
    spatial.uranium.nodes = ["EU uranium"]
    spatial.uranium.locations = ["EU"]

    # coal
    spatial.coal = SimpleNamespace()
    spatial.coal.nodes = ["EU coal"]
    spatial.coal.locations = ["EU"]

    if options["regional_coal_demand"]:
        spatial.coal.demand_locations = nodes
        spatial.coal.industry = nodes + " coal for industry"
    else:
        spatial.coal.demand_locations = ["EU"]
        spatial.coal.industry = ["EU coal for industry"]

    # lignite
    spatial.lignite = SimpleNamespace()
    spatial.lignite.nodes = ["EU lignite"]
    spatial.lignite.locations = ["EU"]

    # deep geothermal
    spatial.geothermal_heat = SimpleNamespace()
    spatial.geothermal_heat.nodes = ["EU enhanced geothermal systems"]
    spatial.geothermal_heat.locations = ["EU"]

    return spatial


spatial = SimpleNamespace()


def determine_emission_sectors(options):
    sectors = ["electricity"]
    if options["transport"]:
        sectors += ["rail non-elec", "road non-elec"]
    if options["heating"]:
        sectors += ["residential non-elec", "services non-elec"]
    if options["industry"]:
        sectors += [
            "industrial non-elec",
            "industrial processes",
            "domestic aviation",
            "international aviation",
            "domestic navigation",
            "international navigation",
        ]
    if options["agriculture"]:
        sectors += ["agriculture"]

    return sectors


def co2_emissions_year(
    countries, input_eurostat, options, emissions_scope, input_co2, year
):
    """
    Calculate CO2 emissions in one specific year (e.g. 1990 or 2018).
    """
    eea_co2 = build_eea_co2(input_co2, year, emissions_scope)

    eurostat = build_eurostat(input_eurostat, countries)

    # this only affects the estimation of CO2 emissions for BA, RS, AL, ME, MK, XK
    eurostat_co2 = build_eurostat_co2(eurostat, year)

    co2_totals = build_co2_totals(countries, eea_co2, eurostat_co2)

    sectors = determine_emission_sectors(options)

    co2_emissions = co2_totals.loc[countries, sectors].sum().sum()

    # convert MtCO2 to GtCO2
    co2_emissions *= 0.001

    return co2_emissions


# TODO: move to own rule with sector-opts wildcard?
def build_carbon_budget(o, input_eurostat, fn, emissions_scope, input_co2, options):
    """
    Distribute carbon budget following beta or exponential transition path.
    """

    if "be" in o:
        # beta decay
        carbon_budget = float(o[o.find("cb") + 2 : o.find("be")])
        be = float(o[o.find("be") + 2 :])
    if "ex" in o:
        # exponential decay
        carbon_budget = float(o[o.find("cb") + 2 : o.find("ex")])
        r = float(o[o.find("ex") + 2 :])

    countries = snakemake.params.countries

    e_1990 = co2_emissions_year(
        countries,
        input_eurostat,
        options,
        emissions_scope,
        input_co2,
        year=1990,
    )

    # emissions at the beginning of the path (last year available 2018)
    e_0 = co2_emissions_year(
        countries,
        input_eurostat,
        options,
        emissions_scope,
        input_co2,
        year=2018,
    )

    planning_horizons = snakemake.params.planning_horizons
    if not isinstance(planning_horizons, list):
        planning_horizons = [planning_horizons]
    t_0 = planning_horizons[0]

    if "be" in o:
        # final year in the path
        t_f = t_0 + (2 * carbon_budget / e_0).round(0)

        def beta_decay(t):
            cdf_term = (t - t_0) / (t_f - t_0)
            return (e_0 / e_1990) * (1 - beta.cdf(cdf_term, be, be))

        # emissions (relative to 1990)
        co2_cap = pd.Series({t: beta_decay(t) for t in planning_horizons}, name=o)

    if "ex" in o:
        T = carbon_budget / e_0
        m = (1 + np.sqrt(1 + r * T)) / T

        def exponential_decay(t):
            return (e_0 / e_1990) * (1 + (m + r) * (t - t_0)) * np.exp(-m * (t - t_0))

        co2_cap = pd.Series(
            {t: exponential_decay(t) for t in planning_horizons}, name=o
        )

    # TODO log in Snakefile
    csvs_folder = fn.rsplit("/", 1)[0]
    if not os.path.exists(csvs_folder):
        os.makedirs(csvs_folder)
    co2_cap.to_csv(fn, float_format="%.3f")


def add_lifetime_wind_solar(n, costs):
    """
    Add lifetime for solar and wind generators.
    """
    for carrier in ["solar", "onwind", "offwind"]:
        gen_i = n.generators.index.str.contains(carrier)
        n.generators.loc[gen_i, "lifetime"] = costs.at[carrier, "lifetime"]


def haversine(p):
    coord0 = n.buses.loc[p.bus0, ["x", "y"]].values
    coord1 = n.buses.loc[p.bus1, ["x", "y"]].values
    return 1.5 * haversine_pts(coord0, coord1)


def create_network_topology(
    n, prefix, carriers=["DC"], connector=" -> ", bidirectional=True
):
    """
    Create a network topology from transmission lines and link carrier
    selection.

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

    candidates = pd.concat(
        [n.lines[ln_attrs], n.links.loc[n.links.carrier.isin(carriers), lk_attrs]]
    ).fillna(0)

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


def update_wind_solar_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    line_length_factor: int | float = 1,
    landfall_lengths: dict = None,
) -> None:
    """
    Update costs for wind and solar generators added with pypsa-eur to those
    cost in the planning year.
    """

    if landfall_lengths is None:
        landfall_lengths = {}

    # NB: solar costs are also manipulated for rooftop
    # when distribution grid is inserted
    n.generators.loc[n.generators.carrier == "solar", "capital_cost"] = costs.at[
        "solar-utility", "fixed"
    ]

    n.generators.loc[n.generators.carrier == "onwind", "capital_cost"] = costs.at[
        "onwind", "fixed"
    ]

    # for offshore wind, need to calculated connection costs
    for connection in ["dc", "ac", "float"]:
        tech = "offwind-" + connection
        landfall_length = landfall_lengths.get(tech, 0.0)
        if tech not in n.generators.carrier.values:
            continue
        profile = snakemake.input["profile_offwind-" + connection]
        with xr.open_dataset(profile) as ds:
            # if-statement for compatibility with old profiles
            if "year" in ds.indexes:
                ds = ds.sel(year=ds.year.min(), drop=True)

            distance = ds["average_distance"].to_pandas()
            submarine_cost = costs.at[tech + "-connection-submarine", "fixed"]
            underground_cost = costs.at[tech + "-connection-underground", "fixed"]
            connection_cost = line_length_factor * (
                distance * submarine_cost + landfall_length * underground_cost
            )

            capital_cost = (
                costs.at["offwind", "fixed"]
                + costs.at[tech + "-station", "fixed"]
                + connection_cost
            )

            logger.info(
                "Added connection cost of {:0.0f}-{:0.0f} Eur/MW/a to {}".format(
                    connection_cost.min(), connection_cost.max(), tech
                )
            )

            n.generators.loc[n.generators.carrier == tech, "capital_cost"] = (
                capital_cost.rename(index=lambda node: node + " " + tech)
            )


def add_carrier_buses(n, carrier, nodes=None):
    """
    Add buses to connect e.g. coal, nuclear and oil plants.
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

    unit = "MWh_LHV" if carrier == "gas" else "MWh_th"
    # preliminary value for non-gas carriers to avoid zeros
    if carrier == "gas":
        capital_cost = costs.at["gas storage", "fixed"]
    elif carrier == "oil":
        # based on https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html
        mwh_per_m3 = 44.9 * 724 * 0.278 * 1e-3  # MJ/kg * kg/m3 * kWh/MJ * MWh/kWh
        capital_cost = (
            costs.at["General liquid hydrocarbon storage (product)", "fixed"]
            / mwh_per_m3
        )
    elif carrier == "methanol":
        # based on https://www.engineeringtoolbox.com/fossil-fuels-energy-content-d_1298.html
        mwh_per_m3 = 5.54 * 791 * 1e-3  # kWh/kg * kg/m3 * MWh/kWh
        capital_cost = (
            costs.at["General liquid hydrocarbon storage (product)", "fixed"]
            / mwh_per_m3
        )
    else:
        capital_cost = 0.1

    n.add("Bus", nodes, location=location, carrier=carrier, unit=unit)

    n.add(
        "Store",
        nodes + " Store",
        bus=nodes,
        e_nom_extendable=True,
        e_cyclic=True,
        carrier=carrier,
        capital_cost=capital_cost,
    )

    fossils = ["coal", "gas", "oil", "lignite"]
    if options["fossil_fuels"] and carrier in fossils:
        suffix = ""

        if carrier == "oil" and cf_industry["oil_refining_emissions"] > 0:
            n.add(
                "Bus",
                nodes + " primary",
                location=location,
                carrier=carrier + " primary",
                unit=unit,
            )

            n.add(
                "Link",
                nodes + " refining",
                bus0=nodes + " primary",
                bus1=nodes,
                bus2="co2 atmosphere",
                location=location,
                carrier=carrier + " refining",
                p_nom=1e6,
                efficiency=1
                - (
                    cf_industry["oil_refining_emissions"]
                    / costs.at[carrier, "CO2 intensity"]
                ),
                efficiency2=cf_industry["oil_refining_emissions"],
            )

            suffix = " primary"

        n.add(
            "Generator",
            nodes + suffix,
            bus=nodes + suffix,
            p_nom_extendable=True,
            carrier=carrier + suffix,
            marginal_cost=costs.at[carrier, "fuel"],
        )


# TODO: PyPSA-Eur merge issue
def remove_elec_base_techs(n):
    """
    Remove conventional generators (e.g. OCGT) and storage units (e.g.
    batteries and H2) from base electricity-only network, since they're added
    here differently using links.
    """
    for c in n.iterate_components(snakemake.params.pypsa_eur):
        to_keep = snakemake.params.pypsa_eur[c.name]
        to_remove = pd.Index(c.df.carrier.unique()).symmetric_difference(to_keep)
        if to_remove.empty:
            continue
        logger.info(f"Removing {c.list_name} with carrier {list(to_remove)}")
        names = c.df.index[c.df.carrier.isin(to_remove)]
        n.remove(c.name, names)
        n.carriers.drop(to_remove, inplace=True, errors="ignore")


# TODO: PyPSA-Eur merge issue
def remove_non_electric_buses(n):
    """
    Remove buses from pypsa-eur with carriers which are not AC buses.
    """
    if to_drop := list(n.buses.query("carrier not in ['AC', 'DC']").carrier.unique()):
        logger.info(f"Drop buses from PyPSA-Eur with carrier: {to_drop}")
        n.buses = n.buses[n.buses.carrier.isin(["AC", "DC"])]


def patch_electricity_network(n, costs, landfall_lengths):
    remove_elec_base_techs(n)
    remove_non_electric_buses(n)
    update_wind_solar_costs(n, costs, landfall_lengths=landfall_lengths)
    n.loads["carrier"] = "electricity"
    n.buses["location"] = n.buses.index
    n.buses["unit"] = "MWh_el"
    # remove trailing white space of load index until new PyPSA version after v0.18.
    n.loads.rename(lambda x: x.strip(), inplace=True)
    n.loads_t.p_set.rename(lambda x: x.strip(), axis=1, inplace=True)


def add_eu_bus(n, x=-5.5, y=46):
    """
    Add EU bus to the network.

    This cosmetic bus serves as a reference point for the location of
    the EU buses in the plots and summaries.
    """
    n.add("Bus", "EU", location="EU", x=x, y=y, carrier="none")
    n.add("Carrier", "none")


def add_co2_tracking(n, costs, options):
    # minus sign because opposite to how fossil fuels used:
    # CH4 burning puts CH4 down, atmosphere up
    n.add("Carrier", "co2", co2_emissions=-1.0)

    # this tracks CO2 in the atmosphere
    n.add("Bus", "co2 atmosphere", location="EU", carrier="co2", unit="t_co2")

    # can also be negative
    n.add(
        "Store",
        "co2 atmosphere",
        e_nom_extendable=True,
        e_min_pu=-1,
        carrier="co2",
        bus="co2 atmosphere",
    )

    # add CO2 tanks
    n.add(
        "Bus",
        spatial.co2.nodes,
        location=spatial.co2.locations,
        carrier="co2 stored",
        unit="t_co2",
    )

    n.add(
        "Store",
        spatial.co2.nodes,
        e_nom_extendable=True,
        capital_cost=costs.at["CO2 storage tank", "fixed"],
        carrier="co2 stored",
        e_cyclic=True,
        bus=spatial.co2.nodes,
    )
    n.add("Carrier", "co2 stored")

    # this tracks CO2 sequestered, e.g. underground
    sequestration_buses = pd.Index(spatial.co2.nodes).str.replace(
        " stored", " sequestered"
    )
    n.add(
        "Bus",
        sequestration_buses,
        location=spatial.co2.locations,
        carrier="co2 sequestered",
        unit="t_co2",
    )

    n.add(
        "Link",
        sequestration_buses,
        bus0=spatial.co2.nodes,
        bus1=sequestration_buses,
        carrier="co2 sequestered",
        efficiency=1.0,
        p_nom_extendable=True,
    )

    if options["regional_co2_sequestration_potential"]["enable"]:
        upper_limit = (
            options["regional_co2_sequestration_potential"]["max_size"] * 1e3
        )  # Mt
        annualiser = options["regional_co2_sequestration_potential"]["years_of_storage"]
        e_nom_max = pd.read_csv(
            snakemake.input.sequestration_potential, index_col=0
        ).squeeze()
        e_nom_max = (
            e_nom_max.reindex(spatial.co2.locations)
            .fillna(0.0)
            .clip(upper=upper_limit)
            .mul(1e6)
            / annualiser
        )  # t
        e_nom_max = e_nom_max.rename(index=lambda x: x + " co2 sequestered")
    else:
        e_nom_max = np.inf

    n.add(
        "Store",
        sequestration_buses,
        e_nom_extendable=True,
        e_nom_max=e_nom_max,
        capital_cost=options["co2_sequestration_cost"],
        marginal_cost=-0.1,
        bus=sequestration_buses,
        lifetime=options["co2_sequestration_lifetime"],
        carrier="co2 sequestered",
    )

    n.add("Carrier", "co2 sequestered")

    if options["co2_vent"]:
        n.add(
            "Link",
            spatial.co2.vents,
            bus0=spatial.co2.nodes,
            bus1="co2 atmosphere",
            carrier="co2 vent",
            efficiency=1.0,
            p_nom_extendable=True,
        )


def add_co2_network(n, costs):
    logger.info("Adding CO2 network.")
    co2_links = create_network_topology(n, "CO2 pipeline ")

    cost_onshore = (
        (1 - co2_links.underwater_fraction)
        * costs.at["CO2 pipeline", "fixed"]
        * co2_links.length
    )
    cost_submarine = (
        co2_links.underwater_fraction
        * costs.at["CO2 submarine pipeline", "fixed"]
        * co2_links.length
    )
    capital_cost = cost_onshore + cost_submarine
    cost_factor = snakemake.config["sector"]["co2_network_cost_factor"]
    capital_cost *= cost_factor

    n.add(
        "Link",
        co2_links.index,
        bus0=co2_links.bus0.values + " co2 stored",
        bus1=co2_links.bus1.values + " co2 stored",
        p_min_pu=-1,
        p_nom_extendable=True,
        length=co2_links.length.values,
        capital_cost=capital_cost.values,
        carrier="CO2 pipeline",
        lifetime=costs.at["CO2 pipeline", "lifetime"],
    )


def add_allam_gas(n, costs):
    logger.info("Adding Allam cycle gas power plants.")

    nodes = pop_layout.index

    n.add(
        "Link",
        nodes,
        suffix=" allam gas",
        bus0=spatial.gas.df.loc[nodes, "nodes"].values,
        bus1=nodes,
        bus2=spatial.co2.df.loc[nodes, "nodes"].values,
        bus3="co2 atmosphere",
        carrier="allam gas",
        p_nom_extendable=True,
        capital_cost=costs.at["allam", "fixed"] * costs.at["allam", "efficiency"],
        marginal_cost=costs.at["allam", "VOM"] * costs.at["allam", "efficiency"],
        efficiency=costs.at["allam", "efficiency"],
        efficiency2=0.98 * costs.at["gas", "CO2 intensity"],
        efficiency3=0.02 * costs.at["gas", "CO2 intensity"],
        lifetime=costs.at["allam", "lifetime"],
    )


def add_biomass_to_methanol(n, costs):
    n.add(
        "Link",
        spatial.biomass.nodes,
        suffix=" biomass-to-methanol",
        bus0=spatial.biomass.nodes,
        bus1=spatial.methanol.nodes,
        bus2="co2 atmosphere",
        carrier="biomass-to-methanol",
        lifetime=costs.at["biomass-to-methanol", "lifetime"],
        efficiency=costs.at["biomass-to-methanol", "efficiency"],
        efficiency2=-costs.at["solid biomass", "CO2 intensity"]
        + costs.at["biomass-to-methanol", "CO2 stored"],
        p_nom_extendable=True,
        capital_cost=costs.at["biomass-to-methanol", "fixed"]
        / costs.at["biomass-to-methanol", "efficiency"],
        marginal_cost=costs.loc["biomass-to-methanol", "VOM"]
        / costs.at["biomass-to-methanol", "efficiency"],
    )


def add_biomass_to_methanol_cc(n, costs):
    n.add(
        "Link",
        spatial.biomass.nodes,
        suffix=" biomass-to-methanol CC",
        bus0=spatial.biomass.nodes,
        bus1=spatial.methanol.nodes,
        bus2="co2 atmosphere",
        bus3=spatial.co2.nodes,
        carrier="biomass-to-methanol CC",
        lifetime=costs.at["biomass-to-methanol", "lifetime"],
        efficiency=costs.at["biomass-to-methanol", "efficiency"],
        efficiency2=-costs.at["solid biomass", "CO2 intensity"]
        + costs.at["biomass-to-methanol", "CO2 stored"]
        * (1 - costs.at["biomass-to-methanol", "capture rate"]),
        efficiency3=costs.at["biomass-to-methanol", "CO2 stored"]
        * costs.at["biomass-to-methanol", "capture rate"],
        p_nom_extendable=True,
        capital_cost=costs.at["biomass-to-methanol", "fixed"]
        / costs.at["biomass-to-methanol", "efficiency"]
        + costs.at["biomass CHP capture", "fixed"]
        * costs.at["biomass-to-methanol", "CO2 stored"],
        marginal_cost=costs.loc["biomass-to-methanol", "VOM"]
        / costs.at["biomass-to-methanol", "efficiency"],
    )


def add_methanol_to_power(n, costs, types=None):
    if types is None:
        types = {}

    nodes = pop_layout.index

    if types["allam"]:
        logger.info("Adding Allam cycle methanol power plants.")

        n.add(
            "Link",
            nodes,
            suffix=" allam methanol",
            bus0=spatial.methanol.nodes,
            bus1=nodes,
            bus2=spatial.co2.df.loc[nodes, "nodes"].values,
            bus3="co2 atmosphere",
            carrier="allam methanol",
            p_nom_extendable=True,
            capital_cost=costs.at["allam", "fixed"] * costs.at["allam", "efficiency"],
            marginal_cost=costs.at["allam", "VOM"] * costs.at["allam", "efficiency"],
            efficiency=costs.at["allam", "efficiency"],
            efficiency2=0.98 * costs.at["methanolisation", "carbondioxide-input"],
            efficiency3=0.02 * costs.at["methanolisation", "carbondioxide-input"],
            lifetime=25,
        )

    if types["ccgt"]:
        logger.info("Adding methanol CCGT power plants.")

        # efficiency * EUR/MW * (annuity + FOM)
        capital_cost = costs.at["CCGT", "efficiency"] * costs.at["CCGT", "fixed"]

        n.add(
            "Link",
            nodes,
            suffix=" CCGT methanol",
            bus0=spatial.methanol.nodes,
            bus1=nodes,
            bus2="co2 atmosphere",
            carrier="CCGT methanol",
            p_nom_extendable=True,
            capital_cost=capital_cost,
            marginal_cost=costs.at["CCGT", "VOM"],
            efficiency=costs.at["CCGT", "efficiency"],
            efficiency2=costs.at["methanolisation", "carbondioxide-input"],
            lifetime=costs.at["CCGT", "lifetime"],
        )

    if types["ccgt_cc"]:
        logger.info(
            "Adding methanol CCGT power plants with post-combustion carbon capture."
        )

        # TODO consider efficiency changes / energy inputs for CC

        # efficiency * EUR/MW * (annuity + FOM)
        capital_cost = costs.at["CCGT", "efficiency"] * costs.at["CCGT", "fixed"]

        capital_cost_cc = (
            capital_cost
            + costs.at["cement capture", "fixed"]
            * costs.at["methanolisation", "carbondioxide-input"]
        )

        n.add(
            "Link",
            nodes,
            suffix=" CCGT methanol CC",
            bus0=spatial.methanol.nodes,
            bus1=nodes,
            bus2=spatial.co2.df.loc[nodes, "nodes"].values,
            bus3="co2 atmosphere",
            carrier="CCGT methanol CC",
            p_nom_extendable=True,
            capital_cost=capital_cost_cc,
            marginal_cost=costs.at["CCGT", "VOM"],
            efficiency=costs.at["CCGT", "efficiency"],
            efficiency2=costs.at["cement capture", "capture_rate"]
            * costs.at["methanolisation", "carbondioxide-input"],
            efficiency3=(1 - costs.at["cement capture", "capture_rate"])
            * costs.at["methanolisation", "carbondioxide-input"],
            lifetime=costs.at["CCGT", "lifetime"],
        )

    if types["ocgt"]:
        logger.info("Adding methanol OCGT power plants.")

        n.add(
            "Link",
            nodes,
            suffix=" OCGT methanol",
            bus0=spatial.methanol.nodes,
            bus1=nodes,
            bus2="co2 atmosphere",
            carrier="OCGT methanol",
            p_nom_extendable=True,
            capital_cost=costs.at["OCGT", "fixed"] * costs.at["OCGT", "efficiency"],
            marginal_cost=costs.at["OCGT", "VOM"] * costs.at["OCGT", "efficiency"],
            efficiency=costs.at["OCGT", "efficiency"],
            efficiency2=costs.at["methanolisation", "carbondioxide-input"],
            lifetime=costs.at["OCGT", "lifetime"],
        )


def add_methanol_to_kerosene(n, costs):
    tech = "methanol-to-kerosene"

    logger.info(f"Adding {tech}.")

    capital_cost = costs.at[tech, "fixed"] / costs.at[tech, "methanol-input"]

    n.add(
        "Link",
        spatial.h2.locations,
        suffix=f" {tech}",
        carrier=tech,
        capital_cost=capital_cost,
        marginal_cost=costs.at[tech, "VOM"] / costs.at[tech, "methanol-input"],
        bus0=spatial.methanol.nodes,
        bus1=spatial.oil.kerosene,
        bus2=spatial.h2.nodes,
        bus3="co2 atmosphere",
        efficiency=1 / costs.at[tech, "methanol-input"],
        efficiency2=-costs.at[tech, "hydrogen-input"]
        / costs.at[tech, "methanol-input"],
        efficiency3=costs.at["oil", "CO2 intensity"] / costs.at[tech, "methanol-input"],
        p_nom_extendable=True,
        lifetime=costs.at[tech, "lifetime"],
    )


def add_methanol_reforming(n, costs):
    logger.info("Adding methanol steam reforming.")

    tech = "Methanol steam reforming"

    capital_cost = costs.at[tech, "fixed"] / costs.at[tech, "methanol-input"]

    n.add(
        "Link",
        spatial.h2.locations,
        suffix=f" {tech}",
        bus0=spatial.methanol.nodes,
        bus1=spatial.h2.nodes,
        bus2="co2 atmosphere",
        p_nom_extendable=True,
        capital_cost=capital_cost,
        efficiency=1 / costs.at[tech, "methanol-input"],
        efficiency2=costs.at["methanolisation", "carbondioxide-input"],
        carrier=tech,
        lifetime=costs.at[tech, "lifetime"],
    )


def add_methanol_reforming_cc(n, costs):
    logger.info("Adding methanol steam reforming with carbon capture.")

    tech = "Methanol steam reforming"

    # TODO: heat release and electricity demand for process and carbon capture
    # but the energy demands for carbon capture have not yet been added for other CC processes
    # 10.1016/j.rser.2020.110171: 0.129 kWh_e/kWh_H2, -0.09 kWh_heat/kWh_H2

    capital_cost = costs.at[tech, "fixed"] / costs.at[tech, "methanol-input"]

    capital_cost_cc = (
        capital_cost
        + costs.at["cement capture", "fixed"]
        * costs.at["methanolisation", "carbondioxide-input"]
    )

    n.add(
        "Link",
        spatial.h2.locations,
        suffix=f" {tech} CC",
        bus0=spatial.methanol.nodes,
        bus1=spatial.h2.nodes,
        bus2="co2 atmosphere",
        bus3=spatial.co2.nodes,
        p_nom_extendable=True,
        capital_cost=capital_cost_cc,
        efficiency=1 / costs.at[tech, "methanol-input"],
        efficiency2=(1 - costs.at["cement capture", "capture_rate"])
        * costs.at["methanolisation", "carbondioxide-input"],
        efficiency3=costs.at["cement capture", "capture_rate"]
        * costs.at["methanolisation", "carbondioxide-input"],
        carrier=f"{tech} CC",
        lifetime=costs.at[tech, "lifetime"],
    )


def add_dac(n, costs):
    heat_carriers = ["urban central heat", "services urban decentral heat"]
    heat_buses = n.buses.index[n.buses.carrier.isin(heat_carriers)]
    locations = n.buses.location[heat_buses]

    electricity_input = (
        costs.at["direct air capture", "electricity-input"]
        + costs.at["direct air capture", "compression-electricity-input"]
    )  # MWh_el / tCO2
    heat_input = (
        costs.at["direct air capture", "heat-input"]
        - costs.at["direct air capture", "compression-heat-output"]
    )  # MWh_th / tCO2

    n.add(
        "Link",
        heat_buses.str.replace(" heat", " DAC"),
        bus0=locations.values,
        bus1=heat_buses,
        bus2="co2 atmosphere",
        bus3=spatial.co2.df.loc[locations, "nodes"].values,
        carrier="DAC",
        capital_cost=costs.at["direct air capture", "fixed"] / electricity_input,
        efficiency=-heat_input / electricity_input,
        efficiency2=-1 / electricity_input,
        efficiency3=1 / electricity_input,
        p_nom_extendable=True,
        lifetime=costs.at["direct air capture", "lifetime"],
    )


def add_co2limit(n, options, nyears=1.0, limit=0.0):
    logger.info(f"Adding CO2 budget limit as per unit of 1990 levels of {limit}")

    countries = snakemake.params.countries

    sectors = determine_emission_sectors(options)

    # convert Mt to tCO2
    co2_totals = 1e6 * pd.read_csv(snakemake.input.co2_totals_name, index_col=0)

    co2_limit = co2_totals.loc[countries, sectors].sum().sum()

    co2_limit *= limit * nyears

    n.add(
        "GlobalConstraint",
        "CO2Limit",
        carrier_attribute="co2_emissions",
        sense="<=",
        type="co2_atmosphere",
        constant=co2_limit,
    )


def cycling_shift(df, steps=1):
    """
    Cyclic shift on index of pd.Series|pd.DataFrame by number of steps.
    """
    df = df.copy()
    new_index = np.roll(df.index, steps)
    df.values[:] = df.reindex(index=new_index).values
    return df


def prepare_costs(cost_file, params, nyears):
    # set all asset costs and other parameters
    costs = pd.read_csv(cost_file, index_col=[0, 1]).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3

    # min_count=1 is important to generate NaNs which are then filled by fillna
    costs = (
        costs.loc[:, "value"].unstack(level=1).groupby("technology").sum(min_count=1)
    )

    costs = costs.fillna(params["fill_values"])

    def annuity_factor(v):
        return calculate_annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100

    costs["fixed"] = [
        annuity_factor(v) * v["investment"] * nyears for i, v in costs.iterrows()
    ]

    return costs


def add_generation(n, costs):
    logger.info("Adding electricity generation")

    nodes = pop_layout.index

    conventionals = options["conventional_generation"]

    for generator, carrier in conventionals.items():
        carrier_nodes = vars(spatial)[carrier].nodes

        add_carrier_buses(n, carrier, carrier_nodes)

        n.add(
            "Link",
            nodes + " " + generator,
            bus0=carrier_nodes,
            bus1=nodes,
            bus2="co2 atmosphere",
            marginal_cost=costs.at[generator, "efficiency"]
            * costs.at[generator, "VOM"],  # NB: VOM is per MWel
            capital_cost=costs.at[generator, "efficiency"]
            * costs.at[generator, "fixed"],  # NB: fixed cost is per MWel
            p_nom_extendable=True,
            carrier=generator,
            efficiency=costs.at[generator, "efficiency"],
            efficiency2=costs.at[carrier, "CO2 intensity"],
            lifetime=costs.at[generator, "lifetime"],
        )


def add_ammonia(n, costs):
    logger.info("Adding ammonia carrier with synthesis, cracking and storage")

    nodes = pop_layout.index

    n.add("Carrier", "NH3")

    n.add(
        "Bus", spatial.ammonia.nodes, location=spatial.ammonia.locations, carrier="NH3"
    )

    n.add(
        "Link",
        nodes,
        suffix=" Haber-Bosch",
        bus0=nodes,
        bus1=spatial.ammonia.nodes,
        bus2=nodes + " H2",
        p_nom_extendable=True,
        carrier="Haber-Bosch",
        efficiency=1 / costs.at["Haber-Bosch", "electricity-input"],
        efficiency2=-costs.at["Haber-Bosch", "hydrogen-input"]
        / costs.at["Haber-Bosch", "electricity-input"],
        capital_cost=costs.at["Haber-Bosch", "fixed"]
        / costs.at["Haber-Bosch", "electricity-input"],
        marginal_cost=costs.at["Haber-Bosch", "VOM"]
        / costs.at["Haber-Bosch", "electricity-input"],
        lifetime=costs.at["Haber-Bosch", "lifetime"],
    )

    n.add(
        "Link",
        nodes,
        suffix=" ammonia cracker",
        bus0=spatial.ammonia.nodes,
        bus1=nodes + " H2",
        p_nom_extendable=True,
        carrier="ammonia cracker",
        efficiency=1 / cf_industry["MWh_NH3_per_MWh_H2_cracker"],
        capital_cost=costs.at["Ammonia cracker", "fixed"]
        / cf_industry["MWh_NH3_per_MWh_H2_cracker"],  # given per MW_H2
        lifetime=costs.at["Ammonia cracker", "lifetime"],
    )

    # Ammonia Storage
    n.add(
        "Store",
        spatial.ammonia.nodes,
        suffix=" ammonia store",
        bus=spatial.ammonia.nodes,
        e_nom_extendable=True,
        e_cyclic=True,
        carrier="ammonia store",
        capital_cost=costs.at["NH3 (l) storage tank incl. liquefaction", "fixed"],
        lifetime=costs.at["NH3 (l) storage tank incl. liquefaction", "lifetime"],
    )


def insert_electricity_distribution_grid(n, costs):
    # TODO pop_layout?
    # TODO options?

    nodes = pop_layout.index

    n.add(
        "Bus",
        nodes + " low voltage",
        location=nodes,
        carrier="low voltage",
        unit="MWh_el",
    )

    n.add(
        "Link",
        nodes + " electricity distribution grid",
        bus0=nodes,
        bus1=nodes + " low voltage",
        p_nom_extendable=True,
        p_min_pu=-1,
        carrier="electricity distribution grid",
        efficiency=1,
        lifetime=costs.at["electricity distribution grid", "lifetime"],
        capital_cost=costs.at["electricity distribution grid", "fixed"],
    )

    # deduct distribution losses from electricity demand as these are included in total load
    # https://nbviewer.org/github/Open-Power-System-Data/datapackage_timeseries/blob/2020-10-06/main.ipynb
    if (
        efficiency := options["transmission_efficiency"]
        .get("electricity distribution grid", {})
        .get("efficiency_static")
    ):
        logger.info(
            f"Deducting distribution losses from electricity demand: {np.around(100*(1-efficiency), decimals=2)}%"
        )
        n.loads_t.p_set.loc[:, n.loads.carrier == "electricity"] *= efficiency

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
    n.generators.loc[solar, "capital_cost"] = costs.at["solar-utility", "fixed"]
    pop_solar = pop_layout.total.rename(index=lambda x: x + " solar")

    # add max solar rooftop potential assuming 0.1 kW/m2 and 20 m2/person,
    # i.e. 2 kW/person (population data is in thousands of people) so we get MW
    potential = 0.1 * 20 * pop_solar

    n.add(
        "Generator",
        solar,
        suffix=" rooftop",
        bus=n.generators.loc[solar, "bus"] + " low voltage",
        carrier="solar rooftop",
        p_nom_extendable=True,
        p_nom_max=potential.loc[solar],
        marginal_cost=n.generators.loc[solar, "marginal_cost"],
        capital_cost=costs.at["solar-rooftop", "fixed"],
        efficiency=n.generators.loc[solar, "efficiency"],
        p_max_pu=n.generators_t.p_max_pu[solar],
        lifetime=costs.at["solar-rooftop", "lifetime"],
    )

    n.add("Carrier", "home battery")

    n.add(
        "Bus",
        nodes + " home battery",
        location=nodes,
        carrier="home battery",
        unit="MWh_el",
    )

    n.add(
        "Store",
        nodes + " home battery",
        bus=nodes + " home battery",
        location=nodes,
        e_cyclic=True,
        e_nom_extendable=True,
        carrier="home battery",
        capital_cost=costs.at["home battery storage", "fixed"],
        lifetime=costs.at["battery storage", "lifetime"],
    )

    n.add(
        "Link",
        nodes + " home battery charger",
        bus0=nodes + " low voltage",
        bus1=nodes + " home battery",
        carrier="home battery charger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        capital_cost=costs.at["home battery inverter", "fixed"],
        p_nom_extendable=True,
        lifetime=costs.at["battery inverter", "lifetime"],
    )

    n.add(
        "Link",
        nodes + " home battery discharger",
        bus0=nodes + " home battery",
        bus1=nodes + " low voltage",
        carrier="home battery discharger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        marginal_cost=options["marginal_cost_storage"],
        p_nom_extendable=True,
        lifetime=costs.at["battery inverter", "lifetime"],
    )


def insert_gas_distribution_costs(n, costs):
    # TODO options?

    f_costs = options["gas_distribution_grid_cost_factor"]

    logger.info(
        f"Inserting gas distribution grid with investment cost factor of {f_costs}"
    )

    capital_cost = costs.at["electricity distribution grid", "fixed"] * f_costs

    # gas boilers
    gas_b = n.links.index[
        n.links.carrier.str.contains("gas boiler")
        & (~n.links.carrier.str.contains("urban central"))
    ]
    n.links.loc[gas_b, "capital_cost"] += capital_cost

    # micro CHPs
    mchp = n.links.index[n.links.carrier.str.contains("micro gas")]
    n.links.loc[mchp, "capital_cost"] += capital_cost


def add_electricity_grid_connection(n, costs):
    carriers = ["onwind", "solar", "solar-hsat"]

    gens = n.generators.index[n.generators.carrier.isin(carriers)]

    n.generators.loc[gens, "capital_cost"] += costs.at[
        "electricity grid connection", "fixed"
    ]


def add_storage_and_grids(n, costs):
    logger.info("Add hydrogen storage")

    nodes = pop_layout.index

    n.add("Carrier", "H2")

    n.add("Bus", nodes + " H2", location=nodes, carrier="H2", unit="MWh_LHV")

    n.add(
        "Link",
        nodes + " H2 Electrolysis",
        bus1=nodes + " H2",
        bus0=nodes,
        p_nom_extendable=True,
        carrier="H2 Electrolysis",
        efficiency=costs.at["electrolysis", "efficiency"],
        capital_cost=costs.at["electrolysis", "fixed"],
        lifetime=costs.at["electrolysis", "lifetime"],
    )

    if options["hydrogen_fuel_cell"]:
        logger.info("Adding hydrogen fuel cell for re-electrification.")

        n.add(
            "Link",
            nodes + " H2 Fuel Cell",
            bus0=nodes + " H2",
            bus1=nodes,
            p_nom_extendable=True,
            carrier="H2 Fuel Cell",
            efficiency=costs.at["fuel cell", "efficiency"],
            capital_cost=costs.at["fuel cell", "fixed"]
            * costs.at["fuel cell", "efficiency"],  # NB: fixed cost is per MWel
            lifetime=costs.at["fuel cell", "lifetime"],
        )

    if options["hydrogen_turbine"]:
        logger.info(
            "Adding hydrogen turbine for re-electrification. Assuming OCGT technology costs."
        )
        # TODO: perhaps replace with hydrogen-specific technology assumptions.

        n.add(
            "Link",
            nodes + " H2 turbine",
            bus0=nodes + " H2",
            bus1=nodes,
            p_nom_extendable=True,
            carrier="H2 turbine",
            efficiency=costs.at["OCGT", "efficiency"],
            capital_cost=costs.at["OCGT", "fixed"]
            * costs.at["OCGT", "efficiency"],  # NB: fixed cost is per MWel
            marginal_cost=costs.at["OCGT", "VOM"],
            lifetime=costs.at["OCGT", "lifetime"],
        )

    cavern_types = snakemake.params.sector["hydrogen_underground_storage_locations"]
    h2_caverns = pd.read_csv(snakemake.input.h2_cavern, index_col=0)

    if (
        not h2_caverns.empty
        and options["hydrogen_underground_storage"]
        and set(cavern_types).intersection(h2_caverns.columns)
    ):
        h2_caverns = h2_caverns[cavern_types].sum(axis=1)

        # only use sites with at least 2 TWh potential
        h2_caverns = h2_caverns[h2_caverns > 2]

        # convert TWh to MWh
        h2_caverns = h2_caverns * 1e6

        # clip at 1000 TWh for one location
        h2_caverns.clip(upper=1e9, inplace=True)

        logger.info("Add hydrogen underground storage")

        h2_capital_cost = costs.at["hydrogen storage underground", "fixed"]

        n.add(
            "Store",
            h2_caverns.index + " H2 Store",
            bus=h2_caverns.index + " H2",
            e_nom_extendable=True,
            e_nom_max=h2_caverns.values,
            e_cyclic=True,
            carrier="H2 Store",
            capital_cost=h2_capital_cost,
            lifetime=costs.at["hydrogen storage underground", "lifetime"],
        )

    # hydrogen stored overground (where not already underground)
    h2_capital_cost = costs.at[
        "hydrogen storage tank type 1 including compressor", "fixed"
    ]
    nodes_overground = h2_caverns.index.symmetric_difference(nodes)

    n.add(
        "Store",
        nodes_overground + " H2 Store",
        bus=nodes_overground + " H2",
        e_nom_extendable=True,
        e_cyclic=True,
        carrier="H2 Store",
        capital_cost=h2_capital_cost,
    )

    if options["gas_network"] or options["H2_retrofit"]:
        fn = snakemake.input.clustered_gas_network
        gas_pipes = pd.read_csv(fn, index_col=0)

    if options["gas_network"]:
        logger.info(
            "Add natural gas infrastructure, incl. LNG terminals, production, storage and entry-points."
        )

        if options["H2_retrofit"]:
            gas_pipes["p_nom_max"] = gas_pipes.p_nom
            gas_pipes["p_nom_min"] = 0.0
            # 0.1 EUR/MWkm/a to prefer decommissioning to address degeneracy
            gas_pipes["capital_cost"] = 0.1 * gas_pipes.length
            gas_pipes["p_nom_extendable"] = True
        else:
            gas_pipes["p_nom_max"] = np.inf
            gas_pipes["p_nom_min"] = gas_pipes.p_nom
            gas_pipes["capital_cost"] = (
                gas_pipes.length * costs.at["CH4 (g) pipeline", "fixed"]
            )
            gas_pipes["p_nom_extendable"] = False

        n.add(
            "Link",
            gas_pipes.index,
            bus0=gas_pipes.bus0 + " gas",
            bus1=gas_pipes.bus1 + " gas",
            p_min_pu=gas_pipes.p_min_pu,
            p_nom=gas_pipes.p_nom,
            p_nom_extendable=gas_pipes.p_nom_extendable,
            p_nom_max=gas_pipes.p_nom_max,
            p_nom_min=gas_pipes.p_nom_min,
            length=gas_pipes.length,
            capital_cost=gas_pipes.capital_cost,
            tags=gas_pipes.name,
            carrier="gas pipeline",
            lifetime=np.inf,
        )

        # remove fossil generators where there is neither
        # production, LNG terminal, nor entry-point beyond system scope

        fn = snakemake.input.gas_input_nodes_simplified
        gas_input_nodes = pd.read_csv(fn, index_col=0)

        unique = gas_input_nodes.index.unique()
        gas_i = n.generators.carrier == "gas"
        internal_i = ~n.generators.bus.map(n.buses.location).isin(unique)

        remove_i = n.generators[gas_i & internal_i].index
        n.generators.drop(remove_i, inplace=True)

        input_types = ["lng", "pipeline", "production"]
        p_nom = gas_input_nodes[input_types].sum(axis=1).rename(lambda x: x + " gas")
        n.generators.loc[gas_i, "p_nom_extendable"] = False
        n.generators.loc[gas_i, "p_nom"] = p_nom

        # add existing gas storage capacity
        gas_i = n.stores.carrier == "gas"
        e_nom = (
            gas_input_nodes["storage"]
            .rename(lambda x: x + " gas Store")
            .reindex(n.stores.index)
            .fillna(0.0)
            * 1e3
        )  # MWh_LHV
        e_nom.clip(
            upper=e_nom.quantile(0.98), inplace=True
        )  # limit extremely large storage
        n.stores.loc[gas_i, "e_nom_min"] = e_nom

        # add candidates for new gas pipelines to achieve full connectivity

        G = nx.Graph()

        gas_buses = n.buses.loc[n.buses.carrier == "gas", "location"]
        G.add_nodes_from(np.unique(gas_buses.values))

        sel = gas_pipes.p_nom > 1500
        attrs = ["bus0", "bus1", "length"]
        G.add_weighted_edges_from(gas_pipes.loc[sel, attrs].values)

        # find all complement edges
        complement_edges = pd.DataFrame(complement(G).edges, columns=["bus0", "bus1"])
        complement_edges["length"] = complement_edges.apply(haversine, axis=1)

        # apply k_edge_augmentation weighted by length of complement edges
        k_edge = options["gas_network_connectivity_upgrade"]
        if augmentation := list(
            k_edge_augmentation(G, k_edge, avail=complement_edges.values)
        ):
            new_gas_pipes = pd.DataFrame(augmentation, columns=["bus0", "bus1"])
            new_gas_pipes["length"] = new_gas_pipes.apply(haversine, axis=1)

            new_gas_pipes.index = new_gas_pipes.apply(
                lambda x: f"gas pipeline new {x.bus0} <-> {x.bus1}", axis=1
            )

            n.add(
                "Link",
                new_gas_pipes.index,
                bus0=new_gas_pipes.bus0 + " gas",
                bus1=new_gas_pipes.bus1 + " gas",
                p_min_pu=-1,  # new gas pipes are bidirectional
                p_nom_extendable=True,
                length=new_gas_pipes.length,
                capital_cost=new_gas_pipes.length
                * costs.at["CH4 (g) pipeline", "fixed"],
                carrier="gas pipeline new",
                lifetime=costs.at["CH4 (g) pipeline", "lifetime"],
            )

    if options["H2_retrofit"]:
        logger.info("Add retrofitting options of existing CH4 pipes to H2 pipes.")

        fr = "gas pipeline"
        to = "H2 pipeline retrofitted"
        h2_pipes = gas_pipes.rename(index=lambda x: x.replace(fr, to))

        n.add(
            "Link",
            h2_pipes.index,
            bus0=h2_pipes.bus0 + " H2",
            bus1=h2_pipes.bus1 + " H2",
            p_min_pu=-1.0,  # allow that all H2 retrofit pipelines can be used in both directions
            p_nom_max=h2_pipes.p_nom * options["H2_retrofit_capacity_per_CH4"],
            p_nom_extendable=True,
            length=h2_pipes.length,
            capital_cost=costs.at["H2 (g) pipeline repurposed", "fixed"]
            * h2_pipes.length,
            tags=h2_pipes.name,
            carrier="H2 pipeline retrofitted",
            lifetime=costs.at["H2 (g) pipeline repurposed", "lifetime"],
        )

    if options["H2_network"]:
        logger.info("Add options for new hydrogen pipelines.")

        h2_pipes = create_network_topology(
            n, "H2 pipeline ", carriers=["DC", "gas pipeline"]
        )

        # TODO Add efficiency losses
        n.add(
            "Link",
            h2_pipes.index,
            bus0=h2_pipes.bus0.values + " H2",
            bus1=h2_pipes.bus1.values + " H2",
            p_min_pu=-1,
            p_nom_extendable=True,
            length=h2_pipes.length.values,
            capital_cost=costs.at["H2 (g) pipeline", "fixed"] * h2_pipes.length.values,
            carrier="H2 pipeline",
            lifetime=costs.at["H2 (g) pipeline", "lifetime"],
        )

    n.add("Carrier", "battery")

    n.add("Bus", nodes + " battery", location=nodes, carrier="battery", unit="MWh_el")

    n.add(
        "Store",
        nodes + " battery",
        bus=nodes + " battery",
        e_cyclic=True,
        e_nom_extendable=True,
        carrier="battery",
        capital_cost=costs.at["battery storage", "fixed"],
        lifetime=costs.at["battery storage", "lifetime"],
    )

    n.add(
        "Link",
        nodes + " battery charger",
        bus0=nodes,
        bus1=nodes + " battery",
        carrier="battery charger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        capital_cost=costs.at["battery inverter", "fixed"],
        p_nom_extendable=True,
        lifetime=costs.at["battery inverter", "lifetime"],
    )

    n.add(
        "Link",
        nodes + " battery discharger",
        bus0=nodes + " battery",
        bus1=nodes,
        carrier="battery discharger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        marginal_cost=options["marginal_cost_storage"],
        p_nom_extendable=True,
        lifetime=costs.at["battery inverter", "lifetime"],
    )

    if options["methanation"]:
        n.add(
            "Link",
            spatial.nodes,
            suffix=" Sabatier",
            bus0=nodes + " H2",
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            p_nom_extendable=True,
            carrier="Sabatier",
            p_min_pu=options["min_part_load_methanation"],
            efficiency=costs.at["methanation", "efficiency"],
            efficiency2=-costs.at["methanation", "efficiency"]
            * costs.at["gas", "CO2 intensity"],
            capital_cost=costs.at["methanation", "fixed"]
            * costs.at["methanation", "efficiency"],  # costs given per kW_gas
            lifetime=costs.at["methanation", "lifetime"],
        )

    if options["coal_cc"]:
        n.add(
            "Link",
            spatial.nodes,
            suffix=" coal CC",
            bus0=spatial.coal.nodes,
            bus1=spatial.nodes,
            bus2="co2 atmosphere",
            bus3=spatial.co2.nodes,
            marginal_cost=costs.at["coal", "efficiency"]
            * costs.at["coal", "VOM"],  # NB: VOM is per MWel
            capital_cost=costs.at["coal", "efficiency"] * costs.at["coal", "fixed"]
            + costs.at["biomass CHP capture", "fixed"]
            * costs.at["coal", "CO2 intensity"],  # NB: fixed cost is per MWel
            p_nom_extendable=True,
            carrier="coal",
            efficiency=costs.at["coal", "efficiency"],
            efficiency2=costs.at["coal", "CO2 intensity"]
            * (1 - costs.at["biomass CHP capture", "capture_rate"]),
            efficiency3=costs.at["coal", "CO2 intensity"]
            * costs.at["biomass CHP capture", "capture_rate"],
            lifetime=costs.at["coal", "lifetime"],
        )

    if options["SMR_cc"]:
        n.add(
            "Link",
            spatial.nodes,
            suffix=" SMR CC",
            bus0=spatial.gas.nodes,
            bus1=nodes + " H2",
            bus2="co2 atmosphere",
            bus3=spatial.co2.nodes,
            p_nom_extendable=True,
            carrier="SMR CC",
            efficiency=costs.at["SMR CC", "efficiency"],
            efficiency2=costs.at["gas", "CO2 intensity"] * (1 - options["cc_fraction"]),
            efficiency3=costs.at["gas", "CO2 intensity"] * options["cc_fraction"],
            capital_cost=costs.at["SMR CC", "fixed"],
            lifetime=costs.at["SMR CC", "lifetime"],
        )

    if options["SMR"]:
        n.add(
            "Link",
            nodes + " SMR",
            bus0=spatial.gas.nodes,
            bus1=nodes + " H2",
            bus2="co2 atmosphere",
            p_nom_extendable=True,
            carrier="SMR",
            efficiency=costs.at["SMR", "efficiency"],
            efficiency2=costs.at["gas", "CO2 intensity"],
            capital_cost=costs.at["SMR", "fixed"],
            lifetime=costs.at["SMR", "lifetime"],
        )


def check_land_transport_shares(shares):
    # Sums up the shares, ignoring None values
    total_share = sum(filter(None, shares))
    if total_share != 1:
        logger.warning(
            f"Total land transport shares sum up to {total_share:.2%},"
            "corresponding to increased or decreased demand assumptions."
        )


def get_temp_efficency(
    car_efficiency,
    temperature,
    deadband_lw,
    deadband_up,
    degree_factor_lw,
    degree_factor_up,
):
    """
    Correct temperature depending on heating and cooling for respective car
    type.
    """
    # temperature correction for EVs
    dd = transport_degree_factor(
        temperature,
        deadband_lw,
        deadband_up,
        degree_factor_lw,
        degree_factor_up,
    )

    temp_eff = 1 / (1 + dd)

    return car_efficiency * temp_eff


def add_EVs(
    n,
    avail_profile,
    dsm_profile,
    p_set,
    electric_share,
    number_cars,
    temperature,
):
    n.add("Carrier", "EV battery")

    n.add(
        "Bus",
        spatial.nodes,
        suffix=" EV battery",
        location=spatial.nodes,
        carrier="EV battery",
        unit="MWh_el",
    )

    car_efficiency = options["transport_electric_efficiency"]

    # temperature corrected efficiency
    efficiency = get_temp_efficency(
        car_efficiency,
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["EV_lower_degree_factor"],
        options["EV_upper_degree_factor"],
    )

    p_shifted = (p_set + cycling_shift(p_set, 1) + cycling_shift(p_set, 2)) / 3

    cyclic_eff = p_set.div(p_shifted)

    efficiency *= cyclic_eff

    profile = electric_share * p_set.div(efficiency)

    n.add(
        "Load",
        spatial.nodes,
        suffix=" land transport EV",
        bus=spatial.nodes + " EV battery",
        carrier="land transport EV",
        p_set=profile.loc[n.snapshots],
    )

    p_nom = number_cars * options["bev_charge_rate"] * electric_share

    n.add(
        "Link",
        spatial.nodes,
        suffix=" BEV charger",
        bus0=spatial.nodes,
        bus1=spatial.nodes + " EV battery",
        p_nom=p_nom,
        carrier="BEV charger",
        p_max_pu=avail_profile.loc[n.snapshots, spatial.nodes],
        lifetime=1,
        efficiency=options["bev_charge_efficiency"],
    )

    if options["bev_dsm"]:
        e_nom = (
            number_cars
            * options["bev_energy"]
            * options["bev_dsm_availability"]
            * electric_share
        )

        n.add(
            "Store",
            spatial.nodes,
            suffix=" EV battery",
            bus=spatial.nodes + " EV battery",
            carrier="EV battery",
            e_cyclic=True,
            e_nom=e_nom,
            e_max_pu=1,
            e_min_pu=dsm_profile.loc[n.snapshots, spatial.nodes],
        )

        if options["v2g"]:
            n.add(
                "Link",
                spatial.nodes,
                suffix=" V2G",
                bus1=spatial.nodes,
                bus0=spatial.nodes + " EV battery",
                p_nom=p_nom * options["bev_dsm_availability"],
                carrier="V2G",
                p_max_pu=avail_profile.loc[n.snapshots, spatial.nodes],
                lifetime=1,
                efficiency=options["bev_charge_efficiency"],
            )


def add_fuel_cell_cars(n, p_set, fuel_cell_share, temperature):
    car_efficiency = options["transport_fuel_cell_efficiency"]

    # temperature corrected efficiency
    efficiency = get_temp_efficency(
        car_efficiency,
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["ICE_lower_degree_factor"],
        options["ICE_upper_degree_factor"],
    )

    profile = fuel_cell_share * p_set.div(efficiency)

    n.add(
        "Load",
        spatial.nodes,
        suffix=" land transport fuel cell",
        bus=spatial.h2.nodes,
        carrier="land transport fuel cell",
        p_set=profile.loc[n.snapshots],
    )


def add_ice_cars(n, p_set, ice_share, temperature):
    add_carrier_buses(n, "oil")

    car_efficiency = options["transport_ice_efficiency"]

    # temperature corrected efficiency
    efficiency = get_temp_efficency(
        car_efficiency,
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["ICE_lower_degree_factor"],
        options["ICE_upper_degree_factor"],
    )

    profile = ice_share * p_set.div(efficiency).rename(
        columns=lambda x: x + " land transport oil"
    )

    if not options["regional_oil_demand"]:
        profile = profile.sum(axis=1).to_frame(name="EU land transport oil")

    n.add(
        "Bus",
        spatial.oil.land_transport,
        location=spatial.oil.demand_locations,
        carrier="land transport oil",
        unit="land transport",
    )

    n.add(
        "Load",
        spatial.oil.land_transport,
        bus=spatial.oil.land_transport,
        carrier="land transport oil",
        p_set=profile.loc[n.snapshots],
    )

    n.add(
        "Link",
        spatial.oil.land_transport,
        bus0=spatial.oil.nodes,
        bus1=spatial.oil.land_transport,
        bus2="co2 atmosphere",
        carrier="land transport oil",
        efficiency2=costs.at["oil", "CO2 intensity"],
        p_nom_extendable=True,
    )


def add_land_transport(n, costs):
    logger.info("Add land transport")

    # read in transport demand in units driven km [100 km]
    transport = pd.read_csv(
        snakemake.input.transport_demand, index_col=0, parse_dates=True
    )
    number_cars = pd.read_csv(snakemake.input.transport_data, index_col=0)[
        "number cars"
    ]
    avail_profile = pd.read_csv(
        snakemake.input.avail_profile, index_col=0, parse_dates=True
    )
    dsm_profile = pd.read_csv(
        snakemake.input.dsm_profile, index_col=0, parse_dates=True
    )

    # exogenous share of passenger car type
    engine_types = ["fuel_cell", "electric", "ice"]
    shares = pd.Series()
    for engine in engine_types:
        shares[engine] = get(options[f"land_transport_{engine}_share"], investment_year)
        logger.info(f"{engine} share: {shares[engine]*100}%")

    check_land_transport_shares(shares)

    p_set = transport[spatial.nodes]

    # temperature for correction factor for heating/cooling
    temperature = xr.open_dataarray(snakemake.input.temp_air_total).to_pandas()

    if shares["electric"] > 0:
        add_EVs(
            n,
            avail_profile,
            dsm_profile,
            p_set,
            shares["electric"],
            number_cars,
            temperature,
        )

    if shares["fuel_cell"] > 0:
        add_fuel_cell_cars(n, p_set, shares["fuel_cell"], temperature)

    if shares["ice"] > 0:
        add_ice_cars(n, p_set, shares["ice"], temperature)


def build_heat_demand(n):
    heat_demand_shape = (
        xr.open_dataset(snakemake.input.hourly_heat_demand_total)
        .to_dataframe()
        .unstack(level=1)
    )

    sectors = [sector.value for sector in HeatSector]
    uses = ["water", "space"]

    heat_demand = {}
    electric_heat_supply = {}
    for sector, use in product(sectors, uses):
        name = f"{sector} {use}"

        # efficiency for final energy to thermal energy service
        eff = pop_weighted_energy_totals.index.str[:2].map(
            heating_efficiencies[f"total {sector} {use} efficiency"]
        )

        heat_demand[name] = (
            heat_demand_shape[name] / heat_demand_shape[name].sum()
        ).multiply(pop_weighted_energy_totals[f"total {sector} {use}"] * eff) * 1e6
        electric_heat_supply[name] = (
            heat_demand_shape[name] / heat_demand_shape[name].sum()
        ).multiply(pop_weighted_energy_totals[f"electricity {sector} {use}"]) * 1e6

    heat_demand = pd.concat(heat_demand, axis=1)
    electric_heat_supply = pd.concat(electric_heat_supply, axis=1)

    # subtract from electricity load since heat demand already in heat_demand
    electric_nodes = n.loads.index[n.loads.carrier == "electricity"]
    n.loads_t.p_set[electric_nodes] = (
        n.loads_t.p_set[electric_nodes]
        - electric_heat_supply.T.groupby(level=1).sum().T[electric_nodes]
    )

    return heat_demand


def add_heat(
    n: pypsa.Network,
    costs: pd.DataFrame,
    cop: xr.DataArray,
    direct_heat_source_utilisation_profile: xr.DataArray,
):
    """
    Add heat sector to the network.

    Parameters:
        n (pypsa.Network): The PyPSA network object.
        costs (pd.DataFrame): DataFrame containing cost information.
        cop (xr.DataArray): DataArray containing coefficient of performance (COP) values.

    Returns:
        None
    """
    logger.info("Add heat sector")

    sectors = [sector.value for sector in HeatSector]

    heat_demand = build_heat_demand(n)

    district_heat_info = pd.read_csv(snakemake.input.district_heat_share, index_col=0)
    dist_fraction = district_heat_info["district fraction of node"]
    urban_fraction = district_heat_info["urban fraction"]

    # NB: must add costs of central heating afterwards (EUR 400 / kWpeak, 50a, 1% FOM from Fraunhofer ISE)

    # exogenously reduce space heat demand
    if options["reduce_space_heat_exogenously"]:
        dE = get(options["reduce_space_heat_exogenously_factor"], investment_year)
        logger.info(f"Assumed space heat reduction of {dE:.2%}")
        for sector in sectors:
            heat_demand[sector + " space"] = (1 - dE) * heat_demand[sector + " space"]

    if options["solar_thermal"]:
        solar_thermal = (
            xr.open_dataarray(snakemake.input.solar_thermal_total)
            .to_pandas()
            .reindex(index=n.snapshots)
        )
        # 1e3 converts from W/m^2 to MW/(1000m^2) = kW/m^2
        solar_thermal = options["solar_cf_correction"] * solar_thermal / 1e3

    for (
        heat_system
    ) in (
        HeatSystem
    ):  # this loops through all heat systems defined in _entities.HeatSystem
        overdim_factor = options["overdimension_heat_generators"][
            heat_system.central_or_decentral
        ]
        if heat_system == HeatSystem.URBAN_CENTRAL:
            nodes = dist_fraction.index[dist_fraction > 0]
        else:
            nodes = pop_layout.index

        n.add("Carrier", f"{heat_system} heat")

        n.add(
            "Bus",
            nodes + f" {heat_system.value} heat",
            location=nodes,
            carrier=f"{heat_system.value} heat",
            unit="MWh_th",
        )

        if heat_system == HeatSystem.URBAN_CENTRAL and options["central_heat_vent"]:
            n.add(
                "Generator",
                nodes + f" {heat_system} heat vent",
                bus=nodes + f" {heat_system} heat",
                location=nodes,
                carrier=f"{heat_system} heat vent",
                p_nom_extendable=True,
                p_max_pu=0,
                p_min_pu=-1,
                unit="MWh_th",
            )

        ## Add heat load
        factor = heat_system.heat_demand_weighting(
            urban_fraction=urban_fraction[nodes], dist_fraction=dist_fraction[nodes]
        )
        if not heat_system == HeatSystem.URBAN_CENTRAL:
            heat_load = (
                heat_demand[
                    [
                        heat_system.sector.value + " water",
                        heat_system.sector.value + " space",
                    ]
                ]
                .T.groupby(level=1)
                .sum()
                .T[nodes]
                .multiply(factor)
            )

        if heat_system == HeatSystem.URBAN_CENTRAL:
            heat_load = (
                heat_demand.T.groupby(level=1)
                .sum()
                .T[nodes]
                .multiply(
                    factor * (1 + options["district_heating"]["district_heating_loss"])
                )
            )

        n.add(
            "Load",
            nodes,
            suffix=f" {heat_system} heat",
            bus=nodes + f" {heat_system} heat",
            carrier=f"{heat_system} heat",
            p_set=heat_load.loc[n.snapshots],
        )

        ## Add heat pumps
        for heat_source in snakemake.params.heat_pump_sources[
            heat_system.system_type.value
        ]:
            costs_name_heat_pump = heat_system.heat_pump_costs_name(heat_source)
            cop_heat_pump = (
                cop.sel(
                    heat_system=heat_system.system_type.value,
                    heat_source=heat_source,
                    name=nodes,
                )
                .to_pandas()
                .reindex(index=n.snapshots)
                if options["time_dep_hp_cop"]
                else costs.at[costs_name_heat_pump, "efficiency"]
            )

            if heat_source in snakemake.params.heat_utilisation_potentials:
                # get potential
                p_max_source = pd.read_csv(
                    snakemake.input[heat_source],
                    index_col=0,
                ).squeeze()[nodes]

                # add resource
                heat_carrier = f"{heat_system} {heat_source} heat"
                n.add("Carrier", heat_carrier)
                n.add(
                    "Bus",
                    nodes,
                    suffix=f" {heat_carrier}",
                    carrier=heat_carrier,
                )

                costs_name_heat_source = heat_system.heat_source_costs_name(heat_source)
                if heat_source in snakemake.params.direct_utilisation_heat_sources:
                    capital_cost = (
                        costs.at[
                            heat_system.heat_source_costs_name(heat_source), "fixed"
                        ]
                        * overdim_factor
                    )
                    lifetime = costs.at[
                        heat_system.heat_source_costs_name(heat_source), "lifetime"
                    ]
                else:
                    capital_cost = 0.0
                    lifetime = np.inf
                n.add(
                    "Generator",
                    nodes,
                    suffix=f" {heat_carrier}",
                    bus=nodes + f" {heat_carrier}",
                    carrier=heat_carrier,
                    p_nom_extendable=True,
                    capital_cost=capital_cost,
                    lifetime=lifetime,
                    p_nom_max=p_max_source,
                )

                # add heat pump converting source heat + electricity to urban central heat
                n.add(
                    "Link",
                    nodes,
                    suffix=f" {heat_system} {heat_source} heat pump",
                    bus0=nodes,
                    bus1=nodes + f" {heat_carrier}",
                    bus2=nodes + f" {heat_system} heat",
                    carrier=f"{heat_system} {heat_source} heat pump",
                    efficiency=-(cop_heat_pump - 1),
                    efficiency2=cop_heat_pump,
                    capital_cost=costs.at[costs_name_heat_pump, "efficiency"]
                    * costs.at[costs_name_heat_pump, "fixed"]
                    * overdim_factor,
                    p_nom_extendable=True,
                    lifetime=costs.at[costs_name_heat_pump, "lifetime"],
                )

                if heat_source in snakemake.params.direct_utilisation_heat_sources:
                    # 1 if source temperature exceeds forward temperature, 0 otherwise:
                    efficiency_direct_utilisation = (
                        direct_heat_source_utilisation_profile.sel(
                            heat_source=heat_source,
                            name=nodes,
                        )
                        .to_pandas()
                        .reindex(index=n.snapshots)
                    )
                    # add link for direct usage of heat source when source temperature exceeds forward temperature
                    n.add(
                        "Link",
                        nodes,
                        suffix=f" {heat_system} {heat_source} heat direct utilisation",
                        bus0=nodes + f" {heat_carrier}",
                        bus1=nodes + f" {heat_system} heat",
                        efficiency=efficiency_direct_utilisation,
                        carrier=f"{heat_system} {heat_source} heat direct utilisation",
                        p_nom_extendable=True,
                    )
            else:
                n.add(
                    "Link",
                    nodes,
                    suffix=f" {heat_system} {heat_source} heat pump",
                    bus0=nodes,
                    bus1=nodes + f" {heat_system} heat",
                    carrier=f"{heat_system} {heat_source} heat pump",
                    efficiency=cop_heat_pump,
                    capital_cost=costs.at[costs_name_heat_pump, "efficiency"]
                    * costs.at[costs_name_heat_pump, "fixed"]
                    * overdim_factor,
                    p_nom_extendable=True,
                    lifetime=costs.at[costs_name_heat_pump, "lifetime"],
                )

        if options["tes"]:
            n.add("Carrier", f"{heat_system} water tanks")

            n.add(
                "Bus",
                nodes + f" {heat_system} water tanks",
                location=nodes,
                carrier=f"{heat_system} water tanks",
                unit="MWh_th",
            )

            n.add(
                "Link",
                nodes + f" {heat_system} water tanks charger",
                bus0=nodes + f" {heat_system} heat",
                bus1=nodes + f" {heat_system} water tanks",
                efficiency=costs.at["water tank charger", "efficiency"],
                carrier=f"{heat_system} water tanks charger",
                p_nom_extendable=True,
            )

            n.add(
                "Link",
                nodes + f" {heat_system} water tanks discharger",
                bus0=nodes + f" {heat_system} water tanks",
                bus1=nodes + f" {heat_system} heat",
                carrier=f"{heat_system} water tanks discharger",
                efficiency=costs.at["water tank discharger", "efficiency"],
                p_nom_extendable=True,
            )

            tes_time_constant_days = options["tes_tau"][
                heat_system.central_or_decentral
            ]

            n.add(
                "Store",
                nodes + f" {heat_system} water tanks",
                bus=nodes + f" {heat_system} water tanks",
                e_cyclic=True,
                e_nom_extendable=True,
                carrier=f"{heat_system} water tanks",
                standing_loss=1 - np.exp(-1 / 24 / tes_time_constant_days),
                capital_cost=costs.at[
                    heat_system.central_or_decentral + " water tank storage", "fixed"
                ],
                lifetime=costs.at[
                    heat_system.central_or_decentral + " water tank storage", "lifetime"
                ],
            )

        if options["resistive_heaters"]:
            key = f"{heat_system.central_or_decentral} resistive heater"

            n.add(
                "Link",
                nodes + f" {heat_system} resistive heater",
                bus0=nodes,
                bus1=nodes + f" {heat_system} heat",
                carrier=f"{heat_system} resistive heater",
                efficiency=costs.at[key, "efficiency"],
                capital_cost=costs.at[key, "efficiency"]
                * costs.at[key, "fixed"]
                * overdim_factor,
                p_nom_extendable=True,
                lifetime=costs.at[key, "lifetime"],
            )

        if options["boilers"]:
            key = f"{heat_system.central_or_decentral} gas boiler"

            n.add(
                "Link",
                nodes + f" {heat_system} gas boiler",
                p_nom_extendable=True,
                bus0=spatial.gas.df.loc[nodes, "nodes"].values,
                bus1=nodes + f" {heat_system} heat",
                bus2="co2 atmosphere",
                carrier=f"{heat_system} gas boiler",
                efficiency=costs.at[key, "efficiency"],
                efficiency2=costs.at["gas", "CO2 intensity"],
                capital_cost=costs.at[key, "efficiency"]
                * costs.at[key, "fixed"]
                * overdim_factor,
                lifetime=costs.at[key, "lifetime"],
            )

        if options["solar_thermal"]:
            n.add("Carrier", f"{heat_system} solar thermal")

            n.add(
                "Generator",
                nodes,
                suffix=f" {heat_system} solar thermal collector",
                bus=nodes + f" {heat_system} heat",
                carrier=f"{heat_system} solar thermal",
                p_nom_extendable=True,
                capital_cost=costs.at[
                    heat_system.central_or_decentral + " solar thermal", "fixed"
                ]
                * overdim_factor,
                p_max_pu=solar_thermal[nodes],
                lifetime=costs.at[
                    heat_system.central_or_decentral + " solar thermal", "lifetime"
                ],
            )

        if options["chp"]["enable"] and heat_system == HeatSystem.URBAN_CENTRAL:
            # add non-biomass CHP; biomass CHP is added in biomass section
            for fuel in options["chp"]["fuel"]:
                if fuel == "solid biomass":
                    # Solid biomass CHP is added in add_biomass
                    continue
                fuel_nodes = getattr(spatial, fuel).df
                n.add(
                    "Link",
                    nodes + f" urban central {fuel} CHP",
                    bus0=fuel_nodes.loc[nodes, "nodes"].values,
                    bus1=nodes,
                    bus2=nodes + " urban central heat",
                    bus3="co2 atmosphere",
                    carrier="urban central CHP",
                    p_nom_extendable=True,
                    capital_cost=costs.at["central gas CHP", "fixed"]
                    * costs.at["central gas CHP", "efficiency"],
                    marginal_cost=costs.at["central gas CHP", "VOM"],
                    efficiency=costs.at["central gas CHP", "efficiency"],
                    efficiency2=costs.at["central gas CHP", "efficiency"]
                    / costs.at["central gas CHP", "c_b"],
                    efficiency3=costs.at[fuel, "CO2 intensity"],
                    lifetime=costs.at["central gas CHP", "lifetime"],
                )

                n.add(
                    "Link",
                    nodes + f" urban central {fuel} CHP CC",
                    bus0=fuel_nodes.loc[nodes, "nodes"].values,
                    bus1=nodes,
                    bus2=nodes + " urban central heat",
                    bus3="co2 atmosphere",
                    bus4=spatial.co2.df.loc[nodes, "nodes"].values,
                    carrier="urban central CHP CC",
                    p_nom_extendable=True,
                    capital_cost=costs.at["central gas CHP", "fixed"]
                    * costs.at["central gas CHP", "efficiency"]
                    + costs.at["biomass CHP capture", "fixed"]
                    * costs.at[fuel, "CO2 intensity"],
                    marginal_cost=costs.at["central gas CHP", "VOM"],
                    efficiency=costs.at["central gas CHP", "efficiency"]
                    - costs.at[fuel, "CO2 intensity"]
                    * (
                        costs.at["biomass CHP capture", "electricity-input"]
                        + costs.at[
                            "biomass CHP capture", "compression-electricity-input"
                        ]
                    ),
                    efficiency2=costs.at["central gas CHP", "efficiency"]
                    / costs.at["central gas CHP", "c_b"]
                    + costs.at[fuel, "CO2 intensity"]
                    * (
                        costs.at["biomass CHP capture", "heat-output"]
                        + costs.at["biomass CHP capture", "compression-heat-output"]
                        - costs.at["biomass CHP capture", "heat-input"]
                    ),
                    efficiency3=costs.at[fuel, "CO2 intensity"]
                    * (1 - costs.at["biomass CHP capture", "capture_rate"]),
                    efficiency4=costs.at[fuel, "CO2 intensity"]
                    * costs.at["biomass CHP capture", "capture_rate"],
                    lifetime=costs.at["central gas CHP", "lifetime"],
                )

        if (
            options["chp"]["enable"]
            and options["chp"]["micro_chp"]
            and heat_system.value != "urban central"
        ):
            n.add(
                "Link",
                nodes + f" {heat_system} micro gas CHP",
                p_nom_extendable=True,
                bus0=spatial.gas.df.loc[nodes, "nodes"].values,
                bus1=nodes,
                bus2=nodes + f" {heat_system} heat",
                bus3="co2 atmosphere",
                carrier=heat_system.value + " micro gas CHP",
                efficiency=costs.at["micro CHP", "efficiency"],
                efficiency2=costs.at["micro CHP", "efficiency-heat"],
                efficiency3=costs.at["gas", "CO2 intensity"],
                capital_cost=costs.at["micro CHP", "fixed"],
                lifetime=costs.at["micro CHP", "lifetime"],
            )

    if options["retrofitting"]["retro_endogen"]:
        logger.info("Add retrofitting endogenously")

        # retrofitting data 'retro_data' with 'costs' [EUR/m^2] and heat
        # demand 'dE' [per unit of original heat demand] for each country and
        # different retrofitting strengths [additional insulation thickness in m]
        retro_data = pd.read_csv(
            snakemake.input.retro_cost,
            index_col=[0, 1],
            skipinitialspace=True,
            header=[0, 1],
        )
        # heated floor area [10^6 * m^2] per country
        floor_area = pd.read_csv(snakemake.input.floor_area, index_col=[0, 1])

        n.add("Carrier", "retrofitting")

        # share of space heat demand 'w_space' of total heat demand
        w_space = {}
        for sector in sectors:
            w_space[sector] = heat_demand[sector + " space"] / (
                heat_demand[sector + " space"] + heat_demand[sector + " water"]
            )
        w_space["tot"] = (
            heat_demand["services space"] + heat_demand["residential space"]
        ) / heat_demand.T.groupby(level=[1]).sum().T

        for name in n.loads[
            n.loads.carrier.isin([x + " heat" for x in HeatSystem])
        ].index:
            node = n.buses.loc[name, "location"]
            ct = pop_layout.loc[node, "ct"]

            # weighting 'f' depending on the size of the population at the node
            if "urban central" in name:
                f = dist_fraction[node]
            elif "urban decentral" in name:
                f = urban_fraction[node] - dist_fraction[node]
            else:
                f = 1 - urban_fraction[node]
            if f == 0:
                continue
            # get sector name ("residential"/"services"/or both "tot" for urban central)
            if "urban central" in name:
                sec = "tot"
            if "residential" in name:
                sec = "residential"
            if "services" in name:
                sec = "services"

            # get floor aread at node and region (urban/rural) in m^2
            floor_area_node = (
                pop_layout.loc[node].fraction * floor_area.loc[ct, "value"] * 10**6
            ).loc[sec] * f
            # total heat demand at node [MWh]
            demand = n.loads_t.p_set[name]

            # space heat demand at node [MWh]
            space_heat_demand = demand * w_space[sec][node]
            # normed time profile of space heat demand 'space_pu' (values between 0-1),
            # p_max_pu/p_min_pu of retrofitting generators
            space_pu = (
                (space_heat_demand / space_heat_demand.max())
                .to_frame(name=node)
                .fillna(0)
            )

            # minimum heat demand 'dE' after retrofitting in units of original heat demand (values between 0-1)
            dE = retro_data.loc[(ct, sec), ("dE")]
            # get additional energy savings 'dE_diff' between the different retrofitting strengths/generators at one node
            dE_diff = abs(dE.diff()).fillna(1 - dE.iloc[0])
            # convert costs Euro/m^2 -> Euro/MWh
            capital_cost = (
                retro_data.loc[(ct, sec), ("cost")]
                * floor_area_node
                / ((1 - dE) * space_heat_demand.max())
            )
            if space_heat_demand.max() == 0:
                capital_cost = capital_cost.apply(lambda b: 0 if b == np.inf else b)

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
            space_pu = space_pu.reindex(index=heat_demand.index).ffill()

            # add for each retrofitting strength a generator with heat generation profile following the profile of the heat demand
            for strength in strengths:
                node_name = " ".join(name.split(" ")[2::])
                n.add(
                    "Generator",
                    [node],
                    suffix=" retrofitting " + strength + " " + node_name,
                    bus=name,
                    carrier="retrofitting",
                    p_nom_extendable=True,
                    p_nom_max=dE_diff[strength]
                    * space_heat_demand.max(),  # maximum energy savings for this renovation strength
                    p_max_pu=space_pu,
                    p_min_pu=space_pu,
                    country=ct,
                    capital_cost=capital_cost[strength]
                    * options["retrofitting"]["cost_factor"],
                )


def add_methanol(n, costs):
    methanol_options = options["methanol"]
    if not any(methanol_options.values()):
        return

    logger.info("Add methanol")
    add_carrier_buses(n, "methanol")

    if options["biomass"]:
        if methanol_options["biomass_to_methanol"]:
            add_biomass_to_methanol(n, costs)

        if methanol_options["biomass_to_methanol"]:
            add_biomass_to_methanol_cc(n, costs)

    if methanol_options["methanol_to_power"]:
        add_methanol_to_power(n, costs, types=methanol_options["methanol_to_power"])

    if methanol_options["methanol_reforming"]:
        add_methanol_reforming(n, costs)

    if methanol_options["methanol_reforming_cc"]:
        add_methanol_reforming_cc(n, costs)


def add_biomass(n, costs):
    logger.info("Add biomass")

    biomass_potentials = pd.read_csv(snakemake.input.biomass_potentials, index_col=0)

    # need to aggregate potentials if gas not nodally resolved
    if options["gas_network"]:
        biogas_potentials_spatial = biomass_potentials["biogas"].rename(
            index=lambda x: x + " biogas"
        )
        unsustainable_biogas_potentials_spatial = biomass_potentials[
            "unsustainable biogas"
        ].rename(index=lambda x: x + " biogas")
    else:
        biogas_potentials_spatial = biomass_potentials["biogas"].sum()
        unsustainable_biogas_potentials_spatial = biomass_potentials[
            "unsustainable biogas"
        ].sum()

    if options.get("biomass_spatial", options["biomass_transport"]):
        solid_biomass_potentials_spatial = biomass_potentials["solid biomass"].rename(
            index=lambda x: x + " solid biomass"
        )
        msw_biomass_potentials_spatial = biomass_potentials[
            "municipal solid waste"
        ].rename(index=lambda x: x + " municipal solid waste")
        unsustainable_solid_biomass_potentials_spatial = biomass_potentials[
            "unsustainable solid biomass"
        ].rename(index=lambda x: x + " unsustainable solid biomass")
        unsustainable_liquid_biofuel_potentials_spatial = biomass_potentials[
            "unsustainable bioliquids"
        ].rename(index=lambda x: x + " unsustainable bioliquids")

    else:
        solid_biomass_potentials_spatial = biomass_potentials["solid biomass"].sum()
        msw_biomass_potentials_spatial = biomass_potentials[
            "municipal solid waste"
        ].sum()
        unsustainable_solid_biomass_potentials_spatial = biomass_potentials[
            "unsustainable solid biomass"
        ].sum()
        unsustainable_liquid_biofuel_potentials_spatial = biomass_potentials[
            "unsustainable bioliquids"
        ].sum()

    n.add("Carrier", "biogas")
    n.add("Carrier", "solid biomass")

    if (
        options["municipal_solid_waste"]
        and not options["industry"]
        and not (cf_industry["waste_to_energy"] or cf_industry["waste_to_energy_cc"])
    ):
        logger.warning(
            "Flag municipal_solid_waste can be only used with industry "
            "sector waste to energy."
            "Setting municipal_solid_waste=False."
        )
        options["municipal_solid_waste"] = False

    if options["municipal_solid_waste"]:
        n.add("Carrier", "municipal solid waste")

        n.add(
            "Bus",
            spatial.msw.nodes,
            location=spatial.msw.locations,
            carrier="municipal solid waste",
        )

        n.add(
            "Generator",
            spatial.msw.nodes,
            bus=spatial.msw.nodes,
            carrier="municipal solid waste",
            p_nom=msw_biomass_potentials_spatial,
            marginal_cost=0,  # costs.at["municipal solid waste", "fuel"],
            e_sum_min=msw_biomass_potentials_spatial,
            e_sum_max=msw_biomass_potentials_spatial,
        )

    n.add(
        "Bus",
        spatial.gas.biogas,
        location=spatial.gas.locations,
        carrier="biogas",
        unit="MWh_LHV",
    )

    n.add(
        "Bus",
        spatial.biomass.nodes,
        location=spatial.biomass.locations,
        carrier="solid biomass",
        unit="MWh_LHV",
    )

    n.add(
        "Generator",
        spatial.gas.biogas,
        bus=spatial.gas.biogas,
        carrier="biogas",
        p_nom=biogas_potentials_spatial,
        marginal_cost=costs.at["biogas", "fuel"],
        e_sum_min=0,
        e_sum_max=biogas_potentials_spatial,
    )

    n.add(
        "Generator",
        spatial.biomass.nodes,
        bus=spatial.biomass.nodes,
        carrier="solid biomass",
        p_nom=solid_biomass_potentials_spatial,
        marginal_cost=costs.at["solid biomass", "fuel"],
        e_sum_min=0,
        e_sum_max=solid_biomass_potentials_spatial,
    )

    if options["solid_biomass_import"].get("enable", False):
        biomass_import_price = options["solid_biomass_import"]["price"]
        # convert TWh in MWh
        biomass_import_max_amount = options["solid_biomass_import"]["max_amount"] * 1e6
        biomass_import_upstream_emissions = options["solid_biomass_import"][
            "upstream_emissions_factor"
        ]

        logger.info(
            "Adding biomass import with cost %.2f EUR/MWh, a limit of %.2f TWh, and embedded emissions of %.2f%%",
            biomass_import_price,
            options["solid_biomass_import"]["max_amount"],
            biomass_import_upstream_emissions * 100,
        )

        n.add("Carrier", "solid biomass import")

        n.add(
            "Bus",
            ["EU solid biomass import"],
            location="EU",
            carrier="solid biomass import",
        )

        n.add(
            "Store",
            ["solid biomass import"],
            bus=["EU solid biomass import"],
            carrier="solid biomass import",
            e_nom=biomass_import_max_amount,
            marginal_cost=biomass_import_price,
            e_initial=biomass_import_max_amount,
        )

        n.add(
            "Link",
            spatial.biomass.nodes,
            suffix=" solid biomass import",
            bus0=["EU solid biomass import"],
            bus1=spatial.biomass.nodes,
            bus2="co2 atmosphere",
            carrier="solid biomass import",
            efficiency=1.0,
            efficiency2=biomass_import_upstream_emissions
            * costs.at["solid biomass", "CO2 intensity"],
            p_nom_extendable=True,
        )

    if biomass_potentials.filter(like="unsustainable").sum().sum() > 0:
        n.add(
            "Generator",
            spatial.gas.biogas,
            suffix=" unsustainable",
            bus=spatial.gas.biogas,
            carrier="unsustainable biogas",
            p_nom=unsustainable_biogas_potentials_spatial,
            p_nom_extendable=False,
            marginal_cost=costs.at["biogas", "fuel"],
            e_sum_min=unsustainable_biogas_potentials_spatial,
            e_sum_max=unsustainable_biogas_potentials_spatial,
        )

        n.add(
            "Generator",
            spatial.biomass.nodes_unsustainable,
            bus=spatial.biomass.nodes,
            carrier="unsustainable solid biomass",
            p_nom=unsustainable_solid_biomass_potentials_spatial,
            p_nom_extendable=False,
            marginal_cost=costs.at["fuelwood", "fuel"],
            e_sum_min=unsustainable_solid_biomass_potentials_spatial,
            e_sum_max=unsustainable_solid_biomass_potentials_spatial,
        )

        n.add(
            "Bus",
            spatial.biomass.bioliquids,
            location=spatial.biomass.locations,
            carrier="unsustainable bioliquids",
            unit="MWh_LHV",
        )

        n.add(
            "Generator",
            spatial.biomass.bioliquids,
            bus=spatial.biomass.bioliquids,
            carrier="unsustainable bioliquids",
            p_nom=unsustainable_liquid_biofuel_potentials_spatial,
            p_nom_extendable=False,
            marginal_cost=costs.at["biodiesel crops", "fuel"],
            e_sum_min=unsustainable_liquid_biofuel_potentials_spatial,
            e_sum_max=unsustainable_liquid_biofuel_potentials_spatial,
        )

        add_carrier_buses(n, "oil")

        n.add(
            "Link",
            spatial.biomass.bioliquids,
            bus0=spatial.biomass.bioliquids,
            bus1=spatial.oil.nodes,
            bus2="co2 atmosphere",
            carrier="unsustainable bioliquids",
            efficiency=1,
            efficiency2=-costs.at["oil", "CO2 intensity"],
            p_nom=unsustainable_liquid_biofuel_potentials_spatial,
            marginal_cost=costs.at["BtL", "VOM"],
        )

    n.add(
        "Link",
        spatial.gas.biogas_to_gas,
        bus0=spatial.gas.biogas,
        bus1=spatial.gas.nodes,
        bus2="co2 atmosphere",
        carrier="biogas to gas",
        capital_cost=costs.at["biogas", "fixed"]
        + costs.at["biogas upgrading", "fixed"],
        marginal_cost=costs.at["biogas upgrading", "VOM"],
        efficiency=costs.at["biogas", "efficiency"],
        efficiency2=-costs.at["gas", "CO2 intensity"],
        p_nom_extendable=True,
        lifetime=costs.at["biogas", "lifetime"],
    )

    if options["biogas_upgrading_cc"]:
        # Assuming for costs that the CO2 from upgrading is pure, such as in amine scrubbing. I.e., with and without CC is
        # equivalent. Adding biomass CHP capture because biogas is often small-scale and decentral so further
        # from e.g. CO2 grid or buyers. This is a proxy for the added cost for e.g. a raw biogas pipeline to a central upgrading facility
        n.add(
            "Link",
            spatial.gas.biogas_to_gas_cc,
            bus0=spatial.gas.biogas,
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            bus3="co2 atmosphere",
            carrier="biogas to gas CC",
            capital_cost=costs.at["biogas CC", "fixed"]
            + costs.at["biogas upgrading", "fixed"]
            + costs.at["biomass CHP capture", "fixed"]
            * costs.at["biogas CC", "CO2 stored"],
            marginal_cost=costs.at["biogas CC", "VOM"]
            + costs.at["biogas upgrading", "VOM"],
            efficiency=costs.at["biogas CC", "efficiency"],
            efficiency2=costs.at["biogas CC", "CO2 stored"]
            * costs.at["biogas CC", "capture rate"],
            efficiency3=-costs.at["gas", "CO2 intensity"]
            - costs.at["biogas CC", "CO2 stored"]
            * costs.at["biogas CC", "capture rate"],
            p_nom_extendable=True,
            lifetime=costs.at["biogas CC", "lifetime"],
        )

    if options["biomass_transport"]:
        # add biomass transport
        transport_costs = pd.read_csv(
            snakemake.input.biomass_transport_costs, index_col=0
        )
        transport_costs = transport_costs.squeeze()
        biomass_transport = create_network_topology(
            n, "biomass transport ", bidirectional=False
        )

        # costs
        bus0_costs = biomass_transport.bus0.apply(lambda x: transport_costs[x[:2]])
        bus1_costs = biomass_transport.bus1.apply(lambda x: transport_costs[x[:2]])
        biomass_transport["costs"] = pd.concat([bus0_costs, bus1_costs], axis=1).mean(
            axis=1
        )

        n.add(
            "Link",
            biomass_transport.index,
            bus0=biomass_transport.bus0 + " solid biomass",
            bus1=biomass_transport.bus1 + " solid biomass",
            p_nom_extendable=False,
            p_nom=5e4,
            length=biomass_transport.length.values,
            marginal_cost=biomass_transport.costs * biomass_transport.length.values,
            carrier="solid biomass transport",
        )

        if options["municipal_solid_waste"]:
            n.add(
                "Link",
                biomass_transport.index + " municipal solid waste",
                bus0=biomass_transport.bus0.values + " municipal solid waste",
                bus1=biomass_transport.bus1.values + " municipal solid waste",
                p_nom_extendable=False,
                p_nom=5e4,
                length=biomass_transport.length.values,
                marginal_cost=(
                    biomass_transport.costs * biomass_transport.length
                ).values,
                carrier="municipal solid waste transport",
            )

    elif options["biomass_spatial"]:
        # add artificial biomass generators at nodes which include transport costs
        transport_costs = pd.read_csv(
            snakemake.input.biomass_transport_costs, index_col=0
        )
        transport_costs = transport_costs.squeeze()
        bus_transport_costs = spatial.biomass.nodes.to_series().apply(
            lambda x: transport_costs[x[:2]]
        )
        average_distance = 200  # km #TODO: validate this assumption

        n.add(
            "Generator",
            spatial.biomass.nodes,
            suffix=" transported",
            bus=spatial.biomass.nodes,
            carrier="solid biomass",
            p_nom=10000,
            marginal_cost=costs.at["solid biomass", "fuel"]
            + bus_transport_costs * average_distance,
        )
        n.add(
            "GlobalConstraint",
            "biomass limit",
            carrier_attribute="solid biomass",
            sense="<=",
            constant=biomass_potentials["solid biomass"].sum(),
            type="operational_limit",
        )
        if biomass_potentials["unsustainable solid biomass"].sum() > 0:
            n.add(
                "Generator",
                spatial.biomass.nodes_unsustainable,
                suffix=" transported",
                bus=spatial.biomass.nodes,
                carrier="unsustainable solid biomass",
                p_nom=10000,
                marginal_cost=costs.at["fuelwood", "fuel"]
                + bus_transport_costs.rename(
                    dict(
                        zip(spatial.biomass.nodes, spatial.biomass.nodes_unsustainable)
                    )
                )
                * average_distance,
            )
            # Set e_sum_min to 0 to allow for the faux biomass transport
            n.generators.loc[
                n.generators.carrier == "unsustainable solid biomass", "e_sum_min"
            ] = 0

            n.add(
                "GlobalConstraint",
                "unsustainable biomass limit",
                carrier_attribute="unsustainable solid biomass",
                sense="==",
                constant=biomass_potentials["unsustainable solid biomass"].sum(),
                type="operational_limit",
            )

        if options["municipal_solid_waste"]:
            # Add municipal solid waste
            n.add(
                "Generator",
                spatial.msw.nodes,
                suffix=" transported",
                bus=spatial.msw.nodes,
                carrier="municipal solid waste",
                p_nom=10000,
                marginal_cost=0  # costs.at["municipal solid waste", "fuel"]
                + bus_transport_costs.rename(
                    dict(zip(spatial.biomass.nodes, spatial.msw.nodes))
                )
                * average_distance,
            )
            n.generators.loc[
                n.generators.carrier == "municipal solid waste", "e_sum_min"
            ] = 0
            n.add(
                "GlobalConstraint",
                "msw limit",
                carrier_attribute="municipal solid waste",
                sense="==",
                constant=biomass_potentials["municipal solid waste"].sum(),
                type="operational_limit",
            )

    # AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if (
        not urban_central.empty
        and options["chp"]["enable"]
        and ("solid biomass" in options["chp"]["fuel"])
    ):
        urban_central = urban_central.str[: -len(" urban central heat")]

        key = "central solid biomass CHP"

        n.add(
            "Link",
            urban_central + " urban central solid biomass CHP",
            bus0=spatial.biomass.df.loc[urban_central, "nodes"].values,
            bus1=urban_central,
            bus2=urban_central + " urban central heat",
            carrier="urban central solid biomass CHP",
            p_nom_extendable=True,
            capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"],
            marginal_cost=costs.at[key, "VOM"],
            efficiency=costs.at[key, "efficiency"],
            efficiency2=costs.at[key, "efficiency-heat"],
            lifetime=costs.at[key, "lifetime"],
        )

        n.add(
            "Link",
            urban_central + " urban central solid biomass CHP CC",
            bus0=spatial.biomass.df.loc[urban_central, "nodes"].values,
            bus1=urban_central,
            bus2=urban_central + " urban central heat",
            bus3="co2 atmosphere",
            bus4=spatial.co2.df.loc[urban_central, "nodes"].values,
            carrier="urban central solid biomass CHP CC",
            p_nom_extendable=True,
            capital_cost=costs.at[key + " CC", "fixed"]
            * costs.at[key + " CC", "efficiency"]
            + costs.at["biomass CHP capture", "fixed"]
            * costs.at["solid biomass", "CO2 intensity"],
            marginal_cost=costs.at[key + " CC", "VOM"],
            efficiency=costs.at[key + " CC", "efficiency"]
            - costs.at["solid biomass", "CO2 intensity"]
            * (
                costs.at["biomass CHP capture", "electricity-input"]
                + costs.at["biomass CHP capture", "compression-electricity-input"]
            ),
            efficiency2=costs.at[key + " CC", "efficiency-heat"],
            efficiency3=-costs.at["solid biomass", "CO2 intensity"]
            * costs.at["biomass CHP capture", "capture_rate"],
            efficiency4=costs.at["solid biomass", "CO2 intensity"]
            * costs.at["biomass CHP capture", "capture_rate"],
            lifetime=costs.at[key + " CC", "lifetime"],
        )

    if options["biomass_boiler"]:
        # TODO: Add surcharge for pellets
        nodes = pop_layout.index
        for name in [
            "residential rural",
            "services rural",
            "residential urban decentral",
            "services urban decentral",
        ]:
            n.add(
                "Link",
                nodes + f" {name} biomass boiler",
                p_nom_extendable=True,
                bus0=spatial.biomass.df.loc[nodes, "nodes"].values,
                bus1=nodes + f" {name} heat",
                carrier=name + " biomass boiler",
                efficiency=costs.at["biomass boiler", "efficiency"],
                capital_cost=costs.at["biomass boiler", "efficiency"]
                * costs.at["biomass boiler", "fixed"]
                * options["overdimension_heat_generators"][
                    HeatSystem(name).central_or_decentral
                ],
                marginal_cost=costs.at["biomass boiler", "pelletizing cost"],
                lifetime=costs.at["biomass boiler", "lifetime"],
            )

    # Solid biomass to liquid fuel
    if options["biomass_to_liquid"]:
        add_carrier_buses(n, "oil")
        n.add(
            "Link",
            spatial.biomass.nodes,
            suffix=" biomass to liquid",
            bus0=spatial.biomass.nodes,
            bus1=spatial.oil.nodes,
            bus2="co2 atmosphere",
            carrier="biomass to liquid",
            lifetime=costs.at["BtL", "lifetime"],
            efficiency=costs.at["BtL", "efficiency"],
            efficiency2=-costs.at["solid biomass", "CO2 intensity"]
            + costs.at["BtL", "CO2 stored"],
            p_nom_extendable=True,
            capital_cost=costs.at["BtL", "fixed"] * costs.at["BtL", "efficiency"],
            marginal_cost=costs.at["BtL", "VOM"] * costs.at["BtL", "efficiency"],
        )

    # Solid biomass to liquid fuel with carbon capture
    if options["biomass_to_liquid_cc"]:
        # Assuming that acid gas removal (incl. CO2) from syngas i performed with Rectisol
        # process (Methanol) and that electricity demand for this is included in the base process
        n.add(
            "Link",
            spatial.biomass.nodes,
            suffix=" biomass to liquid CC",
            bus0=spatial.biomass.nodes,
            bus1=spatial.oil.nodes,
            bus2="co2 atmosphere",
            bus3=spatial.co2.nodes,
            carrier="biomass to liquid CC",
            lifetime=costs.at["BtL", "lifetime"],
            efficiency=costs.at["BtL", "efficiency"],
            efficiency2=-costs.at["solid biomass", "CO2 intensity"]
            + costs.at["BtL", "CO2 stored"] * (1 - costs.at["BtL", "capture rate"]),
            efficiency3=costs.at["BtL", "CO2 stored"] * costs.at["BtL", "capture rate"],
            p_nom_extendable=True,
            capital_cost=costs.at["BtL", "fixed"] * costs.at["BtL", "efficiency"]
            + costs.at["biomass CHP capture", "fixed"] * costs.at["BtL", "CO2 stored"],
            marginal_cost=costs.at["BtL", "VOM"] * costs.at["BtL", "efficiency"],
        )

    # Electrobiofuels (BtL with hydrogen addition to make more use of biogenic carbon).
    # Combination of efuels and biomass to liquid, both based on Fischer-Tropsch.
    # Experimental version - use with caution
    if options["electrobiofuels"]:
        add_carrier_buses(n, "oil")
        efuel_scale_factor = costs.at["BtL", "C stored"]
        name = (
            pd.Index(spatial.biomass.nodes)
            + " "
            + pd.Index(spatial.h2.nodes.str.replace(" H2", ""))
        )
        n.add(
            "Link",
            name,
            suffix=" electrobiofuels",
            bus0=spatial.biomass.nodes,
            bus1=spatial.oil.nodes,
            bus2=spatial.h2.nodes,
            bus3="co2 atmosphere",
            carrier="electrobiofuels",
            lifetime=costs.at["electrobiofuels", "lifetime"],
            efficiency=costs.at["electrobiofuels", "efficiency-biomass"],
            efficiency2=-costs.at["electrobiofuels", "efficiency-hydrogen"],
            efficiency3=-costs.at["solid biomass", "CO2 intensity"]
            + costs.at["BtL", "CO2 stored"]
            * (1 - costs.at["Fischer-Tropsch", "capture rate"]),
            p_nom_extendable=True,
            capital_cost=costs.at["BtL", "fixed"] * costs.at["BtL", "efficiency"]
            + efuel_scale_factor
            * costs.at["Fischer-Tropsch", "fixed"]
            * costs.at["Fischer-Tropsch", "efficiency"],
            marginal_cost=costs.at["BtL", "VOM"] * costs.at["BtL", "efficiency"]
            + efuel_scale_factor
            * costs.at["Fischer-Tropsch", "VOM"]
            * costs.at["Fischer-Tropsch", "efficiency"],
        )

    # BioSNG from solid biomass
    if options["biosng"]:
        n.add(
            "Link",
            spatial.biomass.nodes,
            suffix=" solid biomass to gas",
            bus0=spatial.biomass.nodes,
            bus1=spatial.gas.nodes,
            bus3="co2 atmosphere",
            carrier="BioSNG",
            lifetime=costs.at["BioSNG", "lifetime"],
            efficiency=costs.at["BioSNG", "efficiency"],
            efficiency3=-costs.at["solid biomass", "CO2 intensity"]
            + costs.at["BioSNG", "CO2 stored"],
            p_nom_extendable=True,
            capital_cost=costs.at["BioSNG", "fixed"] * costs.at["BioSNG", "efficiency"],
            marginal_cost=costs.at["BioSNG", "VOM"] * costs.at["BioSNG", "efficiency"],
        )

    # BioSNG from solid biomass with carbon capture
    if options["biosng_cc"]:
        # Assuming that acid gas removal (incl. CO2) from syngas i performed with Rectisol
        # process (Methanol) and that electricity demand for this is included in the base process
        n.add(
            "Link",
            spatial.biomass.nodes,
            suffix=" solid biomass to gas CC",
            bus0=spatial.biomass.nodes,
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            bus3="co2 atmosphere",
            carrier="BioSNG CC",
            lifetime=costs.at["BioSNG", "lifetime"],
            efficiency=costs.at["BioSNG", "efficiency"],
            efficiency2=costs.at["BioSNG", "CO2 stored"]
            * costs.at["BioSNG", "capture rate"],
            efficiency3=-costs.at["solid biomass", "CO2 intensity"]
            + costs.at["BioSNG", "CO2 stored"]
            * (1 - costs.at["BioSNG", "capture rate"]),
            p_nom_extendable=True,
            capital_cost=costs.at["BioSNG", "fixed"] * costs.at["BioSNG", "efficiency"]
            + costs.at["biomass CHP capture", "fixed"]
            * costs.at["BioSNG", "CO2 stored"],
            marginal_cost=costs.at["BioSNG", "VOM"] * costs.at["BioSNG", "efficiency"],
        )

    if options["bioH2"]:
        name = (
            pd.Index(spatial.biomass.nodes)
            + " "
            + pd.Index(spatial.h2.nodes.str.replace(" H2", ""))
        )
        n.add(
            "Link",
            name,
            suffix=" solid biomass to hydrogen CC",
            bus0=spatial.biomass.nodes,
            bus1=spatial.h2.nodes,
            bus2=spatial.co2.nodes,
            bus3="co2 atmosphere",
            carrier="solid biomass to hydrogen",
            efficiency=costs.at["solid biomass to hydrogen", "efficiency"],
            efficiency2=costs.at["solid biomass", "CO2 intensity"]
            * options["cc_fraction"],
            efficiency3=-costs.at["solid biomass", "CO2 intensity"]
            * options["cc_fraction"],
            p_nom_extendable=True,
            capital_cost=costs.at["solid biomass to hydrogen", "fixed"]
            * costs.at["solid biomass to hydrogen", "efficiency"]
            + costs.at["biomass CHP capture", "fixed"]
            * costs.at["solid biomass", "CO2 intensity"],
            marginal_cost=0.0,
            lifetime=25,  # TODO: add value to technology-data
        )


def add_industry(n, costs):
    logger.info("Add industrial demand")
    # add oil buses for shipping, aviation and naptha for industry
    add_carrier_buses(n, "oil")
    # add methanol buses for industry
    add_carrier_buses(n, "methanol")

    nodes = pop_layout.index
    nhours = n.snapshot_weightings.generators.sum()
    nyears = nhours / 8760

    # 1e6 to convert TWh to MWh
    industrial_demand = (
        pd.read_csv(snakemake.input.industrial_demand, index_col=0) * 1e6
    ) * nyears

    n.add(
        "Bus",
        spatial.biomass.industry,
        location=spatial.biomass.locations,
        carrier="solid biomass for industry",
        unit="MWh_LHV",
    )

    if options.get("biomass_spatial", options["biomass_transport"]):
        p_set = (
            industrial_demand.loc[spatial.biomass.locations, "solid biomass"].rename(
                index=lambda x: x + " solid biomass for industry"
            )
            / nhours
        )
    else:
        p_set = industrial_demand["solid biomass"].sum() / nhours

    n.add(
        "Load",
        spatial.biomass.industry,
        bus=spatial.biomass.industry,
        carrier="solid biomass for industry",
        p_set=p_set,
    )

    n.add(
        "Link",
        spatial.biomass.industry,
        bus0=spatial.biomass.nodes,
        bus1=spatial.biomass.industry,
        carrier="solid biomass for industry",
        p_nom_extendable=True,
        efficiency=1.0,
    )

    if len(spatial.biomass.industry_cc) <= 1 and len(spatial.co2.nodes) > 1:
        link_names = nodes + " " + spatial.biomass.industry_cc
    else:
        link_names = spatial.biomass.industry_cc

    n.add(
        "Link",
        link_names,
        bus0=spatial.biomass.nodes,
        bus1=spatial.biomass.industry,
        bus2="co2 atmosphere",
        bus3=spatial.co2.nodes,
        carrier="solid biomass for industry CC",
        p_nom_extendable=True,
        capital_cost=costs.at["cement capture", "fixed"]
        * costs.at["solid biomass", "CO2 intensity"],
        efficiency=0.9,  # TODO: make config option
        efficiency2=-costs.at["solid biomass", "CO2 intensity"]
        * costs.at["cement capture", "capture_rate"],
        efficiency3=costs.at["solid biomass", "CO2 intensity"]
        * costs.at["cement capture", "capture_rate"],
        lifetime=costs.at["cement capture", "lifetime"],
    )

    n.add(
        "Bus",
        spatial.gas.industry,
        location=spatial.gas.locations,
        carrier="gas for industry",
        unit="MWh_LHV",
    )

    gas_demand = industrial_demand.loc[nodes, "methane"] / nhours

    if options["gas_network"]:
        spatial_gas_demand = gas_demand.rename(index=lambda x: x + " gas for industry")
    else:
        spatial_gas_demand = gas_demand.sum()

    n.add(
        "Load",
        spatial.gas.industry,
        bus=spatial.gas.industry,
        carrier="gas for industry",
        p_set=spatial_gas_demand,
    )

    n.add(
        "Link",
        spatial.gas.industry,
        bus0=spatial.gas.nodes,
        bus1=spatial.gas.industry,
        bus2="co2 atmosphere",
        carrier="gas for industry",
        p_nom_extendable=True,
        efficiency=1.0,
        efficiency2=costs.at["gas", "CO2 intensity"],
    )

    n.add(
        "Link",
        spatial.gas.industry_cc,
        bus0=spatial.gas.nodes,
        bus1=spatial.gas.industry,
        bus2="co2 atmosphere",
        bus3=spatial.co2.nodes,
        carrier="gas for industry CC",
        p_nom_extendable=True,
        capital_cost=costs.at["cement capture", "fixed"]
        * costs.at["gas", "CO2 intensity"],
        efficiency=0.9,
        efficiency2=costs.at["gas", "CO2 intensity"]
        * (1 - costs.at["cement capture", "capture_rate"]),
        efficiency3=costs.at["gas", "CO2 intensity"]
        * costs.at["cement capture", "capture_rate"],
        lifetime=costs.at["cement capture", "lifetime"],
    )

    n.add(
        "Load",
        nodes,
        suffix=" H2 for industry",
        bus=nodes + " H2",
        carrier="H2 for industry",
        p_set=industrial_demand.loc[nodes, "hydrogen"] / nhours,
    )

    # methanol for industry

    n.add(
        "Bus",
        spatial.methanol.industry,
        carrier="industry methanol",
        location=spatial.methanol.demand_locations,
        unit="MWh_LHV",
    )

    p_set_methanol = (
        industrial_demand["methanol"].rename(lambda x: x + " industry methanol")
        / nhours
    )

    if not options["methanol"]["regional_methanol_demand"]:
        p_set_methanol = p_set_methanol.sum()

    n.add(
        "Load",
        spatial.methanol.industry,
        bus=spatial.methanol.industry,
        carrier="industry methanol",
        p_set=p_set_methanol,
    )

    n.add(
        "Link",
        spatial.methanol.industry,
        bus0=spatial.methanol.nodes,
        bus1=spatial.methanol.industry,
        bus2="co2 atmosphere",
        carrier="industry methanol",
        p_nom_extendable=True,
        efficiency2=1 / options["MWh_MeOH_per_tCO2"],
        # CO2 intensity methanol based on stoichiometric calculation with 22.7 GJ/t methanol (32 g/mol), CO2 (44 g/mol), 277.78 MWh/TJ = 0.218 t/MWh
    )

    n.add(
        "Link",
        spatial.h2.locations + " methanolisation",
        bus0=spatial.h2.nodes,
        bus1=spatial.methanol.nodes,
        bus2=nodes,
        bus3=spatial.co2.nodes,
        carrier="methanolisation",
        p_nom_extendable=True,
        p_min_pu=options["min_part_load_methanolisation"],
        capital_cost=costs.at["methanolisation", "fixed"]
        * options["MWh_MeOH_per_MWh_H2"],  # EUR/MW_H2/a
        marginal_cost=options["MWh_MeOH_per_MWh_H2"]
        * costs.at["methanolisation", "VOM"],
        lifetime=costs.at["methanolisation", "lifetime"],
        efficiency=options["MWh_MeOH_per_MWh_H2"],
        efficiency2=-options["MWh_MeOH_per_MWh_H2"] / options["MWh_MeOH_per_MWh_e"],
        efficiency3=-options["MWh_MeOH_per_MWh_H2"] / options["MWh_MeOH_per_tCO2"],
    )

    shipping_hydrogen_share = get(options["shipping_hydrogen_share"], investment_year)
    shipping_methanol_share = get(options["shipping_methanol_share"], investment_year)
    shipping_oil_share = get(options["shipping_oil_share"], investment_year)

    total_share = shipping_hydrogen_share + shipping_methanol_share + shipping_oil_share
    if total_share != 1:
        logger.warning(
            f"Total shipping shares sum up to {total_share:.2%}, corresponding to increased or decreased demand assumptions."
        )

    domestic_navigation = pop_weighted_energy_totals.loc[
        nodes, ["total domestic navigation"]
    ].squeeze()
    international_navigation = (
        pd.read_csv(snakemake.input.shipping_demand, index_col=0).squeeze(axis=1)
        * nyears
    )
    all_navigation = domestic_navigation + international_navigation
    p_set = all_navigation * 1e6 / nhours

    if shipping_hydrogen_share:
        oil_efficiency = options.get(
            "shipping_oil_efficiency", options.get("shipping_average_efficiency", 0.4)
        )
        efficiency = oil_efficiency / costs.at["fuel cell", "efficiency"]
        shipping_hydrogen_share = get(
            options["shipping_hydrogen_share"], investment_year
        )

        if options["shipping_hydrogen_liquefaction"]:
            n.add(
                "Bus",
                nodes,
                suffix=" H2 liquid",
                carrier="H2 liquid",
                location=nodes,
                unit="MWh_LHV",
            )

            n.add(
                "Link",
                nodes + " H2 liquefaction",
                bus0=nodes + " H2",
                bus1=nodes + " H2 liquid",
                carrier="H2 liquefaction",
                efficiency=costs.at["H2 liquefaction", "efficiency"],
                capital_cost=costs.at["H2 liquefaction", "fixed"],
                p_nom_extendable=True,
                lifetime=costs.at["H2 liquefaction", "lifetime"],
            )

            shipping_bus = nodes + " H2 liquid"
        else:
            shipping_bus = nodes + " H2"

        efficiency = (
            options["shipping_oil_efficiency"] / costs.at["fuel cell", "efficiency"]
        )
        p_set_hydrogen = shipping_hydrogen_share * p_set * efficiency

        n.add(
            "Load",
            nodes,
            suffix=" H2 for shipping",
            bus=shipping_bus,
            carrier="H2 for shipping",
            p_set=p_set_hydrogen,
        )

    if shipping_methanol_share:
        efficiency = (
            options["shipping_oil_efficiency"] / options["shipping_methanol_efficiency"]
        )

        p_set_methanol_shipping = (
            shipping_methanol_share
            * p_set.rename(lambda x: x + " shipping methanol")
            * efficiency
        )

        if not options["methanol"]["regional_methanol_demand"]:
            p_set_methanol_shipping = p_set_methanol_shipping.sum()

        n.add(
            "Bus",
            spatial.methanol.shipping,
            location=spatial.methanol.demand_locations,
            carrier="shipping methanol",
            unit="MWh_LHV",
        )

        n.add(
            "Load",
            spatial.methanol.shipping,
            bus=spatial.methanol.shipping,
            carrier="shipping methanol",
            p_set=p_set_methanol_shipping,
        )

        n.add(
            "Link",
            spatial.methanol.shipping,
            bus0=spatial.methanol.nodes,
            bus1=spatial.methanol.shipping,
            bus2="co2 atmosphere",
            carrier="shipping methanol",
            p_nom_extendable=True,
            efficiency2=1
            / options[
                "MWh_MeOH_per_tCO2"
            ],  # CO2 intensity methanol based on stoichiometric calculation with 22.7 GJ/t methanol (32 g/mol), CO2 (44 g/mol), 277.78 MWh/TJ = 0.218 t/MWh
        )

    if shipping_oil_share:
        p_set_oil = shipping_oil_share * p_set.rename(lambda x: x + " shipping oil")

        if not options["regional_oil_demand"]:
            p_set_oil = p_set_oil.sum()

        n.add(
            "Bus",
            spatial.oil.shipping,
            location=spatial.oil.demand_locations,
            carrier="shipping oil",
            unit="MWh_LHV",
        )

        n.add(
            "Load",
            spatial.oil.shipping,
            bus=spatial.oil.shipping,
            carrier="shipping oil",
            p_set=p_set_oil,
        )

        n.add(
            "Link",
            spatial.oil.shipping,
            bus0=spatial.oil.nodes,
            bus1=spatial.oil.shipping,
            bus2="co2 atmosphere",
            carrier="shipping oil",
            p_nom_extendable=True,
            efficiency2=costs.at["oil", "CO2 intensity"],
        )

    if options["oil_boilers"]:
        nodes = pop_layout.index

        for heat_system in HeatSystem:
            if not heat_system == HeatSystem.URBAN_CENTRAL:
                n.add(
                    "Link",
                    nodes + f" {heat_system} oil boiler",
                    p_nom_extendable=True,
                    bus0=spatial.oil.nodes,
                    bus1=nodes + f" {heat_system} heat",
                    bus2="co2 atmosphere",
                    carrier=f"{heat_system} oil boiler",
                    efficiency=costs.at["decentral oil boiler", "efficiency"],
                    efficiency2=costs.at["oil", "CO2 intensity"],
                    capital_cost=costs.at["decentral oil boiler", "efficiency"]
                    * costs.at["decentral oil boiler", "fixed"]
                    * options["overdimension_heat_generators"][
                        heat_system.central_or_decentral
                    ],
                    lifetime=costs.at["decentral oil boiler", "lifetime"],
                )

    n.add(
        "Link",
        nodes + " Fischer-Tropsch",
        bus0=nodes + " H2",
        bus1=spatial.oil.nodes,
        bus2=spatial.co2.nodes,
        carrier="Fischer-Tropsch",
        efficiency=costs.at["Fischer-Tropsch", "efficiency"],
        capital_cost=costs.at["Fischer-Tropsch", "fixed"]
        * costs.at["Fischer-Tropsch", "efficiency"],  # EUR/MW_H2/a
        marginal_cost=costs.at["Fischer-Tropsch", "efficiency"]
        * costs.at["Fischer-Tropsch", "VOM"],
        efficiency2=-costs.at["oil", "CO2 intensity"]
        * costs.at["Fischer-Tropsch", "efficiency"],
        p_nom_extendable=True,
        p_min_pu=options["min_part_load_fischer_tropsch"],
        lifetime=costs.at["Fischer-Tropsch", "lifetime"],
    )

    # naphtha
    demand_factor = options["HVC_demand_factor"]
    if demand_factor != 1:
        logger.warning(f"Changing HVC demand by {demand_factor*100-100:+.2f}%.")

    p_set_naphtha = (
        demand_factor
        * industrial_demand.loc[nodes, "naphtha"].rename(
            lambda x: x + " naphtha for industry"
        )
        / nhours
    )

    if not options["regional_oil_demand"]:
        p_set_naphtha = p_set_naphtha.sum()

    n.add(
        "Bus",
        spatial.oil.naphtha,
        location=spatial.oil.demand_locations,
        carrier="naphtha for industry",
        unit="MWh_LHV",
    )

    n.add(
        "Load",
        spatial.oil.naphtha,
        bus=spatial.oil.naphtha,
        carrier="naphtha for industry",
        p_set=p_set_naphtha,
    )

    # some CO2 from naphtha are process emissions from steam cracker
    # rest of CO2 released to atmosphere either in waste-to-energy or decay
    process_co2_per_naphtha = (
        industrial_demand.loc[nodes, "process emission from feedstock"].sum()
        / industrial_demand.loc[nodes, "naphtha"].sum()
    )
    emitted_co2_per_naphtha = costs.at["oil", "CO2 intensity"] - process_co2_per_naphtha

    non_sequestered = 1 - get(
        cf_industry["HVC_environment_sequestration_fraction"],
        investment_year,
    )

    if cf_industry["waste_to_energy"] or cf_industry["waste_to_energy_cc"]:
        non_sequestered_hvc_locations = (
            pd.Index(spatial.oil.demand_locations) + " non-sequestered HVC"
        )

        n.add(
            "Bus",
            non_sequestered_hvc_locations,
            location=spatial.oil.demand_locations,
            carrier="non-sequestered HVC",
            unit="MWh_LHV",
        )

        n.add(
            "Link",
            spatial.oil.naphtha,
            bus0=spatial.oil.nodes,
            bus1=spatial.oil.naphtha,
            bus2=non_sequestered_hvc_locations,
            bus3=spatial.co2.process_emissions,
            carrier="naphtha for industry",
            p_nom_extendable=True,
            efficiency2=non_sequestered
            * emitted_co2_per_naphtha
            / costs.at["oil", "CO2 intensity"],
            efficiency3=process_co2_per_naphtha,
        )

        if options["biomass"] and options["municipal_solid_waste"]:
            n.add(
                "Link",
                spatial.msw.locations,
                bus0=spatial.msw.nodes,
                bus1=non_sequestered_hvc_locations,
                bus2="co2 atmosphere",
                carrier="municipal solid waste",
                p_nom_extendable=True,
                efficiency=1.0,
                efficiency2=-costs.at[
                    "oil", "CO2 intensity"
                ],  # because msw is co2 neutral and will be burned in waste CHP or decomposed as oil
            )

        n.add(
            "Link",
            spatial.oil.demand_locations,
            suffix=" HVC to air",
            bus0=non_sequestered_hvc_locations,
            bus1="co2 atmosphere",
            carrier="HVC to air",
            p_nom_extendable=True,
            efficiency=costs.at["oil", "CO2 intensity"],
        )

        if len(non_sequestered_hvc_locations) == 1:
            waste_source = non_sequestered_hvc_locations[0]
        else:
            waste_source = non_sequestered_hvc_locations

        if cf_industry["waste_to_energy"]:
            urban_central = spatial.nodes + " urban central heat"
            existing_urban_central = n.buses.index[
                n.buses.carrier == "urban central heat"
            ]
            urban_central_nodes = urban_central.map(
                lambda x: x if x in existing_urban_central else ""
            )
            n.add(
                "Link",
                spatial.nodes + " waste CHP",
                bus0=waste_source,
                bus1=spatial.nodes,
                bus2=urban_central_nodes,
                bus3="co2 atmosphere",
                carrier="waste CHP",
                p_nom_extendable=True,
                capital_cost=costs.at["waste CHP", "fixed"]
                * costs.at["waste CHP", "efficiency"],
                marginal_cost=costs.at["waste CHP", "VOM"],
                efficiency=costs.at["waste CHP", "efficiency"],
                efficiency2=costs.at["waste CHP", "efficiency-heat"],
                efficiency3=costs.at["oil", "CO2 intensity"],
                lifetime=costs.at["waste CHP", "lifetime"],
            )

        if cf_industry["waste_to_energy_cc"]:
            n.add(
                "Link",
                spatial.nodes + " waste CHP CC",
                bus0=waste_source,
                bus1=spatial.nodes,
                bus2=urban_central_nodes,
                bus3="co2 atmosphere",
                bus4=spatial.co2.nodes,
                carrier="waste CHP CC",
                p_nom_extendable=True,
                capital_cost=costs.at["waste CHP CC", "fixed"]
                * costs.at["waste CHP CC", "efficiency"]
                + costs.at["biomass CHP capture", "fixed"]
                * costs.at["oil", "CO2 intensity"],
                marginal_cost=costs.at["waste CHP CC", "VOM"],
                efficiency=costs.at["waste CHP CC", "efficiency"],
                efficiency2=costs.at["waste CHP CC", "efficiency-heat"],
                efficiency3=costs.at["oil", "CO2 intensity"]
                * (1 - options["cc_fraction"]),
                efficiency4=costs.at["oil", "CO2 intensity"] * options["cc_fraction"],
                lifetime=costs.at["waste CHP CC", "lifetime"],
            )

    else:
        n.add(
            "Link",
            spatial.oil.naphtha,
            bus0=spatial.oil.nodes,
            bus1=spatial.oil.naphtha,
            bus2="co2 atmosphere",
            bus3=spatial.co2.process_emissions,
            carrier="naphtha for industry",
            p_nom_extendable=True,
            efficiency2=emitted_co2_per_naphtha * non_sequestered,
            efficiency3=process_co2_per_naphtha,
        )

    # aviation
    demand_factor = options["aviation_demand_factor"]
    if demand_factor != 1:
        logger.warning(f"Changing aviation demand by {demand_factor*100-100:+.2f}%.")

    all_aviation = ["total international aviation", "total domestic aviation"]

    p_set = (
        demand_factor
        * pop_weighted_energy_totals.loc[nodes, all_aviation].sum(axis=1)
        * 1e6
        / nhours
    ).rename(lambda x: x + " kerosene for aviation")

    if not options["regional_oil_demand"]:
        p_set = p_set.sum()

    n.add(
        "Bus",
        spatial.oil.kerosene,
        location=spatial.oil.demand_locations,
        carrier="kerosene for aviation",
        unit="MWh_LHV",
    )

    n.add(
        "Load",
        spatial.oil.kerosene,
        bus=spatial.oil.kerosene,
        carrier="kerosene for aviation",
        p_set=p_set,
    )

    n.add(
        "Link",
        spatial.oil.kerosene,
        bus0=spatial.oil.nodes,
        bus1=spatial.oil.kerosene,
        bus2="co2 atmosphere",
        carrier="kerosene for aviation",
        p_nom_extendable=True,
        efficiency2=costs.at["oil", "CO2 intensity"],
    )

    if options["methanol"]["methanol_to_kerosene"]:
        add_methanol_to_kerosene(n, costs)

    # TODO simplify bus expression
    n.add(
        "Load",
        nodes,
        suffix=" low-temperature heat for industry",
        bus=[
            (
                node + " urban central heat"
                if node + " urban central heat" in n.buses.index
                else node + " services urban decentral heat"
            )
            for node in nodes
        ],
        carrier="low-temperature heat for industry",
        p_set=industrial_demand.loc[nodes, "low-temperature heat"] / nhours,
    )

    # remove today's industrial electricity demand by scaling down total electricity demand
    for ct in n.buses.country.dropna().unique():
        # TODO map onto n.bus.country

        loads_i = n.loads.index[
            (n.loads.index.str[:2] == ct) & (n.loads.carrier == "electricity")
        ]
        if n.loads_t.p_set[loads_i].empty:
            continue
        factor = (
            1
            - industrial_demand.loc[loads_i, "current electricity"].sum()
            / n.loads_t.p_set[loads_i].sum().sum()
        )
        n.loads_t.p_set[loads_i] *= factor

    n.add(
        "Load",
        nodes,
        suffix=" industry electricity",
        bus=nodes,
        carrier="industry electricity",
        p_set=industrial_demand.loc[nodes, "electricity"] / nhours,
    )

    n.add(
        "Bus",
        spatial.co2.process_emissions,
        location=spatial.co2.locations,
        carrier="process emissions",
        unit="t_co2",
    )

    if options["co2_spatial"] or options["co2_network"]:
        p_set = (
            -industrial_demand.loc[nodes, "process emission"].rename(
                index=lambda x: x + " process emissions"
            )
            / nhours
        )
    else:
        p_set = -industrial_demand.loc[nodes, "process emission"].sum() / nhours

    n.add(
        "Load",
        spatial.co2.process_emissions,
        bus=spatial.co2.process_emissions,
        carrier="process emissions",
        p_set=p_set,
    )

    n.add(
        "Link",
        spatial.co2.process_emissions,
        bus0=spatial.co2.process_emissions,
        bus1="co2 atmosphere",
        carrier="process emissions",
        p_nom_extendable=True,
        efficiency=1.0,
    )

    # assume enough local waste heat for CC
    n.add(
        "Link",
        spatial.co2.locations,
        suffix=" process emissions CC",
        bus0=spatial.co2.process_emissions,
        bus1="co2 atmosphere",
        bus2=spatial.co2.nodes,
        carrier="process emissions CC",
        p_nom_extendable=True,
        capital_cost=costs.at["cement capture", "fixed"],
        efficiency=1 - costs.at["cement capture", "capture_rate"],
        efficiency2=costs.at["cement capture", "capture_rate"],
        lifetime=costs.at["cement capture", "lifetime"],
    )

    if options["ammonia"]:
        if options["ammonia"] == "regional":
            p_set = (
                industrial_demand.loc[spatial.ammonia.locations, "ammonia"].rename(
                    index=lambda x: x + " NH3"
                )
                / nhours
            )
        else:
            p_set = industrial_demand["ammonia"].sum() / nhours

        n.add(
            "Load",
            spatial.ammonia.nodes,
            bus=spatial.ammonia.nodes,
            carrier="NH3",
            p_set=p_set,
        )

    if industrial_demand[["coke", "coal"]].sum().sum() > 0:
        add_carrier_buses(n, "coal")

        mwh_coal_per_mwh_coke = 1.366  # from eurostat energy balance
        p_set = (
            industrial_demand["coal"]
            + mwh_coal_per_mwh_coke * industrial_demand["coke"]
        ) / nhours

        p_set.rename(lambda x: x + " coal for industry", inplace=True)

        if not options["regional_coal_demand"]:
            p_set = p_set.sum()

        n.add(
            "Bus",
            spatial.coal.industry,
            location=spatial.coal.demand_locations,
            carrier="coal for industry",
            unit="MWh_LHV",
        )

        n.add(
            "Load",
            spatial.coal.industry,
            bus=spatial.coal.industry,
            carrier="coal for industry",
            p_set=p_set,
        )

        n.add(
            "Link",
            spatial.coal.industry,
            bus0=spatial.coal.nodes,
            bus1=spatial.coal.industry,
            bus2="co2 atmosphere",
            carrier="coal for industry",
            p_nom_extendable=True,
            efficiency2=costs.at["coal", "CO2 intensity"],
        )


def add_waste_heat(n):
    # TODO options?

    logger.info("Add possibility to use industrial waste heat in district heating")

    # AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if not urban_central.empty:
        urban_central = urban_central.str[: -len(" urban central heat")]

        link_carriers = n.links.carrier.unique()

        # TODO what is the 0.95 and should it be a config option?
        if (
            options["use_fischer_tropsch_waste_heat"]
            and "Fischer-Tropsch" in link_carriers
        ):
            n.links.loc[urban_central + " Fischer-Tropsch", "bus3"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " Fischer-Tropsch", "efficiency3"] = (
                0.95 - n.links.loc[urban_central + " Fischer-Tropsch", "efficiency"]
            ) * options["use_fischer_tropsch_waste_heat"]

        if options["use_methanation_waste_heat"] and "Sabatier" in link_carriers:
            n.links.loc[urban_central + " Sabatier", "bus3"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " Sabatier", "efficiency3"] = (
                0.95 - n.links.loc[urban_central + " Sabatier", "efficiency"]
            ) * options["use_methanation_waste_heat"]

        # DEA quotes 15% of total input (11% of which are high-value heat)
        if options["use_haber_bosch_waste_heat"] and "Haber-Bosch" in link_carriers:
            n.links.loc[urban_central + " Haber-Bosch", "bus3"] = (
                urban_central + " urban central heat"
            )
            total_energy_input = (
                cf_industry["MWh_H2_per_tNH3_electrolysis"]
                + cf_industry["MWh_elec_per_tNH3_electrolysis"]
            ) / cf_industry["MWh_NH3_per_tNH3"]
            electricity_input = (
                cf_industry["MWh_elec_per_tNH3_electrolysis"]
                / cf_industry["MWh_NH3_per_tNH3"]
            )
            n.links.loc[urban_central + " Haber-Bosch", "efficiency3"] = (
                0.15 * total_energy_input / electricity_input
            ) * options["use_haber_bosch_waste_heat"]

        if (
            options["use_methanolisation_waste_heat"]
            and "methanolisation" in link_carriers
        ):
            n.links.loc[urban_central + " methanolisation", "bus4"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " methanolisation", "efficiency4"] = (
                costs.at["methanolisation", "heat-output"]
                / costs.at["methanolisation", "hydrogen-input"]
            ) * options["use_methanolisation_waste_heat"]

        # TODO integrate usable waste heat efficiency into technology-data from DEA
        if (
            options["use_electrolysis_waste_heat"]
            and "H2 Electrolysis" in link_carriers
        ):
            n.links.loc[urban_central + " H2 Electrolysis", "bus2"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " H2 Electrolysis", "efficiency2"] = (
                0.84 - n.links.loc[urban_central + " H2 Electrolysis", "efficiency"]
            ) * options["use_electrolysis_waste_heat"]

        if options["use_fuel_cell_waste_heat"] and "H2 Fuel Cell" in link_carriers:
            n.links.loc[urban_central + " H2 Fuel Cell", "bus2"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " H2 Fuel Cell", "efficiency2"] = (
                0.95 - n.links.loc[urban_central + " H2 Fuel Cell", "efficiency"]
            ) * options["use_fuel_cell_waste_heat"]


def add_agriculture(n, costs):
    logger.info("Add agriculture, forestry and fishing sector.")

    nodes = pop_layout.index
    nhours = n.snapshot_weightings.generators.sum()

    # electricity

    n.add(
        "Load",
        nodes,
        suffix=" agriculture electricity",
        bus=nodes,
        carrier="agriculture electricity",
        p_set=pop_weighted_energy_totals.loc[nodes, "total agriculture electricity"]
        * 1e6
        / nhours,
    )

    # heat

    n.add(
        "Load",
        nodes,
        suffix=" agriculture heat",
        bus=nodes + " services rural heat",
        carrier="agriculture heat",
        p_set=pop_weighted_energy_totals.loc[nodes, "total agriculture heat"]
        * 1e6
        / nhours,
    )

    # machinery

    electric_share = get(
        options["agriculture_machinery_electric_share"], investment_year
    )
    oil_share = get(options["agriculture_machinery_oil_share"], investment_year)

    total_share = electric_share + oil_share
    if total_share != 1:
        logger.warning(
            f"Total agriculture machinery shares sum up to {total_share:.2%}, corresponding to increased or decreased demand assumptions."
        )

    machinery_nodal_energy = (
        pop_weighted_energy_totals.loc[nodes, "total agriculture machinery"] * 1e6
    )

    if electric_share > 0:
        efficiency_gain = (
            options["agriculture_machinery_fuel_efficiency"]
            / options["agriculture_machinery_electric_efficiency"]
        )

        n.add(
            "Load",
            nodes,
            suffix=" agriculture machinery electric",
            bus=nodes,
            carrier="agriculture machinery electric",
            p_set=electric_share / efficiency_gain * machinery_nodal_energy / nhours,
        )

    if oil_share > 0:
        p_set = (
            oil_share
            * machinery_nodal_energy.rename(lambda x: x + " agriculture machinery oil")
            / nhours
        )

        if not options["regional_oil_demand"]:
            p_set = p_set.sum()

        n.add(
            "Bus",
            spatial.oil.agriculture_machinery,
            location=spatial.oil.demand_locations,
            carrier="agriculture machinery oil",
            unit="MWh_LHV",
        )

        n.add(
            "Load",
            spatial.oil.agriculture_machinery,
            bus=spatial.oil.agriculture_machinery,
            carrier="agriculture machinery oil",
            p_set=p_set,
        )

        n.add(
            "Link",
            spatial.oil.agriculture_machinery,
            bus0=spatial.oil.nodes,
            bus1=spatial.oil.agriculture_machinery,
            bus2="co2 atmosphere",
            carrier="agriculture machinery oil",
            p_nom_extendable=True,
            efficiency2=costs.at["oil", "CO2 intensity"],
        )


def decentral(n):
    """
    Removes the electricity transmission system.
    """
    n.lines.drop(n.lines.index, inplace=True)
    n.links.drop(n.links.index[n.links.carrier.isin(["DC", "B2B"])], inplace=True)


def remove_h2_network(n):
    n.links.drop(
        n.links.index[n.links.carrier.str.contains("H2 pipeline")], inplace=True
    )

    if "EU H2 Store" in n.stores.index:
        n.stores.drop("EU H2 Store", inplace=True)


def limit_individual_line_extension(n, maxext):
    logger.info(f"Limiting new HVAC and HVDC extensions to {maxext} MW")
    n.lines["s_nom_max"] = n.lines["s_nom"] + maxext
    hvdc = n.links.index[n.links.carrier == "DC"]
    n.links.loc[hvdc, "p_nom_max"] = n.links.loc[hvdc, "p_nom"] + maxext


aggregate_dict = {
    "p_nom": pd.Series.sum,
    "s_nom": pd.Series.sum,
    "v_nom": "max",
    "v_mag_pu_max": "min",
    "v_mag_pu_min": "max",
    "p_nom_max": pd.Series.sum,
    "s_nom_max": pd.Series.sum,
    "p_nom_min": pd.Series.sum,
    "s_nom_min": pd.Series.sum,
    "v_ang_min": "max",
    "v_ang_max": "min",
    "terrain_factor": "mean",
    "num_parallel": "sum",
    "p_set": "sum",
    "e_initial": "sum",
    "e_nom": pd.Series.sum,
    "e_nom_max": pd.Series.sum,
    "e_nom_min": pd.Series.sum,
    "state_of_charge_initial": "sum",
    "state_of_charge_set": "sum",
    "inflow": "sum",
    "p_max_pu": "first",
    "x": "mean",
    "y": "mean",
}


def cluster_heat_buses(n):
    """
    Cluster residential and service heat buses to one representative bus.

    This can be done to save memory and speed up optimisation
    """

    def define_clustering(attributes, aggregate_dict):
        """Define how attributes should be clustered.
        Input:
            attributes    : pd.Index()
            aggregate_dict: dictionary (key: name of attribute, value
                                        clustering method)

        Returns:
            agg           : clustering dictionary
        """
        keys = attributes.intersection(aggregate_dict.keys())
        agg = dict(
            zip(
                attributes.difference(keys),
                ["first"] * len(df.columns.difference(keys)),
            )
        )
        for key in keys:
            agg[key] = aggregate_dict[key]
        return agg

    logger.info("Cluster residential and service heat buses.")
    components = ["Bus", "Carrier", "Generator", "Link", "Load", "Store"]

    for c in n.iterate_components(components):
        df = c.df
        cols = df.columns[df.columns.str.contains("bus") | (df.columns == "carrier")]

        # rename columns and index
        df[cols] = df[cols].apply(
            lambda x: x.str.replace("residential ", "").str.replace("services ", ""),
            axis=1,
        )
        df = df.rename(
            index=lambda x: x.replace("residential ", "").replace("services ", "")
        )

        # cluster heat nodes
        # static dataframe
        agg = define_clustering(df.columns, aggregate_dict)
        df = df.groupby(level=0).agg(agg, numeric_only=False)
        # time-varying data
        pnl = c.pnl
        agg = define_clustering(pd.Index(pnl.keys()), aggregate_dict)
        for k in pnl.keys():

            def renamer(s):
                return s.replace("residential ", "").replace("services ", "")

            pnl[k] = pnl[k].T.groupby(renamer).agg(agg[k], numeric_only=False).T

        # remove unclustered assets of service/residential
        to_drop = c.df.index.difference(df.index)
        n.remove(c.name, to_drop)
        # add clustered assets
        to_add = df.index.difference(c.df.index)
        n.add(c.name, df.loc[to_add].index, **df.loc[to_add])


def set_temporal_aggregation(n, resolution, snapshot_weightings):
    """
    Aggregate time-varying data to the given snapshots.
    """
    if not resolution:
        logger.info("No temporal aggregation. Using native resolution.")
        return n
    elif "sn" in resolution.lower():
        # Representative snapshots are dealt with directly
        sn = int(resolution[:-2])
        logger.info("Use every %s snapshot as representative", sn)
        n.set_snapshots(n.snapshots[::sn])
        n.snapshot_weightings *= sn
        return n
    else:
        # Otherwise, use the provided snapshots
        snapshot_weightings = pd.read_csv(
            snapshot_weightings, index_col=0, parse_dates=True
        )

        # Define a series used for aggregation, mapping each hour in
        # n.snapshots to the closest previous timestep in
        # snapshot_weightings.index
        aggregation_map = (
            pd.Series(
                snapshot_weightings.index.get_indexer(n.snapshots), index=n.snapshots
            )
            .replace(-1, np.nan)
            .ffill()
            .astype(int)
            .map(lambda i: snapshot_weightings.index[i])
        )

        m = n.copy(with_time=False)
        m.set_snapshots(snapshot_weightings.index)
        m.snapshot_weightings = snapshot_weightings

        # Aggregation all time-varying data.
        for c in n.iterate_components():
            pnl = getattr(m, c.list_name + "_t")
            for k, df in c.pnl.items():
                if not df.empty:
                    if c.list_name == "stores" and k == "e_max_pu":
                        pnl[k] = df.groupby(aggregation_map).min()
                    elif c.list_name == "stores" and k == "e_min_pu":
                        pnl[k] = df.groupby(aggregation_map).max()
                    else:
                        pnl[k] = df.groupby(aggregation_map).mean()

        return m


def lossy_bidirectional_links(n, carrier, efficiencies={}):
    "Split bidirectional links into two unidirectional links to include transmission losses."

    carrier_i = n.links.query("carrier == @carrier").index

    if (
        not any((v != 1.0) or (v >= 0) for v in efficiencies.values())
        or carrier_i.empty
    ):
        return

    efficiency_static = efficiencies.get("efficiency_static", 1)
    efficiency_per_1000km = efficiencies.get("efficiency_per_1000km", 1)
    compression_per_1000km = efficiencies.get("compression_per_1000km", 0)

    logger.info(
        f"Specified losses for {carrier} transmission "
        f"(static: {efficiency_static}, per 1000km: {efficiency_per_1000km}, compression per 1000km: {compression_per_1000km}). "
        "Splitting bidirectional links."
    )

    n.links.loc[carrier_i, "p_min_pu"] = 0
    n.links.loc[carrier_i, "efficiency"] = (
        efficiency_static
        * efficiency_per_1000km ** (n.links.loc[carrier_i, "length"] / 1e3)
    )
    rev_links = (
        n.links.loc[carrier_i].copy().rename({"bus0": "bus1", "bus1": "bus0"}, axis=1)
    )
    rev_links["length_original"] = rev_links["length"]
    rev_links["capital_cost"] = 0
    rev_links["length"] = 0
    rev_links["reversed"] = True
    rev_links.index = rev_links.index.map(lambda x: x + "-reversed")

    n.links = pd.concat([n.links, rev_links], sort=False)
    n.links["reversed"] = n.links["reversed"].fillna(False).infer_objects(copy=False)
    n.links["length_original"] = n.links["length_original"].fillna(n.links.length)

    # do compression losses after concatenation to take electricity consumption at bus0 in either direction
    carrier_i = n.links.query("carrier == @carrier").index
    if compression_per_1000km > 0:
        n.links.loc[carrier_i, "bus2"] = n.links.loc[carrier_i, "bus0"].map(
            n.buses.location
        )  # electricity
        n.links.loc[carrier_i, "efficiency2"] = (
            -compression_per_1000km * n.links.loc[carrier_i, "length_original"] / 1e3
        )


def add_enhanced_geothermal(n, egs_potentials, egs_overlap, costs):
    """
    Adds EGS potential to model.

    Built in scripts/build_egs_potentials.py
    """

    if len(spatial.geothermal_heat.nodes) > 1:
        logger.warning(
            "'add_enhanced_geothermal' not implemented for multiple geothermal nodes."
        )
    logger.info(
        "[EGS] implemented with 2020 CAPEX from Aghahosseini et al 2021: 'From hot rock to...'."
    )
    logger.info(
        "[EGS] Recommended usage scales CAPEX to future cost expectations using config 'adjustments'."
    )
    logger.info("[EGS] During this the relevant carriers are:")
    logger.info("[EGS] drilling part -> 'geothermal heat'")
    logger.info(
        "[EGS] electricity generation part -> 'geothermal organic rankine cycle'"
    )
    logger.info("[EGS] district heat distribution part -> 'geothermal district heat'")

    egs_config = snakemake.params["sector"]["enhanced_geothermal"]
    costs_config = snakemake.config["costs"]

    # matrix defining the overlap between gridded geothermal potential estimation, and bus regions
    overlap = pd.read_csv(egs_overlap, index_col=0)
    overlap.columns = overlap.columns.astype(int)
    egs_potentials = pd.read_csv(egs_potentials, index_col=0)

    Nyears = n.snapshot_weightings.generators.sum() / 8760
    dr = costs_config["fill_values"]["discount rate"]
    lt = costs.at["geothermal", "lifetime"]
    FOM = costs.at["geothermal", "FOM"]

    egs_annuity = calculate_annuity(lt, dr)

    # under egs optimism, the expected cost reductions also cover costs for ORC
    # hence, the ORC costs are no longer taken from technology-data
    orc_capex = costs.at["organic rankine cycle", "investment"]

    # cost for ORC is subtracted, as it is already included in the geothermal cost.
    # The orc cost are attributed to a separate link representing the ORC.
    # also capital_cost conversion Euro/kW -> Euro/MW

    egs_potentials["capital_cost"] = (
        (egs_annuity + FOM / (1.0 + FOM))
        * (egs_potentials["CAPEX"] * 1e3 - orc_capex)
        * Nyears
    )

    assert (
        egs_potentials["capital_cost"] > 0
    ).all(), "Error in EGS cost, negative values found."

    orc_annuity = calculate_annuity(costs.at["organic rankine cycle", "lifetime"], dr)
    orc_capital_cost = (orc_annuity + FOM / (1 + FOM)) * orc_capex * Nyears

    efficiency_orc = costs.at["organic rankine cycle", "efficiency"]
    efficiency_dh = costs.at["geothermal", "district heat-input"]

    # p_nom_max conversion GW -> MW
    egs_potentials["p_nom_max"] = egs_potentials["p_nom_max"] * 1000.0

    # not using add_carrier_buses, as we are not interested in a Store
    n.add("Carrier", "geothermal heat")

    n.add(
        "Bus",
        spatial.geothermal_heat.nodes,
        carrier="geothermal heat",
        unit="MWh_th",
    )

    n.add(
        "Generator",
        spatial.geothermal_heat.nodes,
        bus=spatial.geothermal_heat.nodes,
        carrier="geothermal heat",
        p_nom_extendable=True,
    )

    if egs_config["var_cf"]:
        efficiency = pd.read_csv(
            snakemake.input.egs_capacity_factors, parse_dates=True, index_col=0
        )
        if snakemake.config["clustering"]["temporal"]["resolution_sector"]:
            efficiency = efficiency.resample(
                snakemake.config["clustering"]["temporal"]["resolution_sector"]
            ).mean()
        logger.info("Adding Enhanced Geothermal with time-varying capacity factors.")
    else:
        efficiency = 1.0

    # if urban central heat exists, adds geothermal as CHP
    as_chp = "urban central heat" in n.loads.carrier.unique()

    if as_chp:
        logger.info("Adding EGS as Combined Heat and Power.")

    else:
        logger.info("Adding EGS for Electricity Only.")

    for bus, bus_overlap in overlap.iterrows():
        if not bus_overlap.sum():
            continue

        overlap = bus_overlap.loc[bus_overlap > 0.0]
        bus_egs = egs_potentials.loc[overlap.index]

        if not len(bus_egs):
            continue

        bus_egs["p_nom_max"] = bus_egs["p_nom_max"].multiply(bus_overlap)
        bus_egs = bus_egs.loc[bus_egs.p_nom_max > 0.0]

        appendix = " " + pd.Index(np.arange(len(bus_egs)).astype(str))

        # add surface bus
        n.add(
            "Bus",
            pd.Index([f"{bus} geothermal heat surface"]),
            location=bus,
            unit="MWh_th",
            carrier="geothermal heat",
        )

        bus_egs.index = np.arange(len(bus_egs)).astype(str)
        well_name = f"{bus} enhanced geothermal" + appendix

        if egs_config["var_cf"]:
            bus_eta = pd.concat(
                (efficiency[bus].rename(idx) for idx in well_name),
                axis=1,
            )
        else:
            bus_eta = efficiency

        p_nom_max = bus_egs["p_nom_max"]
        capital_cost = bus_egs["capital_cost"]
        bus1 = pd.Series(f"{bus} geothermal heat surface", well_name)

        # adding geothermal wells as multiple generators to represent supply curve
        n.add(
            "Link",
            well_name,
            bus0=spatial.geothermal_heat.nodes,
            bus1=bus1,
            carrier="geothermal heat",
            p_nom_extendable=True,
            p_nom_max=p_nom_max.set_axis(well_name) / efficiency_orc,
            capital_cost=capital_cost.set_axis(well_name) * efficiency_orc,
            efficiency=bus_eta,
            lifetime=costs.at["geothermal", "lifetime"],
        )

        # adding Organic Rankine Cycle as a single link
        n.add(
            "Link",
            bus + " geothermal organic rankine cycle",
            bus0=f"{bus} geothermal heat surface",
            bus1=bus,
            p_nom_extendable=True,
            carrier="geothermal organic rankine cycle",
            capital_cost=orc_capital_cost * efficiency_orc,
            efficiency=efficiency_orc,
            lifetime=costs.at["organic rankine cycle", "lifetime"],
        )

        if as_chp and bus + " urban central heat" in n.buses.index:
            n.add(
                "Link",
                bus + " geothermal heat district heat",
                bus0=f"{bus} geothermal heat surface",
                bus1=bus + " urban central heat",
                carrier="geothermal district heat",
                capital_cost=orc_capital_cost
                * efficiency_orc
                * costs.at["geothermal", "district heat surcharge"]
                / 100.0,
                efficiency=efficiency_dh,
                p_nom_extendable=True,
                lifetime=costs.at["geothermal", "lifetime"],
            )

        if egs_config["flexible"]:
            # this StorageUnit represents flexible operation using the geothermal reservoir.
            # Hence, it is counter-intuitive to install it at the surface bus,
            # this is however the more lean and computationally efficient solution.

            max_hours = egs_config["max_hours"]
            boost = egs_config["max_boost"]

            n.add(
                "StorageUnit",
                bus + " geothermal reservoir",
                bus=f"{bus} geothermal heat surface",
                carrier="geothermal heat",
                p_nom_extendable=True,
                p_min_pu=-boost,
                max_hours=max_hours,
                cyclic_state_of_charge=True,
            )


# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_sector_network",
            opts="",
            clusters="38",
            ll="vopt",
            sector_opts="",
            planning_horizons="2030",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    options = snakemake.params.sector
    cf_industry = snakemake.params.industry

    investment_year = int(snakemake.wildcards.planning_horizons)

    n = pypsa.Network(snakemake.input.network)

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    nhours = n.snapshot_weightings.generators.sum()
    nyears = nhours / 8760

    costs = prepare_costs(
        snakemake.input.costs,
        snakemake.params.costs,
        nyears,
    )

    pop_weighted_energy_totals = (
        pd.read_csv(snakemake.input.pop_weighted_energy_totals, index_col=0) * nyears
    )
    pop_weighted_heat_totals = (
        pd.read_csv(snakemake.input.pop_weighted_heat_totals, index_col=0) * nyears
    )
    pop_weighted_energy_totals.update(pop_weighted_heat_totals)

    landfall_lengths = {
        tech: settings["landfall_length"]
        for tech, settings in snakemake.params.renewable.items()
        if "landfall_length" in settings.keys()
    }
    patch_electricity_network(n, costs, landfall_lengths)

    fn = snakemake.input.heating_efficiencies
    year = int(snakemake.params["energy_totals_year"])
    heating_efficiencies = pd.read_csv(fn, index_col=[1, 0]).loc[year]

    spatial = define_spatial(pop_layout.index, options)

    if snakemake.params.foresight in ["myopic", "perfect"]:
        add_lifetime_wind_solar(n, costs)

        conventional = snakemake.params.conventional_carriers
        for carrier in conventional:
            add_carrier_buses(n, carrier)

    add_eu_bus(n)

    add_co2_tracking(n, costs, options)

    add_generation(n, costs)

    add_storage_and_grids(n, costs)

    if options["transport"]:
        add_land_transport(n, costs)

    if options["heating"]:
        add_heat(
            n=n,
            costs=costs,
            cop=xr.open_dataarray(snakemake.input.cop_profiles),
            direct_heat_source_utilisation_profile=xr.open_dataarray(
                snakemake.input.direct_heat_source_utilisation_profiles
            ),
        )

    if options["biomass"]:
        add_biomass(n, costs)

    if options["ammonia"]:
        add_ammonia(n, costs)

    if options["methanol"]:
        add_methanol(n, costs)

    if options["industry"]:
        add_industry(n, costs)

    if options["heating"]:
        add_waste_heat(n)

    if options["agriculture"]:  # requires H and I
        add_agriculture(n, costs)

    if options["dac"]:
        add_dac(n, costs)

    if not options["electricity_transmission_grid"]:
        decentral(n)

    if not options["H2_network"]:
        remove_h2_network(n)

    if options["co2_network"]:
        add_co2_network(n, costs)

    if options["allam_cycle_gas"]:
        add_allam_gas(n, costs)

    n = set_temporal_aggregation(
        n, snakemake.params.time_resolution, snakemake.input.snapshot_weightings
    )

    co2_budget = snakemake.params.co2_budget
    if isinstance(co2_budget, str) and co2_budget.startswith("cb"):
        fn = "results/" + snakemake.params.RDIR + "/csvs/carbon_budget_distribution.csv"
        if not os.path.exists(fn):
            emissions_scope = snakemake.params.emissions_scope
            input_co2 = snakemake.input.co2
            build_carbon_budget(
                co2_budget,
                snakemake.input.eurostat,
                fn,
                emissions_scope,
                input_co2,
                options,
            )
        co2_cap = pd.read_csv(fn, index_col=0).squeeze()
        limit = co2_cap.loc[investment_year]
    else:
        limit = get(co2_budget, investment_year)
    add_co2limit(n, options, nyears, limit)

    maxext = snakemake.params["lines"]["max_extension"]
    if maxext is not None:
        limit_individual_line_extension(n, maxext)

    if options["electricity_distribution_grid"]:
        insert_electricity_distribution_grid(n, costs)

    if options["enhanced_geothermal"].get("enable", False):
        logger.info("Adding Enhanced Geothermal Systems (EGS).")
        add_enhanced_geothermal(
            n, snakemake.input["egs_potentials"], snakemake.input["egs_overlap"], costs
        )

    if options["gas_distribution_grid"]:
        insert_gas_distribution_costs(n, costs)

    if options["electricity_grid_connection"]:
        add_electricity_grid_connection(n, costs)

    for k, v in options["transmission_efficiency"].items():
        lossy_bidirectional_links(n, k, v)

    # Workaround: Remove lines with conflicting (and unrealistic) properties
    # cf. https://github.com/PyPSA/pypsa-eur/issues/444
    if snakemake.config["solving"]["options"]["transmission_losses"]:
        idx = n.lines.query("num_parallel == 0").index
        logger.info(
            f"Removing {len(idx)} line(s) with properties conflicting with transmission losses functionality."
        )
        n.remove("Line", idx)

    first_year_myopic = (snakemake.params.foresight in ["myopic", "perfect"]) and (
        snakemake.params.planning_horizons[0] == investment_year
    )

    if options["cluster_heat_buses"] and not first_year_myopic:
        cluster_heat_buses(n)

    maybe_adjust_costs_and_potentials(
        n, snakemake.params["adjustments"], investment_year
    )

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))

    sanitize_carriers(n, snakemake.config)
    sanitize_locations(n)

    n.export_to_netcdf(snakemake.output[0])
