# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Adds all sector-coupling components to the network, including demand and supply
technologies for the buildings, transport and industry sectors.
"""

import logging
import os
import re
import uuid
from itertools import product

import country_converter as coco
import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import generate_periodic_profiles, update_config_with_sector_opts
from add_electricity import calculate_annuity, sanitize_carriers
from build_energy_totals import build_co2_totals, build_eea_co2, build_eurostat_co2
from geopy.extra.rate_limiter import RateLimiter
from geopy.geocoders import Nominatim
from networkx.algorithms import complement
from networkx.algorithms.connectivity.edge_augmentation import k_edge_augmentation
from pypsa.geo import haversine_pts
from pypsa.io import import_components_from_dataframe
from scipy.stats import beta
from shapely.geometry import Point

cc = coco.CountryConverter()

geolocator = Nominatim(user_agent=str(uuid.uuid4()), timeout=10)
geocode = RateLimiter(geolocator.geocode, min_delay_seconds=2)

logger = logging.getLogger(__name__)

from types import SimpleNamespace

spatial = SimpleNamespace()

from packaging.version import Version, parse

pd_version = parse(pd.__version__)
agg_group_kwargs = dict(numeric_only=False) if pd_version >= Version("1.3") else {}

DISTANCE_CRS = 3857


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

    if options.get("biomass_spatial", options["biomass_transport"]):
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
    else:
        spatial.gas.nodes = ["EU gas"]
        spatial.gas.locations = ["EU"]
        spatial.gas.biogas = ["EU biogas"]
        spatial.gas.industry = ["gas for industry"]
        spatial.gas.biogas_to_gas = ["EU biogas to gas"]
        if options.get("co2_spatial", options["co2network"]):
            spatial.gas.industry_cc = nodes + " gas for industry CC"
        else:
            spatial.gas.industry_cc = ["gas for industry CC"]

    spatial.gas.df = pd.DataFrame(vars(spatial.gas), index=nodes)

    # ammonia

    if options.get("ammonia"):
        spatial.ammonia = SimpleNamespace()
        if options.get("ammonia") == "regional":
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
    spatial.methanol = SimpleNamespace()
    spatial.methanol.nodes = ["EU methanol"]
    spatial.methanol.locations = ["EU"]

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
        sectors += ["rail non-elec", "road non-elec"]
    if "H" in opts:
        sectors += ["residential non-elec", "services non-elec"]
    if "I" in opts:
        sectors += [
            "industrial non-elec",
            "industrial processes",
            "domestic aviation",
            "international aviation",
            "domestic navigation",
            "international navigation",
        ]
    if "A" in opts:
        sectors += ["agriculture"]

    return sectors


def get(item, investment_year=None):
    """
    Check whether item depends on investment year.
    """
    if isinstance(item, dict):
        return item[investment_year]
    else:
        return item


def co2_emissions_year(
    countries, input_eurostat, opts, emissions_scope, report_year, input_co2, year
):
    """
    Calculate CO2 emissions in one specific year (e.g. 1990 or 2018).
    """
    eea_co2 = build_eea_co2(input_co2, year, emissions_scope)

    # TODO: read Eurostat data from year > 2014
    # this only affects the estimation of CO2 emissions for BA, RS, AL, ME, MK
    if year > 2014:
        eurostat_co2 = build_eurostat_co2(
            input_eurostat, countries, report_year, year=2014
        )
    else:
        eurostat_co2 = build_eurostat_co2(input_eurostat, countries, report_year, year)

    co2_totals = build_co2_totals(countries, eea_co2, eurostat_co2)

    sectors = emission_sectors_from_opts(opts)

    co2_emissions = co2_totals.loc[countries, sectors].sum().sum()

    # convert MtCO2 to GtCO2
    co2_emissions *= 0.001

    return co2_emissions


# TODO: move to own rule with sector-opts wildcard?
def build_carbon_budget(o, input_eurostat, fn, emissions_scope, report_year, input_co2):
    """
    Distribute carbon budget following beta or exponential transition path.
    """
    # opts?

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
        opts,
        emissions_scope,
        report_year,
        input_co2,
        year=1990,
    )

    # emissions at the beginning of the path (last year available 2018)
    e_0 = co2_emissions_year(
        countries,
        input_eurostat,
        opts,
        emissions_scope,
        report_year,
        input_co2,
        year=2018,
    )

    planning_horizons = snakemake.params.planning_horizons
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


# TODO merge issue with PyPSA-Eur
def update_wind_solar_costs(n, costs):
    """
    Update costs for wind and solar generators added with pypsa-eur to those
    cost in the planning year.
    """
    # NB: solar costs are also manipulated for rooftop
    # when distribution grid is inserted
    n.generators.loc[n.generators.carrier == "solar", "capital_cost"] = costs.at[
        "solar-utility", "fixed"
    ]

    n.generators.loc[n.generators.carrier == "onwind", "capital_cost"] = costs.at[
        "onwind", "fixed"
    ]

    # for offshore wind, need to calculated connection costs

    # assign clustered bus
    # map initial network -> simplified network
    busmap_s = pd.read_csv(snakemake.input.busmap_s, index_col=0).squeeze()
    busmap_s.index = busmap_s.index.astype(str)
    busmap_s = busmap_s.astype(str)
    # map simplified network -> clustered network
    busmap = pd.read_csv(snakemake.input.busmap, index_col=0).squeeze()
    busmap.index = busmap.index.astype(str)
    busmap = busmap.astype(str)
    # map initial network -> clustered network
    clustermaps = busmap_s.map(busmap)

    # code adapted from pypsa-eur/scripts/add_electricity.py
    for connection in ["dc", "ac"]:
        tech = "offwind-" + connection
        profile = snakemake.input["profile_offwind_" + connection]
        with xr.open_dataset(profile) as ds:
            underwater_fraction = ds["underwater_fraction"].to_pandas()
            connection_cost = (
                snakemake.params.length_factor
                * ds["average_distance"].to_pandas()
                * (
                    underwater_fraction
                    * costs.at[tech + "-connection-submarine", "fixed"]
                    + (1.0 - underwater_fraction)
                    * costs.at[tech + "-connection-underground", "fixed"]
                )
            )

            # convert to aggregated clusters with weighting
            weight = ds["weight"].to_pandas()

            # e.g. clusters == 37m means that VRE generators are left
            # at clustering of simplified network, but that they are
            # connected to 37-node network
            if snakemake.wildcards.clusters[-1:] == "m":
                genmap = busmap_s
            else:
                genmap = clustermaps

            connection_cost = (connection_cost * weight).groupby(
                genmap
            ).sum() / weight.groupby(genmap).sum()

            capital_cost = (
                costs.at["offwind", "fixed"]
                + costs.at[tech + "-station", "fixed"]
                + connection_cost
            )

            logger.info(
                "Added connection cost of {:0.0f}-{:0.0f} Eur/MW/a to {}".format(
                    connection_cost[0].min(), connection_cost[0].max(), tech
                )
            )

            n.generators.loc[
                n.generators.carrier == tech, "capital_cost"
            ] = capital_cost.rename(index=lambda node: node + " " + tech)


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
    capital_cost = costs.at["gas storage", "fixed"] if carrier == "gas" else 0.02

    n.madd("Bus", nodes, location=location, carrier=carrier, unit=unit)

    n.madd(
        "Store",
        nodes + " Store",
        bus=nodes,
        e_nom_extendable=True,
        e_cyclic=True,
        carrier=carrier,
        capital_cost=capital_cost,
    )

    if carrier in costs.index and costs.at[carrier, "fuel"] > 0:
        n.madd(
            "Generator",
            nodes,
            bus=nodes,
            p_nom_extendable=True,
            carrier=carrier,
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
        n.mremove(c.name, names)
        n.carriers.drop(to_remove, inplace=True, errors="ignore")


# TODO: PyPSA-Eur merge issue
def remove_non_electric_buses(n):
    """
    Remove buses from pypsa-eur with carriers which are not AC buses.
    """
    to_drop = list(n.buses.query("carrier not in ['AC', 'DC']").carrier.unique())
    if to_drop:
        logger.info(f"Drop buses from PyPSA-Eur with carrier: {to_drop}")
        n.buses = n.buses[n.buses.carrier.isin(["AC", "DC"])]


def patch_electricity_network(n):
    remove_elec_base_techs(n)
    remove_non_electric_buses(n)
    update_wind_solar_costs(n, costs)
    n.loads["carrier"] = "electricity"
    n.buses["location"] = n.buses.index
    n.buses["unit"] = "MWh_el"
    # remove trailing white space of load index until new PyPSA version after v0.18.
    n.loads.rename(lambda x: x.strip(), inplace=True)
    n.loads_t.p_set.rename(lambda x: x.strip(), axis=1, inplace=True)


def add_co2_tracking(n, options):
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

    # this tracks CO2 stored, e.g. underground
    n.madd(
        "Bus",
        spatial.co2.nodes,
        location=spatial.co2.locations,
        carrier="co2 stored",
        unit="t_co2",
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
        e_nom_max = e_nom_max.rename(index=lambda x: x + " co2 stored")
    else:
        e_nom_max = np.inf

    n.madd(
        "Store",
        spatial.co2.nodes,
        e_nom_extendable=True,
        e_nom_max=e_nom_max,
        capital_cost=options["co2_sequestration_cost"],
        carrier="co2 stored",
        bus=spatial.co2.nodes,
    )

    n.add("Carrier", "co2 stored")

    if options["co2_vent"]:
        n.madd(
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

    n.madd(
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


def add_allam_cycle_gas(n, costs):
    logger.info("Adding Allam cycle gas power plants.")

    nodes = pop_layout.index

    n.madd(
        "Link",
        nodes,
        suffix=" allam gas",
        bus0=spatial.gas.df.loc[nodes, "nodes"].values,
        bus1=nodes,
        bus2=spatial.co2.df.loc[nodes, "nodes"].values,
        bus3="co2 atmosphere",
        carrier="allam gas",
        p_nom_extendable=True,
        # TODO: add costs to technology-data
        capital_cost=0.66
        * 1.832e6
        * (
            calculate_annuity(25, 0.07) + 0.0385
        ),  # efficiency * EUR/MW * (annuity + FOM)
        marginal_cost=2,
        efficiency=0.66,
        efficiency2=0.98 * costs.at["gas", "CO2 intensity"],
        efficiency3=0.02 * costs.at["gas", "CO2 intensity"],
        lifetime=25,
    )


def add_methanol_to_power(n, costs, types={}):
    # TODO: add costs to technology-data

    nodes = pop_layout.index

    if types["allam"]:
        logger.info("Adding Allam cycle methanol power plants.")

        n.madd(
            "Link",
            nodes,
            suffix=" allam methanol",
            bus0=spatial.methanol.nodes,
            bus1=nodes,
            bus2=spatial.co2.df.loc[nodes, "nodes"].values,
            bus3="co2 atmosphere",
            carrier="allam methanol",
            p_nom_extendable=True,
            capital_cost=0.59
            * 1.832e6
            * calculate_annuity(25, 0.07),  # efficiency * EUR/MW * annuity
            marginal_cost=2,
            efficiency=0.59,
            efficiency2=0.98 * costs.at["methanolisation", "carbondioxide-input"],
            efficiency3=0.02 * costs.at["methanolisation", "carbondioxide-input"],
            lifetime=25,
        )

    if types["ccgt"]:
        logger.info("Adding methanol CCGT power plants.")

        # efficiency * EUR/MW * (annuity + FOM)
        capital_cost = 0.58 * 916e3 * (calculate_annuity(25, 0.07) + 0.035)

        n.madd(
            "Link",
            nodes,
            suffix=" CCGT methanol",
            bus0=spatial.methanol.nodes,
            bus1=nodes,
            bus2="co2 atmosphere",
            carrier="CCGT methanol",
            p_nom_extendable=True,
            capital_cost=capital_cost,
            marginal_cost=2,
            efficiency=0.58,
            efficiency2=costs.at["methanolisation", "carbondioxide-input"],
            lifetime=25,
        )

    if types["ccgt_cc"]:
        logger.info(
            "Adding methanol CCGT power plants with post-combustion carbon capture."
        )

        # TODO consider efficiency changes / energy inputs for CC

        capital_cost_cc = (
            capital_cost
            + costs.at["cement capture", "fixed"]
            * costs.at["methanolisation", "carbondioxide-input"]
        )

        n.madd(
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
            marginal_cost=2,
            efficiency=0.58,
            efficiency2=costs.at["cement capture", "capture_rate"]
            * costs.at["methanolisation", "carbondioxide-input"],
            efficiency3=(1 - costs.at["cement capture", "capture_rate"])
            * costs.at["methanolisation", "carbondioxide-input"],
            lifetime=25,
        )

    if types["ocgt"]:
        logger.info("Adding methanol OCGT power plants.")

        n.madd(
            "Link",
            nodes,
            suffix=" OCGT methanol",
            bus0=spatial.methanol.nodes,
            bus1=nodes,
            bus2="co2 atmosphere",
            carrier="OCGT methanol",
            p_nom_extendable=True,
            capital_cost=0.35
            * 458e3
            * (
                calculate_annuity(25, 0.07) + 0.035
            ),  # efficiency * EUR/MW * (annuity + FOM)
            marginal_cost=2,
            efficiency=0.35,
            efficiency2=costs.at["methanolisation", "carbondioxide-input"],
            lifetime=25,
        )


def add_methanol_to_kerosene(n, costs):
    logger.info("Adding methanol-to-kerosene.")

    tech = "methanol-to-kerosene"

    n.madd(
        "Link",
        spatial.nodes,
        suffix=f" {tech}",
        bus0=spatial.methanol.nodes,
        bus1=spatial.oil.nodes,
        bus2=spatial.h2.nodes,
        efficiency=costs.at[tech, "methanol-input"],
        efficiency2=-costs.at[tech, "hydrogen-input"]
        / costs.at[tech, "methanol-input"],
        p_nom_extendable=True,
        p_min_pu=1,
        carrier=tech,
    )


def add_methanol_reforming(n, costs):
    logger.info("Adding methanol steam reforming.")

    nodes = pop_layout.index

    tech = "Methanol steam reforming"

    capital_cost = costs.at[tech, "fixed"] / costs.at[tech, "methanol-input"]

    n.madd(
        "Link",
        nodes,
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

    nodes = pop_layout.index

    # TODO: heat release and electricity demand for process and carbon capture
    # but the energy demands for carbon capture have not yet been added for other CC processes
    # 10.1016/j.rser.2020.110171: 0.129 kWh_e/kWh_H2, -0.09 kWh_heat/kWh_H2

    tech = "Methanol steam reforming"

    capital_cost = costs.at[tech, "fixed"] / costs.at[tech, "methanol-input"]

    capital_cost_cc = (
        capital_cost
        + costs.at["cement capture", "fixed"]
        * costs.at["methanolisation", "carbondioxide-input"]
    )

    n.madd(
        "Link",
        nodes,
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

    efficiency2 = -(
        costs.at["direct air capture", "electricity-input"]
        + costs.at["direct air capture", "compression-electricity-input"]
    )
    efficiency3 = -(
        costs.at["direct air capture", "heat-input"]
        - costs.at["direct air capture", "compression-heat-output"]
    )

    n.madd(
        "Link",
        heat_buses.str.replace(" heat", " DAC"),
        bus0="co2 atmosphere",
        bus1=spatial.co2.df.loc[locations, "nodes"].values,
        bus2=locations.values,
        bus3=heat_buses,
        carrier="DAC",
        capital_cost=costs.at["direct air capture", "fixed"],
        efficiency=1.0,
        efficiency2=efficiency2,
        efficiency3=efficiency3,
        p_nom_extendable=True,
        lifetime=costs.at["direct air capture", "lifetime"],
    )


def add_co2limit(n, nyears=1.0, limit=0.0):
    logger.info(f"Adding CO2 budget limit as per unit of 1990 levels of {limit}")

    countries = snakemake.params.countries

    sectors = emission_sectors_from_opts(opts)

    # convert Mt to tCO2
    co2_totals = 1e6 * pd.read_csv(snakemake.input.co2_totals_name, index_col=0)

    co2_limit = co2_totals.loc[countries, sectors].sum().sum()

    co2_limit *= limit * nyears

    n.add(
        "GlobalConstraint",
        "CO2Limit",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=co2_limit,
    )


# TODO PyPSA-Eur merge issue
def average_every_nhours(n, offset):
    logger.info(f"Resampling the network to {offset}")
    m = n.copy(with_time=False)

    snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name + "_t")
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

    fallback = {"OCGT": "gas"}
    conventionals = options.get("conventional_generation", fallback)

    for generator, carrier in conventionals.items():
        carrier_nodes = vars(spatial)[carrier].nodes

        add_carrier_buses(n, carrier, carrier_nodes)

        n.madd(
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

    cf_industry = snakemake.params.industry

    n.add("Carrier", "NH3")

    n.madd(
        "Bus", spatial.ammonia.nodes, location=spatial.ammonia.locations, carrier="NH3"
    )

    n.madd(
        "Link",
        nodes,
        suffix=" Haber-Bosch",
        bus0=nodes,
        bus1=spatial.ammonia.nodes,
        bus2=nodes + " H2",
        p_nom_extendable=True,
        carrier="Haber-Bosch",
        efficiency=1
        / (
            cf_industry["MWh_elec_per_tNH3_electrolysis"]
            / cf_industry["MWh_NH3_per_tNH3"]
        ),  # output: MW_NH3 per MW_elec
        efficiency2=-cf_industry["MWh_H2_per_tNH3_electrolysis"]
        / cf_industry["MWh_elec_per_tNH3_electrolysis"],  # input: MW_H2 per MW_elec
        capital_cost=costs.at["Haber-Bosch", "fixed"],
        lifetime=costs.at["Haber-Bosch", "lifetime"],
    )

    n.madd(
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
    n.madd(
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


def add_wave(n, wave_cost_factor):
    # TODO: handle in Snakefile
    wave_fn = "data/WindWaveWEC_GLTB.xlsx"

    # in kW
    capacity = pd.Series({"Attenuator": 750, "F2HB": 1000, "MultiPA": 600})

    # in EUR/MW
    annuity_factor = calculate_annuity(25, 0.07) + 0.03
    costs = (
        1e6
        * wave_cost_factor
        * annuity_factor
        * pd.Series({"Attenuator": 2.5, "F2HB": 2, "MultiPA": 1.5})
    )

    sheets = pd.read_excel(
        wave_fn,
        sheet_name=["FirthForth", "Hebrides"],
        usecols=["Attenuator", "F2HB", "MultiPA"],
        index_col=0,
        skiprows=[0],
        parse_dates=True,
    )

    wave = pd.concat(
        [sheets[l].divide(capacity, axis=1) for l in locations], keys=locations, axis=1
    )

    for wave_type in costs.index:
        n.add(
            "Generator",
            "Hebrides " + wave_type,
            bus="GB4 0",  # TODO this location is hardcoded
            p_nom_extendable=True,
            carrier="wave",
            capital_cost=costs[wave_type],
            p_max_pu=wave["Hebrides", wave_type],
        )


def insert_electricity_distribution_grid(n, costs):
    # TODO pop_layout?
    # TODO options?

    cost_factor = options["electricity_distribution_grid_cost_factor"]

    logger.info(
        f"Inserting electricity distribution grid with investment cost factor of {cost_factor:.2f}"
    )

    nodes = pop_layout.index

    n.madd(
        "Bus",
        nodes + " low voltage",
        location=nodes,
        carrier="low voltage",
        unit="MWh_el",
    )

    n.madd(
        "Link",
        nodes + " electricity distribution grid",
        bus0=nodes,
        bus1=nodes + " low voltage",
        p_nom_extendable=True,
        p_min_pu=-1,
        carrier="electricity distribution grid",
        efficiency=1,
        lifetime=costs.at["electricity distribution grid", "lifetime"],
        capital_cost=costs.at["electricity distribution grid", "fixed"] * cost_factor,
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
    n.generators.loc[solar, "capital_cost"] = costs.at["solar-utility", "fixed"]
    if snakemake.wildcards.clusters[-1:] == "m":
        simplified_pop_layout = pd.read_csv(
            snakemake.input.simplified_pop_layout, index_col=0
        )
        pop_solar = simplified_pop_layout.total.rename(index=lambda x: x + " solar")
    else:
        pop_solar = pop_layout.total.rename(index=lambda x: x + " solar")

    # add max solar rooftop potential assuming 0.1 kW/m2 and 10 m2/person,
    # i.e. 1 kW/person (population data is in thousands of people) so we get MW
    potential = 0.1 * 10 * pop_solar

    n.madd(
        "Generator",
        solar,
        suffix=" rooftop",
        bus=n.generators.loc[solar, "bus"] + " low voltage",
        carrier="solar rooftop",
        p_nom_extendable=True,
        p_nom_max=potential,
        marginal_cost=n.generators.loc[solar, "marginal_cost"],
        capital_cost=costs.at["solar-rooftop", "fixed"],
        efficiency=n.generators.loc[solar, "efficiency"],
        p_max_pu=n.generators_t.p_max_pu[solar],
        lifetime=costs.at["solar-rooftop", "lifetime"],
    )

    n.add("Carrier", "home battery")

    n.madd(
        "Bus",
        nodes + " home battery",
        location=nodes,
        carrier="home battery",
        unit="MWh_el",
    )

    n.madd(
        "Store",
        nodes + " home battery",
        bus=nodes + " home battery",
        e_cyclic=True,
        e_nom_extendable=True,
        carrier="home battery",
        capital_cost=costs.at["home battery storage", "fixed"],
        lifetime=costs.at["battery storage", "lifetime"],
    )

    n.madd(
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

    n.madd(
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

    capital_cost = costs.loc["electricity distribution grid"]["fixed"] * f_costs

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
    carriers = ["onwind", "solar"]

    gens = n.generators.index[n.generators.carrier.isin(carriers)]

    n.generators.loc[gens, "capital_cost"] += costs.at[
        "electricity grid connection", "fixed"
    ]


def add_storage_and_grids(n, costs):
    logger.info("Add hydrogen storage")

    nodes = pop_layout.index

    n.add("Carrier", "H2")

    n.madd("Bus", nodes + " H2", location=nodes, carrier="H2", unit="MWh_LHV")

    n.madd(
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

        n.madd(
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

        n.madd(
            "Link",
            nodes + " H2 turbine",
            bus0=nodes + " H2",
            bus1=nodes,
            p_nom_extendable=True,
            carrier="H2 turbine",
            efficiency=costs.at["OCGT", "efficiency"],
            capital_cost=costs.at["OCGT", "fixed"]
            * costs.at["OCGT", "efficiency"],  # NB: fixed cost is per MWel
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

        n.madd(
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

    n.madd(
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
        else:
            gas_pipes["p_nom_max"] = np.inf
            gas_pipes["p_nom_min"] = gas_pipes.p_nom
            gas_pipes["capital_cost"] = (
                gas_pipes.length * costs.at["CH4 (g) pipeline", "fixed"]
            )

        n.madd(
            "Link",
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
            lifetime=costs.at["CH4 (g) pipeline", "lifetime"],
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
        k_edge = options.get("gas_network_connectivity_upgrade", 3)
        augmentation = list(
            k_edge_augmentation(G, k_edge, avail=complement_edges.values)
        )

        if augmentation:
            new_gas_pipes = pd.DataFrame(augmentation, columns=["bus0", "bus1"])
            new_gas_pipes["length"] = new_gas_pipes.apply(haversine, axis=1)

            new_gas_pipes.index = new_gas_pipes.apply(
                lambda x: f"gas pipeline new {x.bus0} <-> {x.bus1}", axis=1
            )

            n.madd(
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

        n.madd(
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

    if options.get("H2_network", True):
        logger.info("Add options for new hydrogen pipelines.")

        h2_pipes = create_network_topology(
            n, "H2 pipeline ", carriers=["DC", "gas pipeline"]
        )

        # TODO Add efficiency losses
        n.madd(
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

    n.madd("Bus", nodes + " battery", location=nodes, carrier="battery", unit="MWh_el")

    n.madd(
        "Store",
        nodes + " battery",
        bus=nodes + " battery",
        e_cyclic=True,
        e_nom_extendable=True,
        carrier="battery",
        capital_cost=costs.at["battery storage", "fixed"],
        lifetime=costs.at["battery storage", "lifetime"],
    )

    n.madd(
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

    n.madd(
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
        n.madd(
            "Link",
            spatial.nodes,
            suffix=" Sabatier",
            bus0=nodes + " H2",
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            p_nom_extendable=True,
            carrier="Sabatier",
            efficiency=costs.at["methanation", "efficiency"],
            efficiency2=-costs.at["methanation", "efficiency"]
            * costs.at["gas", "CO2 intensity"],
            capital_cost=costs.at["methanation", "fixed"]
            * costs.at["methanation", "efficiency"],  # costs given per kW_gas
            lifetime=costs.at["methanation", "lifetime"],
        )

    if options["helmeth"]:
        n.madd(
            "Link",
            spatial.nodes,
            suffix=" helmeth",
            bus0=nodes,
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            carrier="helmeth",
            p_nom_extendable=True,
            efficiency=costs.at["helmeth", "efficiency"],
            efficiency2=-costs.at["helmeth", "efficiency"]
            * costs.at["gas", "CO2 intensity"],
            capital_cost=costs.at["helmeth", "fixed"],
            lifetime=costs.at["helmeth", "lifetime"],
        )

    if options.get("coal_cc"):
        n.madd(
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

    if options["SMR"]:
        n.madd(
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

        n.madd(
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


def add_land_transport(n, costs):
    # TODO options?

    logger.info("Add land transport")
    nhours = n.snapshot_weightings.generators.sum()

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

    fuel_cell_share = get(options["land_transport_fuel_cell_share"], investment_year)
    electric_share = get(options["land_transport_electric_share"], investment_year)
    ice_share = get(options["land_transport_ice_share"], investment_year)

    total_share = fuel_cell_share + electric_share + ice_share
    if total_share != 1:
        logger.warning(
            f"Total land transport shares sum up to {total_share:.2%}, corresponding to increased or decreased demand assumptions."
        )

    logger.info(f"FCEV share: {fuel_cell_share*100}%")
    logger.info(f"EV share: {electric_share*100}%")
    logger.info(f"ICEV share: {ice_share*100}%")

    nodes = pop_layout.index

    if electric_share > 0:
        n.add("Carrier", "Li ion")

        n.madd(
            "Bus",
            nodes,
            location=nodes,
            suffix=" EV battery",
            carrier="Li ion",
            unit="MWh_el",
        )

        p_set = (
            electric_share
            * (
                transport[nodes]
                + cycling_shift(transport[nodes], 1)
                + cycling_shift(transport[nodes], 2)
            )
            / 3
        )

        n.madd(
            "Load",
            nodes,
            suffix=" land transport EV",
            bus=nodes + " EV battery",
            carrier="land transport EV",
            p_set=p_set,
        )

        p_nom = number_cars * options.get("bev_charge_rate", 0.011) * electric_share

        n.madd(
            "Link",
            nodes,
            suffix=" BEV charger",
            bus0=nodes,
            bus1=nodes + " EV battery",
            p_nom=p_nom,
            carrier="BEV charger",
            p_max_pu=avail_profile[nodes],
            efficiency=options.get("bev_charge_efficiency", 0.9),
            # These were set non-zero to find LU infeasibility when availability = 0.25
            # p_nom_extendable=True,
            # p_nom_min=p_nom,
            # capital_cost=1e6,  #i.e. so high it only gets built where necessary
        )

    if electric_share > 0 and options["v2g"]:
        n.madd(
            "Link",
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
        e_nom = (
            number_cars
            * options.get("bev_energy", 0.05)
            * options["bev_availability"]
            * electric_share
        )

        n.madd(
            "Store",
            nodes,
            suffix=" battery storage",
            bus=nodes + " EV battery",
            carrier="battery storage",
            e_cyclic=True,
            e_nom=e_nom,
            e_max_pu=1,
            e_min_pu=dsm_profile[nodes],
        )

    if fuel_cell_share > 0:
        n.madd(
            "Load",
            nodes,
            suffix=" land transport fuel cell",
            bus=nodes + " H2",
            carrier="land transport fuel cell",
            p_set=fuel_cell_share
            / options["transport_fuel_cell_efficiency"]
            * transport[nodes],
        )

    if ice_share > 0:
        if "oil" not in n.buses.carrier.unique():
            n.madd(
                "Bus",
                spatial.oil.nodes,
                location=spatial.oil.locations,
                carrier="oil",
                unit="MWh_LHV",
            )

        ice_efficiency = options["transport_internal_combustion_efficiency"]

        n.madd(
            "Load",
            nodes,
            suffix=" land transport oil",
            bus=spatial.oil.nodes,
            carrier="land transport oil",
            p_set=ice_share / ice_efficiency * transport[nodes],
        )

        co2 = (
            ice_share
            / ice_efficiency
            * transport[nodes].sum().sum()
            / nhours
            * costs.at["oil", "CO2 intensity"]
        )

        n.add(
            "Load",
            "land transport oil emissions",
            bus="co2 atmosphere",
            carrier="land transport oil emissions",
            p_set=-co2,
        )


def build_heat_demand(n):
    # copy forward the daily average heat demand into each hour, so it can be multiplied by the intraday profile
    daily_space_heat_demand = (
        xr.open_dataarray(snakemake.input.heat_demand_total)
        .to_pandas()
        .reindex(index=n.snapshots, method="ffill")
    )

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
            weekly_profile=weekly_profile,
        )

        if use == "space":
            heat_demand_shape = daily_space_heat_demand * intraday_year_profile
        else:
            heat_demand_shape = intraday_year_profile

        heat_demand[f"{sector} {use}"] = (
            heat_demand_shape / heat_demand_shape.sum()
        ).multiply(pop_weighted_energy_totals[f"total {sector} {use}"]) * 1e6
        electric_heat_supply[f"{sector} {use}"] = (
            heat_demand_shape / heat_demand_shape.sum()
        ).multiply(pop_weighted_energy_totals[f"electricity {sector} {use}"]) * 1e6

    heat_demand = pd.concat(heat_demand, axis=1)
    electric_heat_supply = pd.concat(electric_heat_supply, axis=1)

    # subtract from electricity load since heat demand already in heat_demand
    electric_nodes = n.loads.index[n.loads.carrier == "electricity"]
    n.loads_t.p_set[electric_nodes] = (
        n.loads_t.p_set[electric_nodes]
        - electric_heat_supply.groupby(level=1, axis=1).sum()[electric_nodes]
    )

    return heat_demand


def add_heat(n, costs):
    logger.info("Add heat sector")

    sectors = ["residential", "services"]

    heat_demand = build_heat_demand(n)

    nodes, dist_fraction, urban_fraction = create_nodes_for_heat_sector()

    # NB: must add costs of central heating afterwards (EUR 400 / kWpeak, 50a, 1% FOM from Fraunhofer ISE)

    # exogenously reduce space heat demand
    if options["reduce_space_heat_exogenously"]:
        dE = get(options["reduce_space_heat_exogenously_factor"], investment_year)
        logger.info(f"Assumed space heat reduction of {dE:.2%}")
        for sector in sectors:
            heat_demand[sector + " space"] = (1 - dE) * heat_demand[sector + " space"]

    heat_systems = [
        "residential rural",
        "services rural",
        "residential urban decentral",
        "services urban decentral",
        "urban central",
    ]

    cop = {
        "air": xr.open_dataarray(snakemake.input.cop_air_total)
        .to_pandas()
        .reindex(index=n.snapshots),
        "ground": xr.open_dataarray(snakemake.input.cop_soil_total)
        .to_pandas()
        .reindex(index=n.snapshots),
    }

    if options["solar_thermal"]:
        solar_thermal = (
            xr.open_dataarray(snakemake.input.solar_thermal_total)
            .to_pandas()
            .reindex(index=n.snapshots)
        )
        # 1e3 converts from W/m^2 to MW/(1000m^2) = kW/m^2
        solar_thermal = options["solar_cf_correction"] * solar_thermal / 1e3

    for name in heat_systems:
        name_type = "central" if name == "urban central" else "decentral"

        n.add("Carrier", name + " heat")

        n.madd(
            "Bus",
            nodes[name] + f" {name} heat",
            location=nodes[name],
            carrier=name + " heat",
            unit="MWh_th",
        )

        ## Add heat load

        for sector in sectors:
            # heat demand weighting
            if "rural" in name:
                factor = 1 - urban_fraction[nodes[name]]
            elif "urban central" in name:
                factor = dist_fraction[nodes[name]]
            elif "urban decentral" in name:
                factor = urban_fraction[nodes[name]] - dist_fraction[nodes[name]]
            else:
                raise NotImplementedError(
                    f" {name} not in " f"heat systems: {heat_systems}"
                )

            if sector in name:
                heat_load = (
                    heat_demand[[sector + " water", sector + " space"]]
                    .groupby(level=1, axis=1)
                    .sum()[nodes[name]]
                    .multiply(factor)
                )

        if name == "urban central":
            heat_load = (
                heat_demand.groupby(level=1, axis=1)
                .sum()[nodes[name]]
                .multiply(
                    factor * (1 + options["district_heating"]["district_heating_loss"])
                )
            )

        n.madd(
            "Load",
            nodes[name],
            suffix=f" {name} heat",
            bus=nodes[name] + f" {name} heat",
            carrier=name + " heat",
            p_set=heat_load,
        )

        ## Add heat pumps

        heat_pump_type = "air" if "urban" in name else "ground"

        costs_name = f"{name_type} {heat_pump_type}-sourced heat pump"
        efficiency = (
            cop[heat_pump_type][nodes[name]]
            if options["time_dep_hp_cop"]
            else costs.at[costs_name, "efficiency"]
        )

        n.madd(
            "Link",
            nodes[name],
            suffix=f" {name} {heat_pump_type} heat pump",
            bus0=nodes[name],
            bus1=nodes[name] + f" {name} heat",
            carrier=f"{name} {heat_pump_type} heat pump",
            efficiency=efficiency,
            capital_cost=costs.at[costs_name, "efficiency"]
            * costs.at[costs_name, "fixed"],
            p_nom_extendable=True,
            lifetime=costs.at[costs_name, "lifetime"],
        )

        if options["tes"]:
            n.add("Carrier", name + " water tanks")

            n.madd(
                "Bus",
                nodes[name] + f" {name} water tanks",
                location=nodes[name],
                carrier=name + " water tanks",
                unit="MWh_th",
            )

            n.madd(
                "Link",
                nodes[name] + f" {name} water tanks charger",
                bus0=nodes[name] + f" {name} heat",
                bus1=nodes[name] + f" {name} water tanks",
                efficiency=costs.at["water tank charger", "efficiency"],
                carrier=name + " water tanks charger",
                p_nom_extendable=True,
            )

            n.madd(
                "Link",
                nodes[name] + f" {name} water tanks discharger",
                bus0=nodes[name] + f" {name} water tanks",
                bus1=nodes[name] + f" {name} heat",
                carrier=name + " water tanks discharger",
                efficiency=costs.at["water tank discharger", "efficiency"],
                p_nom_extendable=True,
            )

            if isinstance(options["tes_tau"], dict):
                tes_time_constant_days = options["tes_tau"][name_type]
            else:
                logger.warning(
                    "Deprecated: a future version will require you to specify 'tes_tau' ",
                    "for 'decentral' and 'central' separately.",
                )
                tes_time_constant_days = (
                    options["tes_tau"] if name_type == "decentral" else 180.0
                )

            n.madd(
                "Store",
                nodes[name] + f" {name} water tanks",
                bus=nodes[name] + f" {name} water tanks",
                e_cyclic=True,
                e_nom_extendable=True,
                carrier=name + " water tanks",
                standing_loss=1 - np.exp(-1 / 24 / tes_time_constant_days),
                capital_cost=costs.at[name_type + " water tank storage", "fixed"],
                lifetime=costs.at[name_type + " water tank storage", "lifetime"],
            )

        if options["boilers"]:
            key = f"{name_type} resistive heater"

            n.madd(
                "Link",
                nodes[name] + f" {name} resistive heater",
                bus0=nodes[name],
                bus1=nodes[name] + f" {name} heat",
                carrier=name + " resistive heater",
                efficiency=costs.at[key, "efficiency"],
                capital_cost=costs.at[key, "efficiency"] * costs.at[key, "fixed"],
                p_nom_extendable=True,
                lifetime=costs.at[key, "lifetime"],
            )

            key = f"{name_type} gas boiler"

            n.madd(
                "Link",
                nodes[name] + f" {name} gas boiler",
                p_nom_extendable=True,
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name] + f" {name} heat",
                bus2="co2 atmosphere",
                carrier=name + " gas boiler",
                efficiency=costs.at[key, "efficiency"],
                efficiency2=costs.at["gas", "CO2 intensity"],
                capital_cost=costs.at[key, "efficiency"] * costs.at[key, "fixed"],
                lifetime=costs.at[key, "lifetime"],
            )

        if options["solar_thermal"]:
            n.add("Carrier", name + " solar thermal")

            n.madd(
                "Generator",
                nodes[name],
                suffix=f" {name} solar thermal collector",
                bus=nodes[name] + f" {name} heat",
                carrier=name + " solar thermal",
                p_nom_extendable=True,
                capital_cost=costs.at[name_type + " solar thermal", "fixed"],
                p_max_pu=solar_thermal[nodes[name]],
                lifetime=costs.at[name_type + " solar thermal", "lifetime"],
            )

        if options["chp"] and name == "urban central":
            # add gas CHP; biomass CHP is added in biomass section
            n.madd(
                "Link",
                nodes[name] + " urban central gas CHP",
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name],
                bus2=nodes[name] + " urban central heat",
                bus3="co2 atmosphere",
                carrier="urban central gas CHP",
                p_nom_extendable=True,
                capital_cost=costs.at["central gas CHP", "fixed"]
                * costs.at["central gas CHP", "efficiency"],
                marginal_cost=costs.at["central gas CHP", "VOM"],
                efficiency=costs.at["central gas CHP", "efficiency"],
                efficiency2=costs.at["central gas CHP", "efficiency"]
                / costs.at["central gas CHP", "c_b"],
                efficiency3=costs.at["gas", "CO2 intensity"],
                lifetime=costs.at["central gas CHP", "lifetime"],
            )

            n.madd(
                "Link",
                nodes[name] + " urban central gas CHP CC",
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name],
                bus2=nodes[name] + " urban central heat",
                bus3="co2 atmosphere",
                bus4=spatial.co2.df.loc[nodes[name], "nodes"].values,
                carrier="urban central gas CHP CC",
                p_nom_extendable=True,
                capital_cost=costs.at["central gas CHP", "fixed"]
                * costs.at["central gas CHP", "efficiency"]
                + costs.at["biomass CHP capture", "fixed"]
                * costs.at["gas", "CO2 intensity"],
                marginal_cost=costs.at["central gas CHP", "VOM"],
                efficiency=costs.at["central gas CHP", "efficiency"]
                - costs.at["gas", "CO2 intensity"]
                * (
                    costs.at["biomass CHP capture", "electricity-input"]
                    + costs.at["biomass CHP capture", "compression-electricity-input"]
                ),
                efficiency2=costs.at["central gas CHP", "efficiency"]
                / costs.at["central gas CHP", "c_b"]
                + costs.at["gas", "CO2 intensity"]
                * (
                    costs.at["biomass CHP capture", "heat-output"]
                    + costs.at["biomass CHP capture", "compression-heat-output"]
                    - costs.at["biomass CHP capture", "heat-input"]
                ),
                efficiency3=costs.at["gas", "CO2 intensity"]
                * (1 - costs.at["biomass CHP capture", "capture_rate"]),
                efficiency4=costs.at["gas", "CO2 intensity"]
                * costs.at["biomass CHP capture", "capture_rate"],
                lifetime=costs.at["central gas CHP", "lifetime"],
            )

        if options["chp"] and options["micro_chp"] and name != "urban central":
            n.madd(
                "Link",
                nodes[name] + f" {name} micro gas CHP",
                p_nom_extendable=True,
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name],
                bus2=nodes[name] + f" {name} heat",
                bus3="co2 atmosphere",
                carrier=name + " micro gas CHP",
                efficiency=costs.at["micro CHP", "efficiency"],
                efficiency2=costs.at["micro CHP", "efficiency-heat"],
                efficiency3=costs.at["gas", "CO2 intensity"],
                capital_cost=costs.at["micro CHP", "fixed"],
                lifetime=costs.at["micro CHP", "lifetime"],
            )

    if options["retrofitting"]["retro_endogen"]:
        logger.info("Add retrofitting endogenously")

        # resample heat demand temporal 'heat_demand_r' depending on in config
        # specified temporal resolution, to not overestimate retrofitting
        hours = list(filter(re.compile(r"^\d+h$", re.IGNORECASE).search, opts))
        if len(hours) == 0:
            hours = [n.snapshots[1] - n.snapshots[0]]
        heat_demand_r = heat_demand.resample(hours[0]).mean()

        # retrofitting data 'retro_data' with 'costs' [EUR/m^2] and heat
        # demand 'dE' [per unit of original heat demand] for each country and
        # different retrofitting strengths [additional insulation thickness in m]
        retro_data = pd.read_csv(
            snakemake.input.retro_cost_energy,
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
            w_space[sector] = heat_demand_r[sector + " space"] / (
                heat_demand_r[sector + " space"] + heat_demand_r[sector + " water"]
            )
        w_space["tot"] = (
            heat_demand_r["services space"] + heat_demand_r["residential space"]
        ) / heat_demand_r.groupby(level=[1], axis=1).sum()

        for name in n.loads[
            n.loads.carrier.isin([x + " heat" for x in heat_systems])
        ].index:
            node = n.buses.loc[name, "location"]
            ct = pop_layout.loc[node, "ct"]

            # weighting 'f' depending on the size of the population at the node
            f = urban_fraction[node] if "urban" in name else (1 - urban_fraction[node])
            if f == 0:
                continue
            # get sector name ("residential"/"services"/or both "tot" for urban central)
            sec = [x if x in name else "tot" for x in sectors][0]

            # get floor aread at node and region (urban/rural) in m^2
            floor_area_node = (
                pop_layout.loc[node].fraction * floor_area.loc[ct, "value"] * 10**6
            ).loc[sec] * f
            # total heat demand at node [MWh]
            demand = n.loads_t.p_set[name].resample(hours[0]).mean()

            # space heat demand at node [MWh]
            space_heat_demand = demand * w_space[sec][node]
            # normed time profile of space heat demand 'space_pu' (values between 0-1),
            # p_max_pu/p_min_pu of retrofitting generators
            space_pu = (space_heat_demand / space_heat_demand.max()).to_frame(name=node)

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
                n.madd(
                    "Generator",
                    [node],
                    suffix=" retrofitting " + strength + " " + name[6::],
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
    central_fraction = options["district_heating"]["potential"]
    # district heating share at each node
    dist_fraction_node = (
        district_heat_share * pop_layout["urban_ct_fraction"] / pop_layout["fraction"]
    )
    nodes["urban central"] = dist_fraction_node.index
    # if district heating share larger than urban fraction -> set urban
    # fraction to district heating share
    urban_fraction = pd.concat([urban_fraction, dist_fraction_node], axis=1).max(axis=1)
    # difference of max potential and today's share of district heating
    diff = (urban_fraction * central_fraction) - dist_fraction_node
    progress = get(options["district_heating"]["progress"], investment_year)
    dist_fraction_node += diff * progress
    logger.info(
        f"Increase district heating share by a progress factor of {progress:.2%} "
        f"resulting in new average share of {dist_fraction_node.mean():.2%}"
    )

    return nodes, dist_fraction_node, urban_fraction


def add_biomass(n, costs):
    logger.info("Add biomass")

    biomass_potentials = pd.read_csv(snakemake.input.biomass_potentials, index_col=0)

    # need to aggregate potentials if gas not nodally resolved
    if options["gas_network"]:
        biogas_potentials_spatial = biomass_potentials["biogas"].rename(
            index=lambda x: x + " biogas"
        )
    else:
        biogas_potentials_spatial = biomass_potentials["biogas"].sum()

    if options.get("biomass_spatial", options["biomass_transport"]):
        solid_biomass_potentials_spatial = biomass_potentials["solid biomass"].rename(
            index=lambda x: x + " solid biomass"
        )
    else:
        solid_biomass_potentials_spatial = biomass_potentials["solid biomass"].sum()

    n.add("Carrier", "biogas")
    n.add("Carrier", "solid biomass")

    n.madd(
        "Bus",
        spatial.gas.biogas,
        location=spatial.gas.locations,
        carrier="biogas",
        unit="MWh_LHV",
    )

    n.madd(
        "Bus",
        spatial.biomass.nodes,
        location=spatial.biomass.locations,
        carrier="solid biomass",
        unit="MWh_LHV",
    )

    n.madd(
        "Store",
        spatial.gas.biogas,
        bus=spatial.gas.biogas,
        carrier="biogas",
        e_nom=biogas_potentials_spatial,
        marginal_cost=costs.at["biogas", "fuel"],
        e_initial=biogas_potentials_spatial,
    )

    n.madd(
        "Store",
        spatial.biomass.nodes,
        bus=spatial.biomass.nodes,
        carrier="solid biomass",
        e_nom=solid_biomass_potentials_spatial,
        marginal_cost=costs.at["solid biomass", "fuel"],
        e_initial=solid_biomass_potentials_spatial,
    )

    n.madd(
        "Link",
        spatial.gas.biogas_to_gas,
        bus0=spatial.gas.biogas,
        bus1=spatial.gas.nodes,
        bus2="co2 atmosphere",
        carrier="biogas to gas",
        capital_cost=costs.loc["biogas upgrading", "fixed"],
        marginal_cost=costs.loc["biogas upgrading", "VOM"],
        efficiency2=-costs.at["gas", "CO2 intensity"],
        p_nom_extendable=True,
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

        n.madd(
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

        n.madd(
            "Generator",
            spatial.biomass.nodes,
            bus=spatial.biomass.nodes,
            carrier="solid biomass",
            p_nom=10000,
            marginal_cost=costs.at["solid biomass", "fuel"]
            + bus_transport_costs * average_distance,
        )

    # AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if not urban_central.empty and options["chp"]:
        urban_central = urban_central.str[: -len(" urban central heat")]

        key = "central solid biomass CHP"

        n.madd(
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

        n.madd(
            "Link",
            urban_central + " urban central solid biomass CHP CC",
            bus0=spatial.biomass.df.loc[urban_central, "nodes"].values,
            bus1=urban_central,
            bus2=urban_central + " urban central heat",
            bus3="co2 atmosphere",
            bus4=spatial.co2.df.loc[urban_central, "nodes"].values,
            carrier="urban central solid biomass CHP CC",
            p_nom_extendable=True,
            capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"]
            + costs.at["biomass CHP capture", "fixed"]
            * costs.at["solid biomass", "CO2 intensity"],
            marginal_cost=costs.at[key, "VOM"],
            efficiency=costs.at[key, "efficiency"]
            - costs.at["solid biomass", "CO2 intensity"]
            * (
                costs.at["biomass CHP capture", "electricity-input"]
                + costs.at["biomass CHP capture", "compression-electricity-input"]
            ),
            efficiency2=costs.at[key, "efficiency-heat"]
            + costs.at["solid biomass", "CO2 intensity"]
            * (
                costs.at["biomass CHP capture", "heat-output"]
                + costs.at["biomass CHP capture", "compression-heat-output"]
                - costs.at["biomass CHP capture", "heat-input"]
            ),
            efficiency3=-costs.at["solid biomass", "CO2 intensity"]
            * costs.at["biomass CHP capture", "capture_rate"],
            efficiency4=costs.at["solid biomass", "CO2 intensity"]
            * costs.at["biomass CHP capture", "capture_rate"],
            lifetime=costs.at[key, "lifetime"],
        )

    if options["biomass_boiler"]:
        # TODO: Add surcharge for pellets
        nodes_heat = create_nodes_for_heat_sector()[0]
        for name in [
            "residential rural",
            "services rural",
            "residential urban decentral",
            "services urban decentral",
        ]:
            n.madd(
                "Link",
                nodes_heat[name] + f" {name} biomass boiler",
                p_nom_extendable=True,
                bus0=spatial.biomass.df.loc[nodes_heat[name], "nodes"].values,
                bus1=nodes_heat[name] + f" {name} heat",
                carrier=name + " biomass boiler",
                efficiency=costs.at["biomass boiler", "efficiency"],
                capital_cost=costs.at["biomass boiler", "efficiency"]
                * costs.at["biomass boiler", "fixed"],
                lifetime=costs.at["biomass boiler", "lifetime"],
            )

    # Solid biomass to liquid fuel
    if options["biomass_to_liquid"]:
        add_carrier_buses(n, "oil")

        n.madd(
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
            capital_cost=costs.at["BtL", "fixed"],
            marginal_cost=costs.loc["BtL", "VOM"] / costs.at["BtL", "efficiency"],
        )

        # TODO: Update with energy penalty
        n.madd(
            "Link",
            spatial.biomass.nodes,
            suffix=" biomass to liquid CC",
            bus0=spatial.biomass.nodes,
            bus1=spatial.oil.nodes,
            bus2="co2 atmosphere",
            bus3=spatial.co2.nodes,
            carrier="biomass to liquid",
            lifetime=costs.at["BtL", "lifetime"],
            efficiency=costs.at["BtL", "efficiency"],
            efficiency2=-costs.at["solid biomass", "CO2 intensity"]
            + costs.at["BtL", "CO2 stored"] * (1 - costs.at["BtL", "capture rate"]),
            efficiency3=costs.at["BtL", "CO2 stored"] * costs.at["BtL", "capture rate"],
            p_nom_extendable=True,
            capital_cost=costs.at["BtL", "fixed"]
            + costs.at["biomass CHP capture", "fixed"] * costs.at["BtL", "CO2 stored"],
            marginal_cost=costs.loc["BtL", "VOM"] / costs.at["BtL", "efficiency"],
        )

    # biomethanol

    if options.get("biomass_to_methanol"):
        add_carrier_buses(n, "methanol")

        n.madd(
            "Link",
            spatial.biomass.nodes,
            suffix=" biomass to methanol",
            bus0=spatial.biomass.nodes,
            bus1=spatial.methanol.nodes,
            bus2="co2 atmosphere",
            carrier="biomass to methanol",
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

        n.madd(
            "Link",
            spatial.biomass.nodes,
            suffix=" biomass to methanol CC",
            bus0=spatial.biomass.nodes,
            bus1=spatial.methanol.nodes,
            bus2="co2 atmosphere",
            bus3=spatial.co2.nodes,
            carrier="biomass to methanol CC",
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

    # BioSNG from solid biomass
    if options["biosng"]:
        n.madd(
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
            capital_cost=costs.at["BioSNG", "fixed"],
            marginal_cost=costs.at["BioSNG", "efficiency"] * costs.loc["BioSNG", "VOM"],
        )

        # TODO: Update with energy penalty for CC
        n.madd(
            "Link",
            spatial.biomass.nodes,
            suffix=" solid biomass to gas CC",
            bus0=spatial.biomass.nodes,
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            bus3="co2 atmosphere",
            carrier="BioSNG",
            lifetime=costs.at["BioSNG", "lifetime"],
            efficiency=costs.at["BioSNG", "efficiency"],
            efficiency2=costs.at["BioSNG", "CO2 stored"]
            * costs.at["BioSNG", "capture rate"],
            efficiency3=-costs.at["solid biomass", "CO2 intensity"]
            + costs.at["BioSNG", "CO2 stored"]
            * (1 - costs.at["BioSNG", "capture rate"]),
            p_nom_extendable=True,
            capital_cost=costs.at["BioSNG", "fixed"]
            + costs.at["biomass CHP capture", "fixed"]
            * costs.at["BioSNG", "CO2 stored"],
            marginal_cost=costs.at["BioSNG", "efficiency"] * costs.loc["BioSNG", "VOM"],
        )


def add_industry(n, costs):
    logger.info("Add industrial demand")

    nodes = pop_layout.index
    nhours = n.snapshot_weightings.generators.sum()
    nyears = nhours / 8760

    add_carrier_buses(n, "oil")
    add_carrier_buses(n, "methanol")

    # 1e6 to convert TWh to MWh
    industrial_demand = (
        pd.read_csv(snakemake.input.industrial_demand, index_col=[0, 1]) * 1e6 * nyears
    )
    industrial_demand_today = (
        pd.read_csv(snakemake.input.industrial_demand_today, index_col=0) * 1e6 * nyears
    )

    industrial_production = (
        pd.read_csv(snakemake.input.industrial_production, index_col=0)
        * 1e3
        * nyears  # kt/a -> t/a
    )

    endogenous_sectors = []
    if options["endogenous_steel"]:
        endogenous_sectors += ["DRI + Electric arc"]
    if options["endogenous_hvc"]:
        endogenous_sectors += ["HVC"]
    sectors_b = ~industrial_demand.index.get_level_values("sector").isin(
        endogenous_sectors
    )

    if options["endogenous_steel"]:
        logger.info("Adding endogenous primary steel demand in tonnes.")

        sector = "DRI + Electric arc"

        n.add(
            "Bus",
            "EU steel",
            location="EU",
            carrier="steel",
            unit="t",
        )

        n.add(
            "Load",
            "EU steel",
            bus="EU steel",
            carrier="steel",
            p_set=industrial_production[sector].sum() / nhours,
        )

        # n.add(
        #     "Store",
        #     "EU steel Store",
        #     bus="EU steel",
        #     e_nom_extendable=True,
        #     e_cyclic=True,
        #     carrier="steel",
        # )

        electricity_input = (
            costs.at["direct iron reduction furnace", "electricity-input"]
            * costs.at["electric arc furnace", "hbi-input"]
            + costs.at["electric arc furnace", "electricity-input"]
        )

        hydrogen_input = (
            costs.at["direct iron reduction furnace", "hydrogen-input"]
            * costs.at["electric arc furnace", "hbi-input"]
        )

        # so that for each region supply matches consumption
        p_nom = industrial_production[sector] * electricity_input / nhours

        marginal_cost = (
            costs.at["iron ore DRI-ready", "commodity"]
            * costs.at["direct iron reduction furnace", "ore-input"]
            * costs.at["electric arc furnace", "hbi-input"]
            / electricity_input
        )

        capital_cost = (
            costs.at["direct iron reduction furnace", "fixed"]
            + costs.at["electric arc furnace", "fixed"]
        ) / electricity_input

        n.madd(
            "Link",
            nodes,
            suffix=f" {sector}",
            carrier=sector,
            capital_cost=capital_cost,
            marginal_cost=marginal_cost,
            p_nom=p_nom,
            p_min_pu=1,
            bus0=nodes,
            bus1="EU steel",
            bus2=nodes + " H2",
            efficiency=1 / electricity_input,
            efficiency2=-hydrogen_input / electricity_input,
        )

    n.madd(
        "Bus",
        spatial.biomass.industry,
        location=spatial.biomass.locations,
        carrier="solid biomass for industry",
        unit="MWh_LHV",
    )

    if options.get("biomass_spatial", options["biomass_transport"]):
        p_set = (
            industrial_demand.loc[
                (spatial.biomass.locations, sectors_b), "solid biomass"
            ]
            .groupby(level="nodes")
            .sum()
            .rename(index=lambda x: x + " solid biomass for industry")
            / nhours
        )
    else:
        p_set = industrial_demand.loc[sectors_b, "solid biomass"].sum() / nhours

    n.madd(
        "Load",
        spatial.biomass.industry,
        bus=spatial.biomass.industry,
        carrier="solid biomass for industry",
        p_set=p_set,
    )

    n.madd(
        "Link",
        spatial.biomass.industry,
        bus0=spatial.biomass.nodes,
        bus1=spatial.biomass.industry,
        carrier="solid biomass for industry",
        p_nom_extendable=True,
        efficiency=1.0,
    )

    n.madd(
        "Link",
        spatial.biomass.industry_cc,
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

    n.madd(
        "Bus",
        spatial.gas.industry,
        location=spatial.gas.locations,
        carrier="gas for industry",
        unit="MWh_LHV",
    )

    gas_demand = (
        industrial_demand.loc[(nodes, sectors_b), "methane"].groupby(level="node").sum()
        / nhours
    )

    if options["gas_network"]:
        spatial_gas_demand = gas_demand.rename(index=lambda x: x + " gas for industry")
    else:
        spatial_gas_demand = gas_demand.sum()

    n.madd(
        "Load",
        spatial.gas.industry,
        bus=spatial.gas.industry,
        carrier="gas for industry",
        p_set=spatial_gas_demand,
    )

    n.madd(
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

    n.madd(
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

    n.madd(
        "Load",
        nodes,
        suffix=" H2 for industry",
        bus=nodes + " H2",
        carrier="H2 for industry",
        p_set=industrial_demand.loc[(nodes, sectors_b), "hydrogen"]
        .groupby(level="node")
        .sum()
        / nhours,
    )

    if options["oil_boilers"]:
        nodes_heat = create_nodes_for_heat_sector()[0]

        for name in [
            "residential rural",
            "services rural",
            "residential urban decentral",
            "services urban decentral",
        ]:
            n.madd(
                "Link",
                nodes_heat[name] + f" {name} oil boiler",
                p_nom_extendable=True,
                bus0=spatial.oil.nodes,
                bus1=nodes_heat[name] + f" {name}  heat",
                bus2="co2 atmosphere",
                carrier=f"{name} oil boiler",
                efficiency=costs.at["decentral oil boiler", "efficiency"],
                efficiency2=costs.at["oil", "CO2 intensity"],
                capital_cost=costs.at["decentral oil boiler", "efficiency"]
                * costs.at["decentral oil boiler", "fixed"],
                lifetime=costs.at["decentral oil boiler", "lifetime"],
            )

    n.madd(
        "Link",
        nodes + " Fischer-Tropsch",
        bus0=nodes + " H2",
        bus1=spatial.oil.nodes,
        bus2=spatial.co2.nodes,
        carrier="Fischer-Tropsch",
        efficiency=costs.at["Fischer-Tropsch", "efficiency"],
        capital_cost=costs.at["Fischer-Tropsch", "fixed"]
        * costs.at["Fischer-Tropsch", "efficiency"],  # EUR/MW_H2/a
        efficiency2=-costs.at["oil", "CO2 intensity"]
        * costs.at["Fischer-Tropsch", "efficiency"],
        p_nom_extendable=True,
        p_min_pu=options.get("min_part_load_fischer_tropsch", 0),
        lifetime=costs.at["Fischer-Tropsch", "lifetime"],
    )

    demand_factor = options.get("HVC_demand_factor", 1)
    if demand_factor != 1:
        logger.warning(f"Changing HVC demand by {demand_factor*100-100:+.2f}%.")

    if options["endogenous_hvc"]:
        logger.info("Adding endogenous primary high-value chemicals demand in tonnes.")

        n.add(
            "Bus",
            "EU HVC",
            location="EU",
            carrier="HVC",
            unit="t",
        )

        p_set = demand_factor * industrial_production["HVC"].sum() / nhours

        n.add(
            "Load",
            "EU HVC",
            bus="EU HVC",
            carrier="HVC",
            p_set=p_set,
        )

        # n.add(
        #     "Store",
        #     "EU HVC Store",
        #     bus="EU HVC",
        #     e_nom_extendable=True,
        #     e_cyclic=True,
        #     carrier="HVC",
        # )

        tech = "methanol-to-olefins/aromatics"

        logger.info(f"Adding {tech}.")

        p_nom_max = (
            demand_factor
            * industrial_production.loc[nodes, "HVC"]
            / nhours
            * costs.at[tech, "methanol-input"]
        )

        co2_release = (
            costs.at[tech, "carbondioxide-output"] / costs.at[tech, "methanol-input"]
            + costs.at["methanolisation", "carbondioxide-input"]
        )

        n.madd(
            "Link",
            nodes,
            suffix=f" {tech}",
            carrier=tech,
            capital_cost=costs.at[tech, "fixed"] / costs.at[tech, "methanol-input"],
            marginal_cost=costs.at[tech, "VOM"] / costs.at[tech, "methanol-input"],
            p_nom_extendable=True,
            bus0=spatial.methanol.nodes,
            bus1="EU HVC",
            bus2=nodes,
            bus3="co2 atmosphere",
            p_min_pu=1,
            p_nom_max=p_nom_max.values,
            efficiency=1 / costs.at[tech, "methanol-input"],
            efficiency2=-costs.at[tech, "electricity-input"]
            / costs.at[tech, "methanol-input"],
            efficiency3=co2_release,
        )

        tech = "electric steam cracker"

        logger.info(f"Adding {tech}.")

        p_nom_max = (
            demand_factor
            * industrial_production.loc[nodes, "HVC"]
            / nhours
            * costs.at[tech, "naphtha-input"]
        )

        co2_release = (
            costs.at[tech, "carbondioxide-output"] / costs.at[tech, "naphtha-input"]
            + costs.at["oil", "CO2 intensity"]
        )

        n.madd(
            "Link",
            nodes,
            suffix=f" {tech}",
            carrier=tech,
            capital_cost=costs.at[tech, "fixed"] / costs.at[tech, "naphtha-input"],
            marginal_cost=costs.at[tech, "VOM"] / costs.at[tech, "naphtha-input"],
            p_nom_extendable=True,
            bus0=spatial.oil.nodes,
            bus1="EU HVC",
            bus2=nodes,
            bus3="co2 atmosphere",
            p_min_pu=1,
            p_nom_max=p_nom_max.values,
            efficiency=1 / costs.at[tech, "naphtha-input"],
            efficiency2=-costs.at[tech, "electricity-input"]
            / costs.at[tech, "naphtha-input"],
            efficiency3=co2_release,
        )

    p_set = (
        demand_factor
        * industrial_demand.loc[(nodes, sectors_b), "naphtha"].sum()
        / nhours
    )

    n.madd(
        "Load",
        ["naphtha for industry"],
        bus=spatial.oil.nodes,
        carrier="naphtha for industry",
        p_set=p_set,
    )

    demand_factor = options.get("aviation_demand_factor", 1)
    all_aviation = ["total international aviation", "total domestic aviation"]
    p_set = (
        demand_factor
        * pop_weighted_energy_totals.loc[nodes, all_aviation].sum(axis=1).sum()
        * 1e6
        / nhours
    )
    if demand_factor != 1:
        logger.warning(f"Changing aviation demand by {demand_factor*100-100:+.2f}%.")

    n.madd(
        "Load",
        ["kerosene for aviation"],
        bus=spatial.oil.nodes,
        carrier="kerosene for aviation",
        p_set=p_set,
    )

    # NB: CO2 gets released again to atmosphere when plastics decay or kerosene is burned
    # except for the process emissions when naphtha is used for petrochemicals, which can be captured with other industry process emissions
    # tco2 per hour
    co2_release = ["naphtha for industry", "kerosene for aviation"]
    co2 = (
        n.loads.loc[co2_release, "p_set"].sum() * costs.at["oil", "CO2 intensity"]
        - industrial_demand.loc[
            (nodes, sectors_b), "process emission from feedstock"
        ].sum()
        / nhours
    )

    n.add(
        "Load",
        "oil emissions",
        bus="co2 atmosphere",
        carrier="oil emissions",
        p_set=-co2,
    )

    # industry methanol demand

    p_set_methanol = (
        industrial_demand.loc[(nodes, sectors_b), "methanol"].sum() / nhours
    )

    n.madd(
        "Load",
        spatial.methanol.locations,
        suffix=" industry methanol",
        bus=spatial.methanol.nodes,
        carrier="industry methanol",
        p_set=p_set_methanol,
    )

    co2 = p_set_methanol * costs.at["methanolisation", "carbondioxide-input"]

    n.add(
        "Load",
        "industry methanol emissions",
        bus="co2 atmosphere",
        carrier="industry methanol emissions",
        p_set=-co2,
    )

    # methanolisation

    n.madd(
        "Link",
        spatial.h2.locations + " methanolisation",
        bus0=spatial.h2.nodes,
        bus1=spatial.methanol.nodes,
        bus2=nodes,
        bus3=spatial.co2.nodes,
        carrier="methanolisation",
        p_nom_extendable=True,
        p_min_pu=options.get("min_part_load_methanolisation", 0),
        capital_cost=costs.at["methanolisation", "fixed"]
        / costs.at["methanolisation", "hydrogen-input"],  # EUR/MW_H2/a
        lifetime=costs.at["methanolisation", "lifetime"],
        efficiency=1 / costs.at["methanolisation", "hydrogen-input"],
        efficiency2=-costs.at["methanolisation", "electricity-input"]
        / costs.at["methanolisation", "hydrogen-input"],
        efficiency3=-costs.at["methanolisation", "carbondioxide-input"]
        / costs.at["methanolisation", "hydrogen-input"],
    )

    # TODO simplify bus expression
    n.madd(
        "Load",
        nodes,
        suffix=" low-temperature heat for industry",
        bus=[
            node + " urban central heat"
            if node + " urban central heat" in n.buses.index
            else node + " services urban decentral heat"
            for node in nodes
        ],
        carrier="low-temperature heat for industry",
        p_set=industrial_demand.loc[(nodes, sectors_b), "low-temperature heat"]
        .groupby(level="node")
        .sum()
        / nhours,
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
            - industrial_demand_today.loc[loads_i, "electricity"].sum()
            / n.loads_t.p_set[loads_i].sum().sum()
        )
        n.loads_t.p_set[loads_i] *= factor

    n.madd(
        "Load",
        nodes,
        suffix=" industry electricity",
        bus=nodes,
        carrier="industry electricity",
        p_set=industrial_demand.loc[(nodes, sectors_b), "electricity"]
        .groupby(level="node")
        .sum()
        / nhours,
    )

    n.madd(
        "Bus",
        spatial.co2.process_emissions,
        location=spatial.co2.locations,
        carrier="process emissions",
        unit="t_co2",
    )

    sel = ["process emission", "process emission from feedstock"]
    if options["co2_spatial"] or options["co2network"]:
        p_set = (
            -industrial_demand.loc[(nodes, sectors_b), sel]
            .groupby(level="node")
            .sum()
            .sum(axis=1)
            .rename(index=lambda x: x + " process emissions")
            / nhours
        )
    else:
        p_set = (
            -industrial_demand.loc[(nodes, sectors_b), sel].sum(axis=1).sum() / nhours
        )

    # this should be process emissions fossil+feedstock
    # then need load on atmosphere for feedstock emissions that are currently going to atmosphere via Link Fischer-Tropsch demand
    n.madd(
        "Load",
        spatial.co2.process_emissions,
        bus=spatial.co2.process_emissions,
        carrier="process emissions",
        p_set=p_set,
    )

    n.madd(
        "Link",
        spatial.co2.process_emissions,
        bus0=spatial.co2.process_emissions,
        bus1="co2 atmosphere",
        carrier="process emissions",
        p_nom_extendable=True,
        efficiency=1.0,
    )

    # assume enough local waste heat for CC
    n.madd(
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

    if options.get("ammonia"):
        if options["ammonia"] == "regional":
            p_set = (
                industrial_demand.loc[(spatial.ammonia.locations, sectors_b), "ammonia"]
                .groupby(level="node")
                .sum()
                .rename(index=lambda x: x + " NH3")
                / nhours
            )
        else:
            p_set = industrial_demand.loc[sectors_b, "ammonia"].sum() / nhours

        n.madd(
            "Load",
            spatial.ammonia.nodes,
            bus=spatial.ammonia.nodes,
            carrier="NH3",
            p_set=p_set,
        )


def add_shipping(n, costs):
    logger.info("Add shipping demand")

    nodes = pop_layout.index

    shipping_hydrogen_share = get(options["shipping_hydrogen_share"], investment_year)
    shipping_methanol_share = get(options["shipping_methanol_share"], investment_year)
    shipping_oil_share = get(options["shipping_oil_share"], investment_year)

    total_share = shipping_hydrogen_share + shipping_methanol_share + shipping_oil_share
    if total_share != 1:
        logger.warning(
            f"Total shipping shares sum up to {total_share:.2%}, corresponding to increased or decreased demand assumptions."
        )

    domestic_navigation = pop_weighted_energy_totals.loc[
        nodes, "total domestic navigation"
    ].squeeze()
    international_navigation = (
        pd.read_csv(snakemake.input.shipping_demand, index_col=0).squeeze() * nyears
    )
    all_navigation = domestic_navigation + international_navigation
    p_set = all_navigation * 1e6 / nhours

    if shipping_hydrogen_share:
        oil_efficiency = options.get(
            "shipping_oil_efficiency", options.get("shipping_average_efficiency", 0.4)
        )
        efficiency = oil_efficiency / options.get("shipping_lh2_efficiency", 0.44)
        shipping_hydrogen_share = get(
            options["shipping_hydrogen_share"], investment_year
        )

        if options["shipping_hydrogen_liquefaction"]:
            n.madd(
                "Bus",
                nodes,
                suffix=" H2 liquid",
                carrier="H2 liquid",
                location=nodes,
                unit="MWh_LHV",
            )

            n.madd(
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

        efficiency = options["shipping_oil_efficiency"] / options.get(
            "shipping_lh2_efficiency", 0.44
        )
        p_set_hydrogen = shipping_hydrogen_share * p_set * efficiency

        n.madd(
            "Load",
            nodes,
            suffix=" H2 for shipping",
            bus=shipping_bus,
            carrier="H2 for shipping",
            p_set=p_set_hydrogen,
        )

    if shipping_methanol_share:
        add_carrier_buses(n, "methanol")

        efficiency = (
            options["shipping_oil_efficiency"] / options["shipping_methanol_efficiency"]
        )
        p_set_methanol = shipping_methanol_share * p_set.sum() * efficiency

        n.madd(
            "Load",
            spatial.methanol.nodes,
            suffix=" shipping methanol",
            bus=spatial.methanol.nodes,
            carrier="shipping methanol",
            p_set=p_set_methanol,
        )

        # CO2 intensity methanol based on stoichiometric calculation with 22.7 GJ/t methanol (32 g/mol), CO2 (44 g/mol), 277.78 MWh/TJ = 0.218 t/MWh
        co2 = p_set_methanol * costs.at["methanolisation", "carbondioxide-input"]

        n.add(
            "Load",
            "shipping methanol emissions",
            bus="co2 atmosphere",
            carrier="shipping methanol emissions",
            p_set=-co2,
        )

    if shipping_oil_share:
        add_carrier_buses(n, "oil")

        p_set_oil = shipping_oil_share * p_set.sum()

        n.madd(
            "Load",
            spatial.oil.nodes,
            suffix=" shipping oil",
            bus=spatial.oil.nodes,
            carrier="shipping oil",
            p_set=p_set_oil,
        )

        co2 = p_set_oil * costs.at["oil", "CO2 intensity"]

        n.add(
            "Load",
            "shipping oil emissions",
            bus="co2 atmosphere",
            carrier="shipping oil emissions",
            p_set=-co2,
        )


def add_shipping_endogenous(n, costs):
    fuels = options["shipping_endogenous"]["fuels"]

    logger.info(f"Add shipping demand endogenously with options {fuels}")

    nodes = pop_layout.index

    domestic_navigation = pop_weighted_energy_totals.loc[
        nodes, "total domestic navigation"
    ].squeeze()
    international_navigation = (
        pd.read_csv(snakemake.input.shipping_demand, index_col=0).squeeze() * nyears
    )
    all_navigation = domestic_navigation + international_navigation

    # determine demand in MW propulsion
    oil_efficiency = options["shipping_oil_efficiency"]
    p_set = oil_efficiency * all_navigation * 1e6 / nhours

    n.madd(
        "Bus",
        nodes,
        suffix=" shipping",
        unit="MWh propulsion",
        location=nodes,
        carrier="shipping",
    )

    n.madd(
        "Load",
        nodes,
        suffix=" shipping",
        bus=nodes + " shipping",
        carrier="shipping",
        p_set=p_set,
    )

    if "oil" in fuels:
        add_carrier_buses(n, "oil")

        n.madd(
            "Link",
            nodes,
            suffix=" shipping oil",
            carrier="shipping oil",
            bus0=spatial.oil.nodes,
            bus1=nodes + " shipping",
            bus2="co2 atmosphere",
            efficiency=options["shipping_oil_efficiency"],
            efficiency2=costs.at["oil", "CO2 intensity"],
            p_min_pu=1,
            p_nom_extendable=True,
        )

    if "methanol" in fuels:
        add_carrier_buses(n, "methanol")

        n.madd(
            "Link",
            nodes,
            suffix=" shipping methanol",
            carrier="shipping methanol",
            bus0=spatial.methanol.nodes,
            bus1=nodes + " shipping",
            bus2="co2 atmosphere",
            efficiency=options["shipping_methanol_efficiency"],
            efficiency2=costs.at["methanolisation", "carbondioxide-input"],
            p_min_pu=1,
            p_nom_extendable=True,
        )

    if "LH2" in fuels:
        n.madd(
            "Link",
            nodes,
            suffix=" shipping LH2",
            bus0=nodes + " H2",
            bus1=nodes + " shipping",
            carrier="shipping LH2",
            efficiency=options["shipping_lh2_efficiency"]
            * costs.at["H2 liquefaction", "efficiency"],
            capital_cost=costs.at["H2 liquefaction", "fixed"],
            lifetime=costs.at["H2 liquefaction", "lifetime"],
            p_min_pu=1,
            p_nom_extendable=True,
        )

    if "ammonia" in fuels:
        if not options["ammonia"]:
            logger.warning(
                "Ammonia not represented as carrier. Skipping it as shipping fuel option."
            )

        n.madd(
            "Link",
            nodes,
            suffix=" shipping ammonia",
            bus0=spatial.ammonia.nodes,
            bus1=nodes + " shipping",
            carrier="shipping ammonia",
            efficiency=options["shipping_ammonia_efficiency"],
            p_min_pu=1,
            p_nom_extendable=True,
        )

    if "LNG" in fuels:
        n.madd(
            "Link",
            nodes,
            suffix=" shipping LNG",
            bus0=spatial.gas.nodes,
            bus1=nodes + " shipping",
            bus2="co2 atmosphere",
            carrier="shipping LNG",
            efficiency=options["shipping_lng_efficiency"],
            efficiency2=costs.at["gas", "CO2 intensity"],
            capital_cost=costs.at["CH4 liquefaction", "fixed"],
            lifetime=costs.at["CH4 liquefaction", "lifetime"],
            p_min_pu=1,
            p_nom_extendable=True,
        )


def add_waste_heat(n):
    # TODO options?

    logger.info("Add possibility to use industrial waste heat in district heating")

    # AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if not urban_central.empty:
        urban_central = urban_central.str[: -len(" urban central heat")]

        # TODO what is the 0.95 and should it be a config option?
        if options["use_fischer_tropsch_waste_heat"]:
            n.links.loc[urban_central + " Fischer-Tropsch", "bus3"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " Fischer-Tropsch", "efficiency3"] = (
                0.95 - n.links.loc[urban_central + " Fischer-Tropsch", "efficiency"]
            )

        if options["use_methanation_waste_heat"]:
            n.links.loc[urban_central + " Sabatier", "bus3"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " Sabatier", "efficiency3"] = (
                0.95 - n.links.loc[urban_central + " Sabatier", "efficiency"]
            )

        if options["use_methanolisation_waste_heat"]:
            n.links.loc[urban_central + " methanolisation", "bus4"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " methanolisation", "efficiency4"] = (
                costs.at["methanolisation", "heat-output"]
                / costs.at["methanolisation", "hydrogen-input"]
            )

        # TODO integrate usable waste heat efficiency into technology-data from DEA
        if options.get("use_electrolysis_waste_heat", False):
            n.links.loc[urban_central + " H2 Electrolysis", "bus2"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " H2 Electrolysis", "efficiency2"] = (
                0.84 - n.links.loc[urban_central + " H2 Electrolysis", "efficiency"]
            )

        if options["use_fuel_cell_waste_heat"]:
            n.links.loc[urban_central + " H2 Fuel Cell", "bus2"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " H2 Fuel Cell", "efficiency2"] = (
                0.95 - n.links.loc[urban_central + " H2 Fuel Cell", "efficiency"]
            )


def add_agriculture(n, costs):
    logger.info("Add agriculture, forestry and fishing sector.")

    nodes = pop_layout.index
    nhours = n.snapshot_weightings.generators.sum()

    # electricity

    n.madd(
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

    n.madd(
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

    machinery_nodal_energy = pop_weighted_energy_totals.loc[
        nodes, "total agriculture machinery"
    ]

    if electric_share > 0:
        efficiency_gain = (
            options["agriculture_machinery_fuel_efficiency"]
            / options["agriculture_machinery_electric_efficiency"]
        )

        n.madd(
            "Load",
            nodes,
            suffix=" agriculture machinery electric",
            bus=nodes,
            carrier="agriculture machinery electric",
            p_set=electric_share
            / efficiency_gain
            * machinery_nodal_energy
            * 1e6
            / nhours,
        )

    if oil_share > 0:
        add_carrier_buses(n, "oil")

        n.madd(
            "Load",
            ["agriculture machinery oil"],
            bus=spatial.oil.nodes,
            carrier="agriculture machinery oil",
            p_set=oil_share * machinery_nodal_energy.sum() * 1e6 / nhours,
        )

        co2 = (
            oil_share
            * machinery_nodal_energy.sum()
            * 1e6
            / nhours
            * costs.at["oil", "CO2 intensity"]
        )

        n.add(
            "Load",
            "agriculture machinery oil emissions",
            bus="co2 atmosphere",
            carrier="agriculture machinery oil emissions",
            p_set=-co2,
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


def add_endogenous_hvdc_import_options(n, cost_factor=1.0):
    logger.info("Add import options: endogenous hvdc-to-elec")
    cf = snakemake.config["sector"]["import"].get("endogenous_hvdc_import", {})
    if not cf["enable"]:
        return

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    p_max_pu = xr.open_dataset(snakemake.input.import_p_max_pu).p_max_pu.sel(
        importer="EUE"
    )

    p_nom_max = (
        xr.open_dataset(snakemake.input.import_p_max_pu)
        .p_nom_max.sel(importer="EUE")
        .to_pandas()
    )

    def _coordinates(ct):
        iso2 = ct.split("-")[0]
        if iso2 in country_centroids.index:
            return country_centroids.loc[iso2, ["longitude", "latitude"]].values
        else:
            query = cc.convert(iso2, to="name")
            loc = geocode(dict(country=query), language="en")
            return [loc.longitude, loc.latitude]

    exporters = pd.DataFrame(
        {ct: _coordinates(ct) for ct in cf["exporters"]}, index=["x", "y"]
    ).T
    geometry = gpd.points_from_xy(exporters.x, exporters.y)
    exporters = gpd.GeoDataFrame(exporters, geometry=geometry, crs=4326)

    import_links = {}
    a = regions.representative_point().to_crs(DISTANCE_CRS)
    for ct in exporters.index:
        b = exporters.to_crs(DISTANCE_CRS).loc[ct].geometry
        d = a.distance(b)
        import_links[ct] = (
            d.where(d < d.quantile(cf["distance_threshold"])).div(1e3).dropna()
        )  # km
    import_links = pd.concat(import_links)

    # xlinks
    xlinks = {}
    for bus0, links in cf["xlinks"].items():
        for link in links:
            landing_point = gpd.GeoSeries(
                [Point(link["x"], link["y"])], crs=4326
            ).to_crs(DISTANCE_CRS)
            bus1 = (
                regions.to_crs(DISTANCE_CRS)
                .geometry.distance(landing_point[0])
                .idxmin()
            )
            xlinks[(bus0, bus1)] = link["length"]

    import_links = pd.concat([import_links, pd.Series(xlinks)], axis=0)
    import_links = import_links.drop_duplicates(keep="first")

    hvdc_cost = (
        import_links.values * cf["length_factor"] * costs.at["HVDC submarine", "fixed"]
        + costs.at["HVDC inverter pair", "fixed"]
    )

    buses_i = exporters.index

    n.madd("Bus", buses_i, **exporters.drop("geometry", axis=1))

    efficiency = cf["efficiency_static"] * cf["efficiency_per_1000km"] ** (
        import_links.values / 1e3
    )

    n.madd(
        "Link",
        ["import hvdc-to-elec " + " ".join(idx).strip() for idx in import_links.index],
        bus0=import_links.index.get_level_values(0),
        bus1=import_links.index.get_level_values(1),
        carrier="import hvdc-to-elec",
        p_min_pu=0,
        p_nom_extendable=True,
        length=import_links.values,
        capital_cost=hvdc_cost * cost_factor,
        efficiency=efficiency,
        p_nom_max=cf["p_nom_max"],
    )

    for tech in ["solar-utility", "onwind"]:
        p_max_pu_tech = p_max_pu.sel(technology=tech).to_pandas().dropna().T

        exporters_tech_i = exporters.index.intersection(p_max_pu_tech.columns)

        grid_connection = costs.at["electricity grid connection", "fixed"]

        n.madd(
            "Generator",
            exporters_tech_i,
            suffix=" " + tech,
            bus=exporters_tech_i,
            carrier=f"external {tech}",
            p_nom_extendable=True,
            capital_cost=(costs.at[tech, "fixed"] + grid_connection) * cost_factor,
            lifetime=costs.at[tech, "lifetime"],
            p_max_pu=p_max_pu_tech.reindex(columns=exporters_tech_i),
            p_nom_max=p_nom_max[tech].reindex(index=exporters_tech_i)
            * cf["share_of_p_nom_max_available"],
        )

    # hydrogen storage

    h2_buses_i = n.madd(
        "Bus", buses_i, suffix=" H2", carrier="external H2", location=buses_i
    )

    n.madd(
        "Store",
        h2_buses_i,
        bus=h2_buses_i,
        carrier="external H2",
        e_nom_extendable=True,
        e_cyclic=True,
        capital_cost=costs.at[
            "hydrogen storage tank type 1 including compressor", "fixed"
        ]
        * cost_factor,
    )

    n.madd(
        "Link",
        h2_buses_i + " Electrolysis",
        bus0=buses_i,
        bus1=h2_buses_i,
        carrier="external H2 Electrolysis",
        p_nom_extendable=True,
        efficiency=costs.at["electrolysis", "efficiency"],
        capital_cost=costs.at["electrolysis", "fixed"] * cost_factor,
        lifetime=costs.at["electrolysis", "lifetime"],
    )

    n.madd(
        "Link",
        h2_buses_i + " H2 Turbine",
        bus0=h2_buses_i,
        bus1=buses_i,
        carrier="external H2 Turbine",
        p_nom_extendable=True,
        efficiency=costs.at["OCGT", "efficiency"],
        capital_cost=costs.at["OCGT", "fixed"]
        * costs.at["OCGT", "efficiency"]
        * cost_factor,
        lifetime=costs.at["OCGT", "lifetime"],
    )

    # battery storage

    b_buses_i = n.madd(
        "Bus", buses_i, suffix=" battery", carrier="external battery", location=buses_i
    )

    n.madd(
        "Store",
        b_buses_i,
        bus=b_buses_i,
        carrier="external battery",
        e_cyclic=True,
        e_nom_extendable=True,
        capital_cost=costs.at["battery storage", "fixed"] * cost_factor,
        lifetime=costs.at["battery storage", "lifetime"],
    )

    n.madd(
        "Link",
        b_buses_i + " charger",
        bus0=buses_i,
        bus1=b_buses_i,
        carrier="external battery charger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        capital_cost=costs.at["battery inverter", "fixed"] * cost_factor,
        p_nom_extendable=True,
        lifetime=costs.at["battery inverter", "lifetime"],
    )

    n.madd(
        "Link",
        b_buses_i + " discharger",
        bus0=b_buses_i,
        bus1=buses_i,
        carrier="external battery discharger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        p_nom_extendable=True,
        lifetime=costs.at["battery inverter", "lifetime"],
    )

    # add extra HVDC connections between MENA countries

    for bus0_bus1 in cf.get("extra_connections", []):
        bus0, bus1 = bus0_bus1.split("-")

        a = exporters.to_crs(DISTANCE_CRS).at[bus0, "geometry"]
        b = exporters.to_crs(DISTANCE_CRS).at[bus1, "geometry"]
        d = a.distance(b) / 1e3  # km

        capital_cost = (
            d * cf["length_factor"] * costs.at["HVDC overhead", "fixed"]
            + costs.at["HVDC inverter pair", "fixed"]
        )

        n.add(
            "Link",
            f"external HVDC {bus0_bus1}",
            bus0=bus0,
            bus1=bus1,
            carrier="external HVDC",
            p_min_pu=-1,
            p_nom_extendable=True,
            capital_cost=capital_cost * cost_factor,
            length=d,
        )


def add_import_options(
    n,
    capacity_boost=3.0,
    import_options=[
        "hvdc-to-elec",
        "pipeline-h2",
        "shipping-lh2",
        "shipping-lch4",
        "shipping-meoh",
        "shipping-ftfuel",
        "shipping-lnh3",
        "shipping-steel",
    ],
    endogenous_hvdc=False,
):
    if not isinstance(import_options, dict):
        import_options = {k: 1.0 for k in import_options}

    logger.info("Add import options: " + " ".join(import_options.keys()))
    fn = snakemake.input.gas_input_nodes_simplified
    import_nodes = pd.read_csv(fn, index_col=0)
    import_nodes["hvdc-to-elec"] = 15000

    import_config = snakemake.config["sector"]["import"]

    ports = pd.read_csv(snakemake.input.ports, index_col=0)

    translate = {
        "pipeline-h2": "pipeline",
        "hvdc-to-elec": "hvdc-to-elec",
        "shipping-lh2": "lng",
        "shipping-lch4": "lng",
        "shipping-lnh3": "lng",
    }

    bus_suffix = {
        "pipeline-h2": " H2",
        "hvdc-to-elec": "",
        "shipping-lh2": " H2",
        "shipping-lch4": " gas",
        "shipping-lnh3": " NH3",
        "shipping-ftfuel": " oil",
        "shipping-meoh": " methanol",
        "shipping-steel": " steel",
    }

    co2_intensity = {
        "shipping-lch4": ("gas", "CO2 intensity"),
        "shipping-ftfuel": ("oil", "CO2 intensity"),
        "shipping-meoh": ("methanolisation", "carbondioxide-input"),
    }

    import_costs = pd.read_csv(snakemake.input.import_costs, delimiter=";")
    cols = ["esc", "exporter", "importer", "value"]
    fields = ["Cost per MWh delivered", "Cost per t delivered"]
    import_costs = import_costs.query("subcategory in @fields")[cols]
    import_costs.rename(columns={"value": "marginal_cost"}, inplace=True)

    for k, v in translate.items():
        import_nodes[k] = import_nodes[v]
        ports[k] = ports.get(v)

    if endogenous_hvdc and "hvdc-to-elec" in import_options:
        add_endogenous_hvdc_import_options(n, import_options.pop("hvdc-to-elec"))

    regionalised_options = {
        "hvdc-to-elec",
        "pipeline-h2",
        "shipping-lh2",
        "shipping-lch4",
    }

    for tech in set(import_options).intersection(regionalised_options):
        import_costs_tech = (
            import_costs.query("esc == @tech").groupby("importer").marginal_cost.min()
        ) * import_options[tech]

        sel = ~import_nodes[tech].isna()
        if tech == "pipeline-h2":
            entrypoints_internal = ["DE", "BE", "FR", "GB"]
            entrypoints_via_RU_BY = ["EE", "LT", "LV", "FI"]  # maybe PL Yamal
            forbidden_pipelines = entrypoints_internal + entrypoints_via_RU_BY
            sel &= ~import_nodes.index.str[:2].isin(forbidden_pipelines)
        import_nodes_tech = import_nodes.loc[sel, [tech]]

        import_nodes_tech = import_nodes_tech.rename(columns={tech: "p_nom"})

        marginal_costs = ports[tech].dropna().map(import_costs_tech)
        import_nodes_tech["marginal_cost"] = marginal_costs

        import_nodes_tech.dropna(inplace=True)

        upper_p_nom_max = import_config["p_nom_max"].get(tech, np.inf)

        suffix = bus_suffix[tech]

        if tech in co2_intensity.keys():
            buses = import_nodes_tech.index + f"{suffix} import {tech}"
            capital_cost = 7018.0 if tech == "shipping-lch4" else 0.0  # €/MW/a

            n.madd("Bus", buses + " bus", carrier=f"import {tech}")

            n.madd(
                "Store",
                buses + " store",
                bus=buses + " bus",
                e_nom_extendable=True,
                e_nom_min=-np.inf,
                e_nom_max=0,
                e_min_pu=1,
                e_max_pu=0,
            )

            n.madd(
                "Link",
                buses,
                bus0=buses + " bus",
                bus1=import_nodes_tech.index + suffix,
                bus2="co2 atmosphere",
                carrier=f"import {tech}",
                efficiency2=-costs.at[co2_intensity[tech][0], co2_intensity[tech][1]],
                marginal_cost=import_nodes_tech.marginal_cost.values,
                p_nom_extendable=True,
                capital_cost=capital_cost,
                p_nom_min=import_nodes_tech.p_nom.clip(upper=upper_p_nom_max).values,
                p_nom_max=import_nodes_tech.p_nom.mul(capacity_boost)
                .clip(upper=upper_p_nom_max)
                .values,
            )

        else:
            location = import_nodes_tech.index
            buses = location if tech == "hvdc-to-elec" else location + suffix

            capital_cost = (
                1.2 * 7018.0 if tech == "shipping-lh2" else 0.0
            )  # €/MW/a, +20% compared to LNG

            n.madd(
                "Generator",
                import_nodes_tech.index + f"{suffix} import {tech}",
                bus=buses,
                carrier=f"import {tech}",
                marginal_cost=import_nodes_tech.marginal_cost.values,
                p_nom_extendable=True,
                capital_cost=capital_cost,
                p_nom_max=import_nodes_tech.p_nom.mul(capacity_boost)
                .clip(upper=upper_p_nom_max)
                .values,
            )

    # need special handling for copperplated imports
    copperplated_carbonaceous_options = {
        "shipping-ftfuel",
        "shipping-meoh",
    }

    for tech in set(import_options).intersection(copperplated_carbonaceous_options):
        marginal_costs = (
            import_costs.query("esc == @tech").marginal_cost.min()
            * import_options[tech]
        )

        suffix = bus_suffix[tech]

        n.add("Bus", f"EU import {tech} bus", carrier=f"import {tech}")

        n.add(
            "Store",
            f"EU import {tech} store",
            bus=f"EU import {tech} bus",
            e_nom_extendable=True,
            e_nom_min=-np.inf,
            e_nom_max=0,
            e_min_pu=1,
            e_max_pu=0,
        )

        n.add(
            "Link",
            f"EU import {tech}",
            bus0=f"EU import {tech} bus",
            bus1="EU" + suffix,
            bus2="co2 atmosphere",
            carrier=f"import {tech}",
            efficiency2=-costs.at[co2_intensity[tech][0], co2_intensity[tech][1]],
            marginal_cost=marginal_costs,
            p_nom=1e7,
        )

    copperplated_carbonfree_options = {
        "shipping-steel",
        "shipping-lnh3",
    }

    for tech in set(import_options).intersection(copperplated_carbonfree_options):
        suffix = bus_suffix[tech]

        marginal_costs = (
            import_costs.query("esc == @tech").marginal_cost.min()
            * import_options[tech]
        )

        n.add(
            "Generator",
            f"EU import {tech}",
            bus="EU" + suffix,
            carrier=f"import {tech}",
            marginal_cost=marginal_costs,
            p_nom=1e7,
        )


def maybe_adjust_costs_and_potentials(n, opts):
    for o in opts:
        flags = ["+e", "+p", "+m"]
        if all(flag not in o for flag in flags):
            continue
        oo = o.split("+")
        carrier_list = np.hstack(
            (
                n.generators.carrier.unique(),
                n.links.carrier.unique(),
                n.stores.carrier.unique(),
                n.storage_units.carrier.unique(),
            )
        )
        suptechs = map(lambda c: c.split("-", 2)[0], carrier_list)
        if oo[0].startswith(tuple(suptechs)):
            carrier = oo[0]
            attr_lookup = {"p": "p_nom_max", "e": "e_nom_max", "c": "capital_cost"}
            attr = attr_lookup[oo[1][0]]
            factor = float(oo[1][1:])
            # beware if factor is 0 and p_nom_max is np.inf, 0*np.inf is nan
            if carrier == "AC":  # lines do not have carrier
                n.lines[attr] *= factor
            else:
                if attr == "p_nom_max":
                    comps = {"Generator", "Link", "StorageUnit"}
                elif attr == "e_nom_max":
                    comps = {"Store"}
                else:
                    comps = {"Generator", "Link", "StorageUnit", "Store"}
                for c in n.iterate_components(comps):
                    if carrier == "solar":
                        sel = c.df.carrier.str.contains(
                            carrier
                        ) & ~c.df.carrier.str.contains("solar rooftop")
                    else:
                        sel = c.df.carrier.str.contains(carrier)
                    c.df.loc[sel, attr] *= factor
            logger.info(f"changing {attr} for {carrier} by factor {factor}")


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
        df = df.groupby(level=0).agg(agg, **agg_group_kwargs)
        # time-varying data
        pnl = c.pnl
        agg = define_clustering(pd.Index(pnl.keys()), aggregate_dict)
        for k in pnl.keys():

            def renamer(s):
                return s.replace("residential ", "").replace("services ", "")

            pnl[k] = pnl[k].groupby(renamer, axis=1).agg(agg[k], **agg_group_kwargs)

        # remove unclustered assets of service/residential
        to_drop = c.df.index.difference(df.index)
        n.mremove(c.name, to_drop)
        # add clustered assets
        to_add = df.index.difference(c.df.index)
        import_components_from_dataframe(n, df.loc[to_add], c.name)


def apply_time_segmentation(
    n, segments, solver_name="cbc", overwrite_time_dependent=True
):
    """
    Aggregating time series to segments with different lengths.

    Input:
        n: pypsa Network
        segments: (int) number of segments in which the typical period should be
                  subdivided
        solver_name: (str) name of solver
        overwrite_time_dependent: (bool) overwrite time dependent data of pypsa network
        with typical time series created by tsam
    """
    try:
        import tsam.timeseriesaggregation as tsam
    except:
        raise ModuleNotFoundError(
            "Optional dependency 'tsam' not found." "Install via 'pip install tsam'"
        )

    # get all time-dependent data
    columns = pd.MultiIndex.from_tuples([], names=["component", "key", "asset"])
    raw = pd.DataFrame(index=n.snapshots, columns=columns)
    for c in n.iterate_components():
        for attr, pnl in c.pnl.items():
            # exclude e_min_pu which is used for SOC of EVs in the morning
            if not pnl.empty and attr != "e_min_pu":
                df = pnl.copy()
                df.columns = pd.MultiIndex.from_product([[c.name], [attr], df.columns])
                raw = pd.concat([raw, df], axis=1)

    # normalise all time-dependent data
    annual_max = raw.max().replace(0, 1)
    raw = raw.div(annual_max, level=0)

    # get representative segments
    agg = tsam.TimeSeriesAggregation(
        raw,
        hoursPerPeriod=len(raw),
        noTypicalPeriods=1,
        noSegments=int(segments),
        segmentation=True,
        solver=solver_name,
    )
    segmented = agg.createTypicalPeriods()

    weightings = segmented.index.get_level_values("Segment Duration")
    offsets = np.insert(np.cumsum(weightings[:-1]), 0, 0)
    timesteps = [raw.index[0] + pd.Timedelta(f"{offset}h") for offset in offsets]
    snapshots = pd.DatetimeIndex(timesteps)
    sn_weightings = pd.Series(
        weightings, index=snapshots, name="weightings", dtype="float64"
    )
    logger.info("Distribution of snapshot durations: ", weightings.value_counts())

    n.set_snapshots(sn_weightings.index)
    n.snapshot_weightings = n.snapshot_weightings.mul(sn_weightings, axis=0)

    # overwrite time-dependent data with timeseries created by tsam
    if overwrite_time_dependent:
        values_t = segmented.mul(annual_max).set_index(snapshots)
        for component, key in values_t.columns.droplevel(2).unique():
            n.pnl(component)[key] = values_t[component, key]

    return n


def set_temporal_aggregation(n, opts, solver_name):
    """
    Aggregate network temporally.
    """
    for o in opts:
        # temporal averaging
        m = re.match(r"^\d+h$", o, re.IGNORECASE)
        if m is not None:
            n = average_every_nhours(n, m.group(0))
            break
        # representative snapshots
        m = re.match(r"(^\d+)sn$", o, re.IGNORECASE)
        if m is not None:
            sn = int(m[1])
            logger.info(f"Use every {sn} snapshot as representative")
            n.set_snapshots(n.snapshots[::sn])
            n.snapshot_weightings *= sn
            break
        # segments with package tsam
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            segments = int(m[1])
            logger.info(f"Use temporal segmentation with {segments} segments")
            n = apply_time_segmentation(n, segments, solver_name=solver_name)
            break
    return n


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
    n.links.loc[
        carrier_i, "efficiency"
    ] = efficiency_static * efficiency_per_1000km ** (
        n.links.loc[carrier_i, "length"] / 1e3
    )
    rev_links = (
        n.links.loc[carrier_i].copy().rename({"bus0": "bus1", "bus1": "bus0"}, axis=1)
    )
    rev_links.capital_cost = 0
    rev_links["reversed"] = True
    rev_links.index = rev_links.index.map(lambda x: x + "-reversed")

    n.links = pd.concat([n.links, rev_links], sort=False)
    n.links["reversed"] = n.links["reversed"].fillna(False)

    # do compression losses after concatenation to take electricity consumption at bus0 in either direction
    carrier_i = n.links.query("carrier == @carrier").index
    if compression_per_1000km > 0:
        n.links.loc[carrier_i, "bus2"] = n.links.loc[carrier_i, "bus0"].map(
            n.buses.location
        )  # electricity
        n.links.loc[carrier_i, "efficiency2"] = (
            -compression_per_1000km * n.links.loc[carrier_i, "length"] / 1e3
        )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_sector_network",
            configfiles="../../../config/config.yaml",
            simpl="",
            opts="",
            clusters="128",
            ll="v1.5",
            sector_opts="CO2L0-3H-T-H-B-I-A",
            planning_horizons="2050",
        )

    logging.basicConfig(level=snakemake.config["logging"]["level"])

    update_config_with_sector_opts(snakemake.config, snakemake.wildcards.sector_opts)

    options = snakemake.params.sector

    opts = snakemake.wildcards.sector_opts.split("-")

    investment_year = int(snakemake.wildcards.planning_horizons[-4:])

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

    country_centroids = pd.read_csv(
        snakemake.input.country_centroids[0], index_col="ISO"
    )

    patch_electricity_network(n)

    spatial = define_spatial(pop_layout.index, options)

    if snakemake.params.foresight == "myopic":
        add_lifetime_wind_solar(n, costs)

        conventional = snakemake.params.conventional_carriers
        for carrier in conventional:
            add_carrier_buses(n, carrier)

    add_co2_tracking(n, options)

    add_generation(n, costs)

    add_storage_and_grids(n, costs)

    # TODO merge with opts cost adjustment below
    for o in opts:
        if o[:4] == "wave":
            wave_cost_factor = float(o[4:].replace("p", ".").replace("m", "-"))
            logger.info(
                f"Including wave generators with cost factor of {wave_cost_factor}"
            )
            add_wave(n, wave_cost_factor)
        if o[:4] == "dist":
            options["electricity_distribution_grid"] = True
            options["electricity_distribution_grid_cost_factor"] = float(
                o[4:].replace("p", ".").replace("m", "-")
            )
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

    if options["ammonia"]:
        add_ammonia(n, costs)

    if "I" in opts:
        add_industry(n, costs)

    if "S" in opts and not options["shipping_endogenous"]["enable"]:
        add_shipping(n, costs)

    if "S" in opts and options["shipping_endogenous"]["enable"]:
        add_shipping_endogenous(n, costs)

    if "I" in opts and "H" in opts:
        add_waste_heat(n)

    if "A" in opts:  # requires H and I
        add_agriculture(n, costs)

    if options["dac"]:
        add_dac(n, costs)

    if "decentral" in opts:
        decentral(n)

    if "noH2network" in opts:
        remove_h2_network(n)

    if options["co2network"]:
        add_co2_network(n, costs)

    translate = dict(
        H2=["pipeline-h2", "shipping-lh2"],
        AC=["hvdc-to-elec"],
        CH4=["shipping-lch4"],
        NH3=["shipping-lnh3"],
        FT=["shipping-ftfuel"],
        MeOH=["shipping-meoh"],
        St=["shipping-steel"],
    )
    for o in opts:
        if not o.startswith("imp"):
            continue
        subsets = o.split("+")[1:]
        if len(subsets):

            def parse_carriers(s):
                prefixes = sorted(translate.keys(), key=lambda k: len(k), reverse=True)
                pattern = rf'({"|".join(prefixes)})(\d+(\.\d+)?)?'
                match = re.search(pattern, s)
                prefix = match.group(1) if match else None
                number = float(match.group(2)) if match and match.group(2) else 1.0
                return {prefix: number}

            carriers = {
                tk: v
                for s in subsets
                for k, v in parse_carriers(s).items()
                for tk in translate.get(k, [])
            }
        else:
            carriers = options["import"]["options"]
        add_import_options(
            n,
            capacity_boost=options["import"]["capacity_boost"],
            import_options=carriers,
            endogenous_hvdc=options["import"]["endogenous_hvdc_import"]["enable"],
        )
        break

    if options["allam_cycle_gas"]:
        add_allam_cycle_gas(n, costs)

    add_methanol_to_power(n, costs, types=options["methanol_to_power"])

    if options["methanol_to_kerosene"]:
        add_methanol_to_kerosene(n, costs)

    if options["methanol_reforming"]:
        add_methanol_reforming(n, costs)

    if options["methanol_reforming_cc"]:
        add_methanol_reforming_cc(n, costs)

    solver_name = snakemake.config["solving"]["solver"]["name"]
    n = set_temporal_aggregation(n, opts, solver_name)

    limit_type = "config"
    limit = get(snakemake.params.co2_budget, investment_year)
    for o in opts:
        if "cb" not in o:
            continue
        limit_type = "carbon budget"
        fn = "results/" + snakemake.params.RDIR + "csvs/carbon_budget_distribution.csv"
        if not os.path.exists(fn):
            emissions_scope = snakemake.params.emissions_scope
            report_year = snakemake.params.eurostat_report_year
            input_co2 = snakemake.input.co2
            build_carbon_budget(
                o,
                snakemake.input.eurostat,
                fn,
                emissions_scope,
                report_year,
                input_co2,
            )
        co2_cap = pd.read_csv(fn, index_col=0).squeeze()
        limit = co2_cap.loc[investment_year]
        break
    for o in opts:
        if "Co2L" not in o:
            continue
        limit_type = "wildcard"
        limit = o[o.find("Co2L") + 4 :]
        limit = float(limit.replace("p", ".").replace("m", "-"))
        break
    logger.info(f"Add CO2 limit from {limit_type}")
    add_co2limit(n, nyears, limit)

    for o in opts:
        if not o[:10] == "linemaxext":
            continue
        maxext = float(o[10:]) * 1e3
        limit_individual_line_extension(n, maxext)
        break

    if options["electricity_distribution_grid"]:
        insert_electricity_distribution_grid(n, costs)

    maybe_adjust_costs_and_potentials(n, opts)

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
        n.mremove("Line", idx)

    first_year_myopic = (snakemake.params.foresight == "myopic") and (
        snakemake.params.planning_horizons[0] == investment_year
    )

    if options.get("cluster_heat_buses", False) and not first_year_myopic:
        cluster_heat_buses(n)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))

    sanitize_carriers(n, snakemake.config)

    n.export_to_netcdf(snakemake.output[0])
