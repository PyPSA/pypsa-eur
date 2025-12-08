# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
# SPDX-FileCopyrightText: Open Energy Transition gGmbH
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
from networkx.algorithms import complement
from networkx.algorithms.connectivity.edge_augmentation import k_edge_augmentation
from pypsa.geo import haversine_pts
from scipy.stats import beta

from scripts._helpers import (
    configure_logging,
    get,
    load_costs,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.add_electricity import (
    calculate_annuity,
    flatten,
    sanitize_carriers,
    sanitize_locations,
)
from scripts.build_energy_totals import (
    build_co2_totals,
    build_eea_co2,
    build_eurostat,
    build_eurostat_co2,
)
from scripts.build_transport_demand import transport_degree_factor
from scripts.definitions.heat_sector import HeatSector
from scripts.definitions.heat_system import HeatSystem
from scripts.prepare_network import maybe_adjust_costs_and_potentials

spatial = SimpleNamespace()
logger = logging.getLogger(__name__)


def define_spatial(nodes, options):
    """
    Namespace for spatial.

    Parameters
    ----------
    nodes : list-like
    """

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
        spatial.oil.non_sequestered_hvc = nodes + " non-sequestered HVC"
        spatial.oil.kerosene = nodes + " kerosene for aviation"
        spatial.oil.shipping = nodes + " shipping oil"
        spatial.oil.agriculture_machinery = nodes + " agriculture machinery oil"
        spatial.oil.land_transport = nodes + " land transport oil"
    else:
        spatial.oil.demand_locations = ["EU"]
        spatial.oil.naphtha = ["EU naphtha for industry"]
        spatial.oil.non_sequestered_hvc = ["EU non-sequestered HVC"]
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
def build_carbon_budget(
    o,
    input_eurostat,
    fn,
    emissions_scope,
    input_co2,
    options,
    countries,
    planning_horizons,
):
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

    elif "ex" in o:
        T = carbon_budget / e_0
        m = (1 + np.sqrt(1 + r * T)) / T

        def exponential_decay(t):
            return (e_0 / e_1990) * (1 + (m + r) * (t - t_0)) * np.exp(-m * (t - t_0))

        co2_cap = pd.Series(
            {t: exponential_decay(t) for t in planning_horizons}, name=o
        )
    else:
        raise ValueError("Transition path must be either beta or exponential decay")

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


def haversine(p, n):
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
    profiles: dict[str, str],
    landfall_lengths: dict = None,
    line_length_factor: int | float = 1,
) -> None:
    """
    Update costs for wind and solar generators added with pypsa-eur to those
    cost in the planning year.

    Parameters
    ----------
    n : pypsa.Network
        Network to update generator costs
    costs : pd.DataFrame
        Cost assumptions DataFrame
    line_length_factor : int | float, optional
        Factor to multiply line lengths by, by default 1
    landfall_lengths : dict, optional
        Dictionary of landfall lengths per technology, by default None
    profiles : dict[str, str]
        Dictionary mapping technology names to profile file paths
        e.g. {'offwind-dc': 'path/to/profile.nc'}
    """

    if landfall_lengths is None:
        landfall_lengths = {}

    # NB: solar costs are also manipulated for rooftop
    # when distribution grid is inserted
    n.generators.loc[n.generators.carrier == "solar", "capital_cost"] = costs.at[
        "solar-utility", "capital_cost"
    ]

    n.generators.loc[n.generators.carrier == "onwind", "capital_cost"] = costs.at[
        "onwind", "capital_cost"
    ]

    # for offshore wind, need to calculated connection costs
    for key, fn in profiles.items():
        tech = key[len("profile_") :]
        landfall_length = landfall_lengths.get(tech, 0.0)

        if tech not in n.generators.carrier.values:
            continue

        with xr.open_dataset(fn) as ds:
            # if-statement for compatibility with old profiles
            if "year" in ds.indexes:
                ds = ds.sel(year=ds.year.min(), drop=True)

            ds = ds.stack(bus_bin=["bus", "bin"])

            distance = ds["average_distance"].to_pandas()
            distance.index = distance.index.map(flatten)
            submarine_cost = costs.at[tech + "-connection-submarine", "capital_cost"]
            underground_cost = costs.at[
                tech + "-connection-underground", "capital_cost"
            ]
            connection_cost = line_length_factor * (
                distance * submarine_cost + landfall_length * underground_cost
            )

            # Take 'offwind-float' capital cost for 'float', and 'offwind' capital cost for the rest ('ac' and 'dc')
            midtech = tech.split("-", 2)[1]
            if midtech == "float":
                capital_cost = (
                    costs.at[tech, "capital_cost"]
                    + costs.at[tech + "-station", "capital_cost"]
                    + connection_cost
                )
            else:
                capital_cost = (
                    costs.at["offwind", "capital_cost"]
                    + costs.at[tech + "-station", "capital_cost"]
                    + connection_cost
                )

            logger.info(
                f"Added connection cost of {connection_cost.min():0.0f}-{connection_cost.max():0.0f} Eur/MW/a to {tech}"
            )

            n.generators.loc[n.generators.carrier == tech, "capital_cost"] = (
                capital_cost.rename(index=lambda node: node + " " + tech)
            )


def add_carrier_buses(
    n: pypsa.Network,
    carrier: str,
    costs: pd.DataFrame,
    spatial: SimpleNamespace,
    options: dict,
    cf_industry: dict | None = None,
    nodes: pd.Index | list | set | None = None,
) -> None:
    """
    Add buses and associated components for a specific carrier to the network.

    Creates a new carrier type in the network and adds corresponding buses, stores,
    and potentially generators depending on the carrier type. Special handling is
    implemented for fossil fuels, particularly oil which may include refining processes.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    carrier : str
        Name of the energy carrier (e.g., 'gas', 'oil', 'coal', 'nuclear')
    costs : pd.DataFrame
        DataFrame containing cost assumptions for different technologies and fuels
    spatial : SimpleNamespace
        Namespace containing spatial information for different carriers, including
        nodes and locations
    options : dict
        Configuration dictionary, must contain 'fossil_fuels' boolean
    cf_industry : dict, optional
        Dictionary of industrial conversion factors, must contain 'oil_refining_emissions'
        if carrier is 'oil'
    nodes : pd.Index or list or set, optional
        Nodes where the carrier should be added. If None, nodes are taken from
        spatial data for the carrier

    Returns
    -------
    None
        Modifies the network object in-place by adding new components

    Notes
    -----
    - For gas carriers, energy is tracked in MWh_LHV (Lower Heating Value)
    - For other carriers, energy is tracked in MWh_th (thermal)
    - Special handling is implemented for oil refining emissions
    - Storage costs are technology-specific and based on volumetric capacity
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

    # Calculate carrier-specific storage costs
    if carrier == "gas":
        capital_cost = costs.at["gas storage", "capital_cost"]
    elif carrier == "oil":
        # based on https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html
        mwh_per_m3 = 44.9 * 724 * 0.278 * 1e-3  # MJ/kg * kg/m3 * kWh/MJ * MWh/kWh
        capital_cost = (
            costs.at["General liquid hydrocarbon storage (product)", "capital_cost"]
            / mwh_per_m3
        )
    elif carrier == "methanol":
        # based on https://www.engineeringtoolbox.com/fossil-fuels-energy-content-d_1298.html
        mwh_per_m3 = 5.54 * 791 * 1e-3  # kWh/kg * kg/m3 * MWh/kWh
        capital_cost = (
            costs.at["General liquid hydrocarbon storage (product)", "capital_cost"]
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
def remove_elec_base_techs(n: pypsa.Network, carriers_to_keep: dict) -> None:
    """
    Remove conventional generators (e.g. OCGT) and storage units (e.g.
    batteries and H2) from base electricity-only network, since they're added
    here differently using links.

    Parameters
    ----------
    n : pypsa.Network
        Network to remove components from
    carriers_to_keep : dict
        Dictionary specifying which carriers to keep for each component type
        e.g. {'Generator': ['hydro'], 'StorageUnit': ['PHS']}
    """
    for c in n.iterate_components(carriers_to_keep):
        to_keep = carriers_to_keep[c.name]
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


def patch_electricity_network(n, costs, carriers_to_keep, profiles, landfall_lengths):
    remove_elec_base_techs(n, carriers_to_keep)
    remove_non_electric_buses(n)
    update_wind_solar_costs(
        n, costs, landfall_lengths=landfall_lengths, profiles=profiles
    )
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


def add_co2_tracking(
    n, costs, options, sequestration_potential_file=None, co2_price: float = 0.0
):
    """
    Add CO2 tracking components to the network including atmospheric CO2,
    CO2 storage, and sequestration infrastructure.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        Cost assumptions for different technologies, must include
        'CO2 storage tank' with 'capital_cost' column
    options : dict
        Configuration options containing at least:
        - regional_co2_sequestration_potential: dict with keys
            - enable: bool
            - max_size: float
            - years_of_storage: float
        - co2_sequestration_cost: float
        - co2_sequestration_lifetime: float
        - co2_vent: bool
    sequestration_potential_file : str, optional
        Path to CSV file containing regional CO2 sequestration potentials.
        Required if options['regional_co2_sequestration_potential']['enable'] is True.
    co2_price : float, optional
        CO2 price that needs to be paid for emitting into the atmosphere and which is
        gained by removing from the atmosphere.

    Returns
    -------
    None
        Modifies the network object in-place by adding CO2-related components.

    Notes
    -----
    Adds several components to track CO2:
    - Atmospheric CO2 store
    - CO2 storage tanks
    - CO2 sequestration infrastructure
    - Optional CO2 venting facilities
    """
    # minus sign because opposite to how fossil fuels used:
    # CH4 burning puts CH4 down, atmosphere up
    n.add("Carrier", "co2", co2_emissions=-1.0)

    # this tracks CO2 in the atmosphere
    n.add("Bus", "co2 atmosphere", location="EU", carrier="co2", unit="t_co2")

    # can also be negative
    n.add(
        "Store",
        "co2 atmosphere",
        e_nom=np.inf,
        e_min_pu=-1,
        carrier="co2",
        bus="co2 atmosphere",
        marginal_cost=-co2_price,
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
        capital_cost=costs.at["CO2 storage tank", "capital_cost"],
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
        if sequestration_potential_file is None:
            raise ValueError(
                "sequestration_potential_file must be provided when "
                "regional_co2_sequestration_potential is enabled"
            )
        upper_limit = (
            options["regional_co2_sequestration_potential"]["max_size"] * 1e3
        )  # Mt
        annualiser = options["regional_co2_sequestration_potential"]["years_of_storage"]
        df = pd.read_csv(sequestration_potential_file, index_col=0)
        if df.shape == (1, 1):
            # if only one value, manually convert to a Series
            e_nom_max = pd.Series(df.iloc[0, 0], index=df.index)
        else:
            e_nom_max = df.squeeze()

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


def add_co2_network(n, costs, co2_network_cost_factor=1.0):
    """
    Add CO2 transport network to the PyPSA network.

    Creates a CO2 pipeline network with both onshore and submarine pipeline segments,
    considering different costs for each type. The network allows bidirectional flow
    and is extendable.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        Cost assumptions for different technologies. Must contain entries for
        'CO2 pipeline' and 'CO2 submarine pipeline' with 'capital_cost' and 'lifetime'
        columns
    co2_network_cost_factor : float, optional
        Factor to scale the capital costs of the CO2 network, default 1.0

    Returns
    -------
    None
        Modifies the network object in-place by adding CO2 pipeline links

    Notes
    -----
    The function creates bidirectional CO2 pipeline links between nodes, with costs
    depending on the underwater fraction of the pipeline. The network topology is
    created using the create_network_topology helper function.
    """
    logger.info("Adding CO2 network.")
    co2_links = create_network_topology(n, "CO2 pipeline ")

    if "underwater_fraction" not in co2_links.columns:
        co2_links["underwater_fraction"] = 0.0

    cost_onshore = (
        (1 - co2_links.underwater_fraction)
        * costs.at["CO2 pipeline", "capital_cost"]
        * co2_links.length
    )
    cost_submarine = (
        co2_links.underwater_fraction
        * costs.at["CO2 submarine pipeline", "capital_cost"]
        * co2_links.length
    )
    capital_cost = cost_onshore + cost_submarine
    capital_cost *= co2_network_cost_factor

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


def add_allam_gas(
    n: pypsa.Network,
    costs: pd.DataFrame,
    pop_layout: pd.DataFrame,
    spatial: SimpleNamespace,
) -> None:
    """
    Add Allam cycle gas power plants to the network as Link components.

    Allam cycle plants are modeled as links with four buses:
    - Input: natural gas
    - Output: electricity
    - Output: CO2 for storage/usage (98% of emissions)
    - Output: CO2 to atmosphere (2% of emissions)

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        Costs and parameters for different technologies. Must contain 'allam' and 'gas'
        entries with columns for 'fixed', 'VOM', 'efficiency', 'lifetime', and
        'CO2 intensity' parameters
    pop_layout : pd.DataFrame
        DataFrame containing population layout data with nodes as index
    spatial : SimpleNamespace
        Container for spatial data with attributes:
        - gas.df: DataFrame with gas network nodes
        - co2.df: DataFrame with CO2 network nodes
        Both DataFrames must have a 'nodes' column

    Returns
    -------
    None
        Modifies the network object in-place by adding Allam cycle plants as Links

    Notes
    -----
    The Allam cycle is a novel natural gas power plant design with integrated
    carbon capture. It captures approximately 98% of CO2 emissions, with 2%
    going to the atmosphere.
    """
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
        capital_cost=costs.at["allam", "capital_cost"]
        * costs.at["allam", "efficiency"],
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
        capital_cost=costs.at["biomass-to-methanol", "capital_cost"]
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
        capital_cost=costs.at["biomass-to-methanol", "capital_cost"]
        / costs.at["biomass-to-methanol", "efficiency"]
        + costs.at["biomass CHP capture", "capital_cost"]
        * costs.at["biomass-to-methanol", "CO2 stored"],
        marginal_cost=costs.loc["biomass-to-methanol", "VOM"]
        / costs.at["biomass-to-methanol", "efficiency"],
    )


def add_methanol_to_power(n, costs, pop_layout, types=None):
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
            capital_cost=costs.at["allam", "capital_cost"]
            * costs.at["allam", "efficiency"],
            marginal_cost=costs.at["allam", "VOM"] * costs.at["allam", "efficiency"],
            efficiency=costs.at["allam", "efficiency"],
            efficiency2=0.98 * costs.at["methanolisation", "carbondioxide-input"],
            efficiency3=0.02 * costs.at["methanolisation", "carbondioxide-input"],
            lifetime=25,
        )

    if types["ccgt"]:
        logger.info("Adding methanol CCGT power plants.")

        # efficiency * EUR/MW * (annuity + FOM)
        capital_cost = costs.at["CCGT", "efficiency"] * costs.at["CCGT", "capital_cost"]

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
        capital_cost = costs.at["CCGT", "efficiency"] * costs.at["CCGT", "capital_cost"]

        capital_cost_cc = (
            capital_cost
            + costs.at["cement capture", "capital_cost"]
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
            capital_cost=costs.at["OCGT", "capital_cost"]
            * costs.at["OCGT", "efficiency"],
            marginal_cost=costs.at["OCGT", "VOM"] * costs.at["OCGT", "efficiency"],
            efficiency=costs.at["OCGT", "efficiency"],
            efficiency2=costs.at["methanolisation", "carbondioxide-input"],
            lifetime=costs.at["OCGT", "lifetime"],
        )


def add_methanol_reforming(n, costs):
    logger.info("Adding methanol steam reforming.")

    tech = "Methanol steam reforming"

    capital_cost = costs.at[tech, "capital_cost"] / costs.at[tech, "methanol-input"]

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

    capital_cost = costs.at[tech, "capital_cost"] / costs.at[tech, "methanol-input"]

    capital_cost_cc = (
        capital_cost
        + costs.at["cement capture", "capital_cost"]
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
        capital_cost=costs.at["direct air capture", "capital_cost"] / electricity_input,
        efficiency=-heat_input / electricity_input,
        efficiency2=-1 / electricity_input,
        efficiency3=1 / electricity_input,
        p_nom_extendable=True,
        lifetime=costs.at["direct air capture", "lifetime"],
    )


def add_co2limit(n, options, co2_totals_file, countries, nyears, limit):
    """
    Add a global CO2 emissions constraint to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    options : dict
        Dictionary of options determining which sectors to consider for emissions
    co2_totals_file : str
        Path to CSV file containing historical CO2 emissions data in Mt
        (megatonnes) per country and sector
    countries : list
        List of country codes to consider for the CO2 limit
    nyears : float, optional
        Number of years for the CO2 budget, by default 1.0
    limit : float, optional
        CO2 limit as a fraction of 1990 levels

    Returns
    -------
    None
        The function modifies the network object in-place by adding a global
        CO2 constraint.

    Notes
    -----
    The function reads historical CO2 emissions data, calculates a total limit
    based on the specified countries and sectors, and adds a global constraint
    to the network. The limit is calculated as a fraction of historical emissions
    multiplied by the number of years.
    """
    if limit is None:
        return

    logger.info(f"Adding CO2 budget limit as per unit of 1990 levels of {limit}")

    sectors = determine_emission_sectors(options)

    # convert Mt to tCO2
    co2_totals = 1e6 * pd.read_csv(co2_totals_file, index_col=0)

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


def add_generation(
    n: pypsa.Network,
    costs: pd.DataFrame,
    pop_layout: pd.DataFrame,
    conventionals: dict[str, str],
    spatial: SimpleNamespace,
    options: dict,
    cf_industry: dict,
    ext_carriers,
    existing_capacities=None,
    existing_efficiencies=None,
) -> None:
    """
    Add conventional electricity generation to the network.

    Creates links between carrier buses and demand nodes for conventional generators,
    including their efficiency, costs, and CO2 emissions.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        DataFrame containing cost and technical parameters for different technologies
    pop_layout : pd.DataFrame
        DataFrame with population layout data, used for demand nodes
    conventionals : Dict[str, str]
        Dictionary mapping generator types to their energy carriers
        e.g., {'OCGT': 'gas', 'CCGT': 'gas', 'coal': 'coal'}
    spatial : SimpleNamespace
        Namespace containing spatial information for different carriers,
        including nodes and locations
    options : dict
        Configuration dictionary containing settings for the model
    cf_industry : dict
        Dictionary of industrial conversion factors, needed for carrier buses

    Returns
    -------
    None
        Modifies the network object in-place by adding generation components

    Notes
    -----
    - Costs (VOM and fixed) are given per MWel and automatically adjusted by efficiency
    - CO2 emissions are tracked through a link to the 'co2 atmosphere' bus
    - Generator lifetimes are considered in the capital cost calculation
    """
    logger.info("Adding electricity generation")

    nodes = pop_layout.index

    for generator, carrier in conventionals.items():
        carrier_nodes = vars(spatial)[carrier].nodes

        add_carrier_buses(
            n=n,
            carrier=carrier,
            costs=costs,
            spatial=spatial,
            options=options,
            cf_industry=cf_industry,
        )

        n.add(
            "Link",
            nodes + " " + generator,
            bus0=carrier_nodes,
            bus1=nodes,
            bus2="co2 atmosphere",
            marginal_cost=costs.at[generator, "efficiency"]
            * costs.at[generator, "VOM"],  # NB: VOM is per MWel
            capital_cost=costs.at[generator, "efficiency"]
            * costs.at[generator, "capital_cost"],  # NB: fixed cost is per MWel
            p_nom_extendable=bool(generator in ext_carriers.get("Generator", [])),
            p_nom=(
                existing_capacities[generator] / existing_efficiencies[generator]
                if existing_capacities is not None
                else 0
            ),  # NB: existing capacities are MWel
            p_max_pu=0.7
            if carrier == "uranium"
            else 1,  # be conservative for nuclear (maintenance or unplanned shut downs)
            p_nom_min=(
                existing_capacities[generator] if existing_capacities is not None else 0
            ),
            carrier=generator,
            efficiency=(
                existing_efficiencies[generator]
                if existing_efficiencies is not None
                else costs.at[generator, "efficiency"]
            ),
            efficiency2=costs.at[carrier, "CO2 intensity"],
            lifetime=costs.at[generator, "lifetime"],
        )


def add_ammonia(
    n: pypsa.Network,
    costs: pd.DataFrame,
    pop_layout: pd.DataFrame,
    spatial: SimpleNamespace,
    cf_industry: dict,
) -> None:
    """
    Add ammonia synthesis, cracking, and storage infrastructure to the network.

    Creates the necessary components for an ammonia economy including Haber-Bosch
    synthesis plants, ammonia crackers, and storage facilities. Links are created
    between electricity, hydrogen, and ammonia buses.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        Technology cost assumptions with MultiIndex columns containing
        'fixed', 'VOM', 'efficiency', 'lifetime', etc.
    pop_layout : pd.DataFrame
        Population layout data with index of location nodes
    spatial : Namespace
        Configuration object containing ammonia-specific spatial information
        with attributes:
        - nodes: list of ammonia bus nodes
        - locations: list of geographical locations
    cf_industry : dict
        Industry-specific conversion factors including
        'MWh_NH3_per_MWh_H2_cracker' for ammonia cracking efficiency
    logger : logging.Logger
        Logger object for output messages

    Returns
    -------
    None
        Modifies the network object in-place by adding ammonia-related components

    Notes
    -----
    The function adds several components:
    - NH3 carrier
    - Ammonia buses at specified locations
    - Haber-Bosch synthesis plants linking electricity, hydrogen, and ammonia
    - Ammonia crackers for converting back to hydrogen
    - Ammonia storage facilities
    """
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
        capital_cost=costs.at["Haber-Bosch", "capital_cost"]
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
        capital_cost=costs.at["Ammonia cracker", "capital_cost"]
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
        capital_cost=costs.at[
            "NH3 (l) storage tank incl. liquefaction", "capital_cost"
        ],
        lifetime=costs.at["NH3 (l) storage tank incl. liquefaction", "lifetime"],
    )


def insert_electricity_distribution_grid(
    n: pypsa.Network,
    costs: pd.DataFrame,
    options: dict,
    pop_layout: pd.DataFrame,
    solar_rooftop_potentials_fn: str,
) -> None:
    """
    Insert electricity distribution grid components into the network.

    Adds low voltage buses, distribution grid links, rooftop solar potential,
    and home battery storage systems to the network. Also adjusts the connection
    points of various loads and distributed energy resources to the low voltage grid.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object to be modified
    costs : pd.DataFrame
        Technology cost assumptions with technologies as index and cost parameters
        as columns, including 'fixed' costs, 'lifetime', and component-specific
        parameters like 'efficiency'
    options : dict
        Configuration options containing at least:
        - transmission_efficiency: dict with distribution grid parameters
        - marginal_cost_storage: float for storage operation costs
    pop_layout : pd.DataFrame
        Population data per node with at least:
        - 'total' column containing population in thousands
        Index should match network nodes

    Returns
    -------
    None
        Modifies the network object in-place by adding components

    Notes
    -----
    Components added to the network:
    - Low voltage buses for each node
    - Distribution grid links connecting high to low voltage
    - Rooftop solar potential based on population density
    - Home battery storage systems with separate charger/discharger links

    The function also adjusts the connection points of loads like:
    - Regular electricity demand
    - Electric vehicles (BEV chargers and V2G)
    - Heat pumps
    - Resistive heaters
    - Micro-CHP units
    """
    nodes = n.buses.query("carrier == 'AC'").index

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
        capital_cost=costs.at["electricity distribution grid", "capital_cost"],
    )

    # deduct distribution losses from electricity demand as these are included in total load
    # https://nbviewer.org/github/Open-Power-System-Data/datapackage_timeseries/blob/2020-10-06/main.ipynb
    if (
        efficiency := options["transmission_efficiency"]
        .get("electricity distribution grid", {})
        .get("efficiency_static")
    ) and "electricity distribution grid" in options["transmission_efficiency"][
        "enable"
    ]:
        logger.info(
            f"Deducting distribution losses from electricity demand: {np.around(100 * (1 - efficiency), decimals=2)}%"
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
    n.links.loc[hps, "bus1"] += " low voltage"

    rh = n.links.index[n.links.carrier.str.contains("resistive heater")]
    n.links.loc[rh, "bus0"] += " low voltage"

    mchp = n.links.index[n.links.carrier.str.contains("micro gas")]
    n.links.loc[mchp, "bus1"] += " low voltage"

    # set existing solar to cost of utility cost rather the 50-50 rooftop-utility
    solar = n.generators.index[n.generators.carrier == "solar"]
    n.generators.loc[solar, "capital_cost"] = costs.at["solar-utility", "capital_cost"]

    fn = solar_rooftop_potentials_fn
    if len(fn) > 0:
        potential = pd.read_csv(fn, index_col=["bus", "bin"]).squeeze(axis=1)
        potential.index = potential.index.map(flatten) + " solar"

        n.add(
            "Generator",
            solar,
            suffix=" rooftop",
            bus=n.generators.loc[solar, "bus"] + " low voltage",
            carrier="solar rooftop",
            p_nom_extendable=True,
            p_nom_max=potential.loc[solar],
            marginal_cost=n.generators.loc[solar, "marginal_cost"],
            capital_cost=costs.at["solar-rooftop", "capital_cost"],
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
        capital_cost=costs.at["home battery storage", "capital_cost"],
        lifetime=costs.at["battery storage", "lifetime"],
    )

    n.add(
        "Link",
        nodes + " home battery charger",
        bus0=nodes + " low voltage",
        bus1=nodes + " home battery",
        carrier="home battery charger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        capital_cost=costs.at["home battery inverter", "capital_cost"],
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
        marginal_cost=costs.at["home battery storage", "marginal_cost"],
        p_nom_extendable=True,
        lifetime=costs.at["battery inverter", "lifetime"],
    )


def insert_gas_distribution_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    options: dict,
) -> None:
    """
    Insert gas distribution grid costs into gas-consuming components.

    Adds distribution grid costs to decentralized gas boilers and micro-CHP units
    by increasing their capital costs. The additional cost is calculated as a factor
    of electricity distribution grid costs.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object to be modified
    costs : pd.DataFrame
        Technology cost assumptions with technologies as index and cost parameters
        as columns, must include 'electricity distribution grid' with 'fixed' costs
    options : dict
        Configuration options containing at least:
        - gas_distribution_grid_cost_factor: float
          Factor to multiply electricity distribution grid costs by

    Returns
    -------
    None
        Modifies the network object in-place by updating capital costs of gas
        components

    Notes
    -----
    The function adds distribution grid costs to:
    - Decentralized gas boilers (excluding urban central heating)
    - Micro-CHP units

    The additional cost is calculated by multiplying the electricity distribution
    grid fixed cost by the gas distribution grid cost factor.
    """
    f_costs = options["gas_distribution_grid_cost_factor"]

    logger.info(
        f"Inserting gas distribution grid with investment cost factor of {f_costs}"
    )

    capital_cost = costs.at["electricity distribution grid", "capital_cost"] * f_costs

    # Add costs to decentralized gas boilers
    gas_b = n.links.index[
        n.links.carrier.str.contains("gas boiler")
        & (~n.links.carrier.str.contains("urban central"))
    ]
    n.links.loc[gas_b, "capital_cost"] += capital_cost

    # Add costs to micro CHPs
    mchp = n.links.index[n.links.carrier.str.contains("micro gas")]
    n.links.loc[mchp, "capital_cost"] += capital_cost


def add_electricity_grid_connection(n, costs):
    carriers = ["onwind", "solar", "solar-hsat"]

    gens = n.generators.index[n.generators.carrier.isin(carriers)]

    n.generators.loc[gens, "capital_cost"] += costs.at[
        "electricity grid connection", "capital_cost"
    ]


def add_storage_and_grids(
    n,
    costs,
    pop_layout,
    h2_cavern_file,
    cavern_types,
    clustered_gas_network_file,
    gas_input_nodes,
    spatial,
    options,
):
    """
    Add storage and grid infrastructure to the network including hydrogen, gas, and battery systems.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        Technology cost assumptions
    pop_layout : pd.DataFrame
        Population layout with index of locations/nodes
    h2_cavern_file : str
        Path to CSV file containing hydrogen cavern storage potentials
    cavern_types : list
        List of underground storage types to consider
    clustered_gas_network_file : str, optional
        Path to CSV file containing gas network data
    gas_input_nodes : pd.DataFrame
        DataFrame containing gas input node information (LNG, pipeline, etc.)
    spatial : object, optional
        Object containing spatial information about nodes and their locations
    options : dict, optional
        Dictionary of configuration options. Defaults to empty dict if not provided.
        Key options include:
        - hydrogen_fuel_cell : bool
        - hydrogen_turbine : bool
        - hydrogen_underground_storage : bool
        - gas_network : bool
        - H2_retrofit : bool
        - H2_network : bool
        - methanation : bool
        - coal_cc : bool
        - SMR_cc : bool
        - SMR : bool
        - min_part_load_methanation : float
        - cc_fraction : float
    logger : logging.Logger, optional
        Logger for output messages. If None, no logging is performed.

    Returns
    -------
    None
        The function modifies the network object in-place by adding various
        storage and grid components.

    Notes
    -----
    This function adds multiple types of storage and grid infrastructure:
    - Hydrogen infrastructure (electrolysis, fuel cells, storage)
    - Gas network infrastructure
    - Battery storage systems
    - Carbon capture and conversion facilities (if enabled in options)
    """
    # Set defaults
    options = options or {}

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
        capital_cost=costs.at["electrolysis", "capital_cost"],
        p_min_pu=options["min_part_load_electrolysis"],
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
            capital_cost=costs.at["fuel cell", "capital_cost"]
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
            capital_cost=costs.at["OCGT", "capital_cost"]
            * costs.at["OCGT", "efficiency"],  # NB: fixed cost is per MWel
            marginal_cost=costs.at["OCGT", "VOM"],
            lifetime=costs.at["OCGT", "lifetime"],
        )

    h2_caverns = pd.read_csv(h2_cavern_file, index_col=0)

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

        h2_capital_cost = costs.at["hydrogen storage underground", "capital_cost"]

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
    tech = "hydrogen storage tank type 1 including compressor"
    nodes_overground = h2_caverns.index.symmetric_difference(nodes)

    n.add(
        "Store",
        nodes_overground + " H2 Store",
        bus=nodes_overground + " H2",
        e_nom_extendable=True,
        e_cyclic=True,
        carrier="H2 Store",
        capital_cost=costs.at[tech, "capital_cost"],
        lifetime=costs.at[tech, "lifetime"],
    )

    if options["H2_retrofit"]:
        gas_pipes = pd.read_csv(clustered_gas_network_file, index_col=0)

    if options["gas_network"]:
        logger.info(
            "Add natural gas infrastructure, incl. LNG terminals, production, storage and entry-points."
        )
        gas_pipes = pd.read_csv(clustered_gas_network_file, index_col=0)

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
                gas_pipes.length * costs.at["CH4 (g) pipeline", "capital_cost"]
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

        # check if network is already fully connected and only add new pipelines if not
        if len(complement_edges) > 0:
            complement_edges["length"] = complement_edges.apply(
                haversine, axis=1, args=(n,)
            )

            # apply k_edge_augmentation weighted by length of complement edges
            k_edge = options["gas_network_connectivity_upgrade"]
            if augmentation := list(
                k_edge_augmentation(G, k_edge, avail=complement_edges.values)
            ):
                new_gas_pipes = pd.DataFrame(augmentation, columns=["bus0", "bus1"])
                new_gas_pipes["length"] = new_gas_pipes.apply(
                    haversine, axis=1, args=(n,)
                )

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
                    * costs.at["CH4 (g) pipeline", "capital_cost"],
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
            capital_cost=costs.at["H2 (g) pipeline repurposed", "capital_cost"]
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
        h2_buses_loc = n.buses.query("carrier == 'H2'").location  # noqa: F841
        h2_pipes = h2_pipes.query("bus0 in @h2_buses_loc and bus1 in @h2_buses_loc")

        # TODO Add efficiency losses
        n.add(
            "Link",
            h2_pipes.index,
            bus0=h2_pipes.bus0.values + " H2",
            bus1=h2_pipes.bus1.values + " H2",
            p_min_pu=-1,
            p_nom_extendable=True,
            length=h2_pipes.length.values,
            capital_cost=costs.at["H2 (g) pipeline", "capital_cost"]
            * h2_pipes.length.values,
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
        capital_cost=costs.at["battery storage", "capital_cost"],
        lifetime=costs.at["battery storage", "lifetime"],
    )

    n.add(
        "Link",
        nodes + " battery charger",
        bus0=nodes,
        bus1=nodes + " battery",
        carrier="battery charger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        capital_cost=costs.at["battery inverter", "capital_cost"],
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
            capital_cost=costs.at["methanation", "capital_cost"]
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
            capital_cost=costs.at["coal", "efficiency"]
            * costs.at["coal", "capital_cost"]
            + costs.at["biomass CHP capture", "capital_cost"]
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
            capital_cost=costs.at["SMR CC", "capital_cost"],
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
            capital_cost=costs.at["SMR", "capital_cost"],
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
    n: pypsa.Network,
    avail_profile: pd.DataFrame,
    dsm_profile: pd.DataFrame,
    p_set: pd.Series,
    electric_share: pd.Series,
    number_cars: pd.Series,
    temperature: pd.DataFrame,
    spatial: SimpleNamespace,
    options: dict,
) -> None:
    """
    Add electric vehicle (EV) infrastructure to the network.

    Creates EV batteries, chargers, and optional vehicle-to-grid (V2G) components
    with temperature-dependent efficiency and demand-side management capabilities.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object to be modified
    avail_profile : pd.DataFrame
        Availability profile for EV charging with snapshots as index and nodes as columns
    dsm_profile : pd.DataFrame
        Demand-side management profile defining minimum state of charge
        with snapshots as index and nodes as columns
    p_set : pd.Series
        Base power demand profile for EVs
    electric_share : pd.Series
        Share of electric vehicles per node
    number_cars : pd.Series
        Number of cars per node
    temperature : pd.DataFrame
        Ambient temperature per node and timestamp
    spatial : SimpleNamespace
        Spatial configuration containing at least:
        - nodes: list or Index of node names
    options : dict
        Configuration options containing at least:
        - transport_electric_efficiency: float
        - transport_heating_deadband_lower: float
        - transport_heating_deadband_upper: float
        - EV_lower_degree_factor: float
        - EV_upper_degree_factor: float
        - bev_charge_rate: float
        - bev_charge_efficiency: float
        - bev_dsm: bool
        - bev_energy: float
        - bev_dsm_availability: float
        - v2g: bool

    Returns
    -------
    None
        Modifies the network object in-place by adding EV components

    Notes
    -----
    Components added to the network:
    - EV battery buses for each node
    - EV loads with temperature-corrected efficiency
    - BEV chargers with availability profiles
    - Optional EV battery storage if DSM is enabled
    - Optional V2G links if V2G is enabled

    The function accounts for temperature effects on efficiency and implements
    a rolling average smoothing for the power profile.
    """
    # Add EV battery carrier and buses
    n.add("Carrier", "EV battery")

    n.add(
        "Bus",
        spatial.nodes,
        suffix=" EV battery",
        location=spatial.nodes,
        carrier="EV battery",
        unit="MWh_el",
    )

    # Calculate temperature-corrected efficiency
    car_efficiency = options["transport_electric_efficiency"]
    efficiency = get_temp_efficency(
        car_efficiency,
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["EV_lower_degree_factor"],
        options["EV_upper_degree_factor"],
    )

    # Apply rolling average smoothing to power profile
    p_shifted = (p_set + cycling_shift(p_set, 1) + cycling_shift(p_set, 2)) / 3
    cyclic_eff = p_set.div(p_shifted)
    efficiency *= cyclic_eff

    # Calculate load profile
    profile = electric_share * p_set.div(efficiency)

    # Add EV load
    n.add(
        "Load",
        spatial.nodes,
        suffix=" land transport EV",
        bus=spatial.nodes + " EV battery",
        carrier="land transport EV",
        p_set=profile.loc[n.snapshots],
    )

    # Add BEV chargers
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

    # Add demand-side management components if enabled
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

        # Add vehicle-to-grid if enabled
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


def add_fuel_cell_cars(
    n: pypsa.Network,
    p_set: pd.Series,
    fuel_cell_share: float,
    temperature: pd.Series,
    options: dict,
    spatial: SimpleNamespace,
) -> None:
    """
    Add hydrogen fuel cell vehicles to the network as hydrogen loads.

    Creates temperature-dependent hydrogen demand profiles for fuel cell vehicles
    based on transport energy demand, fuel cell share, and temperature-dependent
    efficiency factors.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object to be modified
    p_set : pd.Series
        Base transport energy demand profile
    fuel_cell_share : float
        Share of transport demand met by fuel cell vehicles (between 0 and 1)
    temperature : pd.Series
        Temperature time series used for efficiency correction
    options : dict
        Configuration options containing at least:
        - transport_fuel_cell_efficiency: float
          Base efficiency of fuel cell vehicles
        - transport_heating_deadband_lower: float
          Lower temperature threshold for efficiency correction
        - transport_heating_deadband_upper: float
          Upper temperature threshold for efficiency correction
        - ICE_lower_degree_factor: float
          Efficiency correction factor for low temperatures
        - ICE_upper_degree_factor: float
          Efficiency correction factor for high temperatures
    spatial : SimpleNamespace
        Spatial configuration containing at least:
        - nodes: list of network nodes
        - h2.nodes: list of hydrogen bus locations

    Returns
    -------
    None
        Modifies the network object in-place by adding fuel cell vehicle loads

    Notes
    -----
    The hydrogen demand is calculated by:
    1. Applying temperature-dependent efficiency corrections
    2. Converting transport energy demand to hydrogen demand
    3. Scaling by the fuel cell vehicle share
    """
    car_efficiency = options["transport_fuel_cell_efficiency"]

    # Calculate temperature-corrected efficiency
    efficiency = get_temp_efficency(
        car_efficiency,
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["ICE_lower_degree_factor"],
        options["ICE_upper_degree_factor"],
    )

    # Calculate hydrogen demand profile
    profile = fuel_cell_share * p_set.div(efficiency)

    # Add hydrogen load for fuel cell vehicles
    n.add(
        "Load",
        spatial.nodes,
        suffix=" land transport fuel cell",
        bus=spatial.h2.nodes,
        carrier="land transport fuel cell",
        p_set=profile.loc[n.snapshots],
    )


def add_ice_cars(
    n: pypsa.Network,
    costs: pd.DataFrame,
    p_set: pd.DataFrame,
    ice_share: pd.DataFrame,
    temperature: pd.DataFrame,
    cf_industry: pd.DataFrame,
    spatial: SimpleNamespace,
    options: dict,
) -> None:
    """
    Add internal combustion engine (ICE) vehicles to the network.

    Creates the necessary infrastructure for representing ICE vehicles in the
    transport sector, including oil buses, temperature-dependent efficiency,
    and CO2 emissions. Can model demand either regionally or aggregated at EU level.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object to be modified
    costs : pd.DataFrame
        Technology cost assumptions with technologies as index and cost parameters
        as columns, must include 'oil' with 'CO2 intensity'
    p_set : pd.DataFrame
        Transport demand time series
    ice_share : pd.DataFrame
        Share of internal combustion engines in transport demand
    temperature : pd.DataFrame
        Hourly temperature time series per node
    cf_industry : pd.DataFrame
        Industrial capacity factors for oil demand
    spatial : SimpleNamespace
        Spatial resolution configuration containing at least:
        - oil.land_transport: names for transport nodes
        - oil.demand_locations: locations of demand
        - oil.nodes: names of oil supply nodes
    options : dict
        Configuration options containing at least:
        - transport_ice_efficiency: float for baseline ICE efficiency
        - transport_heating_deadband_lower: float for lower temperature threshold
        - transport_heating_deadband_upper: float for upper temperature threshold
        - ICE_lower_degree_factor: float for low temperature efficiency impact
        - ICE_upper_degree_factor: float for high temperature efficiency impact
        - regional_oil_demand: bool for regional vs EU-wide demand modeling

    Returns
    -------
    None
        Modifies the network object in-place by adding ICE-related components

    Notes
    -----
    The function adds:
    - Oil carrier buses
    - Temperature-dependent transport oil demand
    - Links from oil supply to transport with CO2 emissions
    """
    add_carrier_buses(
        n=n,
        carrier="oil",
        costs=costs,
        spatial=spatial,
        options=options,
        cf_industry=cf_industry,
    )

    car_efficiency = options["transport_ice_efficiency"]

    # Calculate temperature-corrected efficiency
    efficiency = get_temp_efficency(
        car_efficiency,
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["ICE_lower_degree_factor"],
        options["ICE_upper_degree_factor"],
    )

    # Calculate oil demand profile
    profile = ice_share * p_set.div(efficiency).rename(
        columns=lambda x: x + " land transport oil"
    )

    if not options["regional_oil_demand"]:
        profile = profile.sum(axis=1).to_frame(name="EU land transport oil")

    # Add transport oil buses
    n.add(
        "Bus",
        spatial.oil.land_transport,
        location=spatial.oil.demand_locations,
        carrier="land transport oil",
        unit="land transport",
    )

    # Add transport oil demand
    n.add(
        "Load",
        spatial.oil.land_transport,
        bus=spatial.oil.land_transport,
        carrier="land transport oil",
        p_set=profile.loc[n.snapshots],
    )

    # Add oil supply links with CO2 emissions
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


def add_land_transport(
    n,
    costs,
    transport_demand_file,
    transport_data_file,
    avail_profile_file,
    dsm_profile_file,
    temp_air_total_file,
    cf_industry,
    options,
    investment_year,
    nodes,
) -> None:
    """
    Add land transport demand and infrastructure to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        Cost assumptions for different technologies
    transport_demand_file : str
        Path to CSV file containing transport demand in driven km [100 km]
    transport_data_file : str
        Path to CSV file containing number of cars per region
    avail_profile_file : str
        Path to CSV file containing availability profiles
    dsm_profile_file : str
        Path to CSV file containing demand-side management profiles
    temp_air_total_file : str
        Path to netCDF file containing air temperature data
    options : dict
        Dictionary containing configuration options, specifically:
        - land_transport_fuel_cell_share
        - land_transport_electric_share
        - land_transport_ice_share
    investment_year : int
        Year for which to get the transport shares
    nodes : list-like
        List of spatial nodes to consider

    Returns
    -------
    None
        Modifies the network object in-place by adding transport-related
        components and their properties.

    Notes
    -----
    The function adds different types of land transport (electric vehicles,
    fuel cell vehicles, and internal combustion engines) to the network
    based on specified shares and profiles.
    """
    if logger:
        logger.info("Add land transport")

    # read in transport demand in units driven km [100 km]
    transport = pd.read_csv(transport_demand_file, index_col=0, parse_dates=True)
    number_cars = pd.read_csv(transport_data_file, index_col=0)["number cars"]
    avail_profile = pd.read_csv(avail_profile_file, index_col=0, parse_dates=True)
    dsm_profile = pd.read_csv(dsm_profile_file, index_col=0, parse_dates=True)

    # exogenous share of passenger car type
    engine_types = ["fuel_cell", "electric", "ice"]
    shares = pd.Series()
    for engine in engine_types:
        share_key = f"land_transport_{engine}_share"
        shares[engine] = get(options[share_key], investment_year)
        if logger:
            logger.info(f"{engine} share: {shares[engine] * 100}%")

    check_land_transport_shares(shares)

    p_set = transport[nodes]

    # temperature for correction factor for heating/cooling
    temperature = xr.open_dataarray(temp_air_total_file).to_pandas()

    if shares["electric"] > 0:
        add_EVs(
            n,
            avail_profile,
            dsm_profile,
            p_set,
            shares["electric"],
            number_cars,
            temperature,
            spatial,
            options,
        )

    if shares["fuel_cell"] > 0:
        add_fuel_cell_cars(
            n=n,
            p_set=p_set,
            fuel_cell_share=shares["fuel_cell"],
            temperature=temperature,
            options=options,
            spatial=spatial,
        )
    if shares["ice"] > 0:
        add_ice_cars(
            n,
            costs,
            p_set,
            shares["ice"],
            temperature,
            cf_industry,
            spatial,
            options,
        )


def build_heat_demand(
    n, hourly_heat_demand_file, pop_weighted_energy_totals, heating_efficiencies
):
    """
    Build heat demand time series and adjust electricity load to account for electric heating.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    hourly_heat_demand_file : str
        Path to netCDF file containing hourly heat demand data
    pop_weighted_energy_totals : pd.DataFrame
        Population-weighted energy totals containing columns for total and
        electricity consumption for different sectors and uses
    heating_efficiencies : dict
        Dictionary mapping sector and use combinations to their heating efficiencies

    Returns
    -------
    pd.DataFrame
        Heat demand time series with hierarchical columns for different sectors
        and uses (residential/services, water/space)

    Notes
    -----
    The function:
    - Constructs heat demand profiles for different sectors and uses
    - Adjusts the electricity load profiles by subtracting electric heating
    - Modifies the network object in-place by updating n.loads_t.p_set
    """
    heat_demand_shape = (
        xr.open_dataset(hourly_heat_demand_file).to_dataframe().unstack(level=1)
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
    cop_profiles_file: str,
    direct_heat_source_utilisation_profile_file: str,
    hourly_heat_demand_total_file: str,
    ptes_e_max_pu_file: str,
    ptes_direct_utilisation_profile: str,
    ates_e_nom_max: str,
    ates_capex_as_fraction_of_geothermal_heat_source: float,
    ates_recovery_factor: float,
    enable_ates: bool,
    ates_marginal_cost_charger: float,
    district_heat_share_file: str,
    solar_thermal_total_file: str,
    retro_cost_file: str,
    floor_area_file: str,
    heat_source_profile_files: dict[str, str],
    heat_dsm_profile_file: str,
    params: dict,
    pop_weighted_energy_totals: pd.DataFrame,
    heating_efficiencies: pd.DataFrame,
    pop_layout: pd.DataFrame,
    spatial: object,
    options: dict,
    investment_year: int,
):
    """
    Add heat sector to the network including heat demand, heat pumps, storage, and conversion technologies.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object
    costs : pd.DataFrame
        DataFrame containing cost information for different technologies
    cop_profiles_file : str
        Path to NetCDF file containing coefficient of performance (COP) values for heat pumps
    direct_heat_source_utilisation_profile_file : str
        Path to NetCDF file containing direct heat source utilisation profiles
    hourly_heat_demand_total_file : str
        Path to CSV file containing hourly heat demand data
    ptes_supplemental_heating_required_file: str
        Path to CSV file indicating when supplemental heating for thermal energy storage (TES) is needed
    district_heat_share_file : str
        Path to CSV file containing district heating share information
    solar_thermal_total_file : str
        Path to NetCDF file containing solar thermal generation data
    retro_cost_file : str
        Path to CSV file containing retrofitting costs
    floor_area_file : str
        Path to CSV file containing floor area data
    heat_source_profile_files : dict[str, str]
        Dictionary mapping heat source names to their data file paths
    heat_dsm_profile_file : str
        Path to CSV file containing demand-side management profiles for heat
    params : dict
        Dictionary containing parameters including:
        - heat_pump_sources
        - heat_utilisation_potentials
        - direct_utilisation_heat_sources
    pop_weighted_energy_totals : pd.DataFrame
        Population-weighted energy totals by region
    heating_efficiencies : pd.DataFrame
        Heating system efficiencies
    pop_layout : pd.DataFrame
        Population layout data with columns for fraction and country
    spatial : object
        Object containing spatial data with attributes for different carriers (gas, co2, etc.)
    options : dict
        Dictionary containing configuration options for heat sector components
    investment_year : int
        Year for which to get the heat sector components and costs

    Returns
    -------
    None
        Modifies the network object in-place by adding heat sector components

    Notes
    -----
    The function adds various heat sector components to the network including:
    - Heat demand for different sectors (residential, services)
    - Heat pumps with different heat sources
    - Thermal energy storage if enabled
    - Gas boilers if enabled
    - Solar thermal if enabled
    - Combined heat and power (CHP) plants if enabled
    - Building retrofitting options if enabled
    """
    logger.info("Add heat sector")

    sectors = [sector.value for sector in HeatSector]

    heat_demand = build_heat_demand(
        n,
        hourly_heat_demand_total_file,
        pop_weighted_energy_totals,
        heating_efficiencies,
    )

    cop = xr.open_dataarray(cop_profiles_file)
    direct_heat_profile = xr.open_dataarray(direct_heat_source_utilisation_profile_file)
    district_heat_info = pd.read_csv(district_heat_share_file, index_col=0)
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
            xr.open_dataarray(solar_thermal_total_file)
            .to_pandas()
            .reindex(index=n.snapshots)
        )
        # 1e3 converts from W/m^2 to MW/(1000m^2) = kW/m^2
        solar_thermal = options["solar_cf_correction"] * solar_thermal / 1e3

    for heat_system in (
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

        # if heat_system == HeatSystem.URBAN_CENTRAL and options["central_heat_vent"]:
        if options["heat_vent"][heat_system.system_type.value]:
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
                marginal_cost=-params["sector"]["marginal_cost_heat_vent"],
            )

        ## Add heat load
        factor = heat_system.heat_demand_weighting(
            urban_fraction=urban_fraction[nodes], dist_fraction=dist_fraction[nodes]
        )
        if heat_system != HeatSystem.URBAN_CENTRAL:
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

        else:
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

        if options["residential_heat"]["dsm"]["enable"] and heat_system in [
            HeatSystem.RESIDENTIAL_RURAL,
            HeatSystem.RESIDENTIAL_URBAN_DECENTRAL,
            HeatSystem.URBAN_CENTRAL,
        ]:
            factor = heat_system.heat_demand_weighting(
                urban_fraction=urban_fraction[nodes], dist_fraction=dist_fraction[nodes]
            )

            heat_dsm_profile = pd.read_csv(
                heat_dsm_profile_file,
                header=1,
                index_col=0,
                parse_dates=True,
            )[nodes].reindex(n.snapshots)

            e_nom = (
                heat_demand[["residential space"]]
                .T.groupby(level=1)
                .sum()
                .T[nodes]
                .multiply(factor)
            )

            heat_dsm_restriction_value = options["residential_heat"]["dsm"][
                "restriction_value"
            ].get(investment_year)
            heat_dsm_profile = heat_dsm_profile * heat_dsm_restriction_value
            e_nom = e_nom.max()

            # Allow to overshoot or undercool the target temperatures / heat demand in dsm
            e_min_pu, e_max_pu = 0, 0
            if "overheat" in options["residential_heat"]["dsm"]["direction"]:
                e_max_pu = heat_dsm_profile
            if "undercool" in options["residential_heat"]["dsm"]["direction"]:
                e_min_pu = (-1) * heat_dsm_profile

            # Thermal (standing) losses of buildings assumed to be the same as decentralized water tanks
            n.add(
                "Store",
                nodes,
                suffix=f" {heat_system} heat dsm",
                bus=nodes + f" {heat_system} heat",
                carrier=f"{heat_system} heat dsm",
                standing_loss=costs.at[
                    "decentral water tank storage", "standing_losses"
                ]
                / 100,  # convert %/hour into unit/hour
                e_cyclic=True,
                e_nom=e_nom,
                e_max_pu=e_max_pu,
                e_min_pu=e_min_pu,
            )

            logger.info(f"Adding DSM in {heat_system} heating.")

        if options["tes"]:
            n.add("Carrier", f"{heat_system} water tanks")

            n.add(
                "Bus",
                nodes + f" {heat_system} water tanks",
                location=nodes,
                carrier=f"{heat_system} water tanks",
                unit="MWh_th",
            )

            energy_to_power_ratio_water_tanks = costs.at[
                heat_system.central_or_decentral + " water tank storage",
                "energy to power ratio",
            ]

            n.add(
                "Link",
                nodes,
                suffix=f" {heat_system} water tanks charger",
                bus0=nodes + f" {heat_system} heat",
                bus1=nodes + f" {heat_system} water tanks",
                efficiency=costs.at[
                    heat_system.central_or_decentral + " water tank charger",
                    "efficiency",
                ],
                carrier=f"{heat_system} water tanks charger",
                p_nom_extendable=True,
                marginal_cost=costs.at["water tank charger", "marginal_cost"],
                lifetime=costs.at[
                    heat_system.central_or_decentral + " water tank storage", "lifetime"
                ],
            )

            n.add(
                "Link",
                nodes,
                suffix=f" {heat_system} water tanks discharger",
                bus0=nodes + f" {heat_system} water tanks",
                bus1=nodes + f" {heat_system} heat",
                carrier=f"{heat_system} water tanks discharger",
                efficiency=costs.at[
                    heat_system.central_or_decentral + " water tank discharger",
                    "efficiency",
                ],
                p_nom_extendable=True,
                lifetime=costs.at[
                    heat_system.central_or_decentral + " water tank storage", "lifetime"
                ],
            )

            n.links.loc[
                nodes + f" {heat_system} water tanks charger", "energy to power ratio"
            ] = energy_to_power_ratio_water_tanks

            n.add(
                "Store",
                nodes,
                suffix=f" {heat_system} water tanks",
                bus=nodes + f" {heat_system} water tanks",
                e_cyclic=True,
                e_nom_extendable=True,
                carrier=f"{heat_system} water tanks",
                standing_loss=costs.at[
                    heat_system.central_or_decentral + " water tank storage",
                    "standing_losses",
                ]
                / 100,  # convert %/hour into unit/hour
                capital_cost=costs.at[
                    heat_system.central_or_decentral + " water tank storage",
                    "capital_cost",
                ],
                lifetime=costs.at[
                    heat_system.central_or_decentral + " water tank storage", "lifetime"
                ],
            )

            if heat_system == HeatSystem.URBAN_CENTRAL:
                n.add("Carrier", f"{heat_system} water pits")

                n.add(
                    "Bus",
                    nodes + f" {heat_system} water pits",
                    location=nodes,
                    carrier=f"{heat_system} water pits",
                    unit="MWh_th",
                )

                energy_to_power_ratio_water_pit = costs.at[
                    "central water pit storage", "energy to power ratio"
                ]

                n.add(
                    "Link",
                    nodes,
                    suffix=f" {heat_system} water pits charger",
                    bus0=nodes + f" {heat_system} heat",
                    bus1=nodes + f" {heat_system} water pits",
                    efficiency=costs.at[
                        "central water pit charger",
                        "efficiency",
                    ],
                    carrier=f"{heat_system} water pits charger",
                    p_nom_extendable=True,
                    lifetime=costs.at["central water pit storage", "lifetime"],
                    marginal_cost=costs.at[
                        "central water pit charger", "marginal_cost"
                    ],
                )

                if options["district_heating"]["ptes"]["supplemental_heating"][
                    "enable"
                ]:
                    ptes_supplemental_heating_required = (
                        xr.open_dataarray(ptes_direct_utilisation_profile)
                        .sel(name=nodes)
                        .to_pandas()
                        .reindex(index=n.snapshots)
                    )
                else:
                    ptes_supplemental_heating_required = 1

                n.add(
                    "Link",
                    nodes,
                    suffix=f" {heat_system} water pits discharger",
                    bus0=nodes + f" {heat_system} water pits",
                    bus1=nodes + f" {heat_system} heat",
                    carrier=f"{heat_system} water pits discharger",
                    efficiency=costs.at[
                        "central water pit discharger",
                        "efficiency",
                    ]
                    * ptes_supplemental_heating_required,
                    p_nom_extendable=True,
                    lifetime=costs.at["central water pit storage", "lifetime"],
                )
                n.links.loc[
                    nodes + f" {heat_system} water pits charger",
                    "energy to power ratio",
                ] = energy_to_power_ratio_water_pit

                if options["district_heating"]["ptes"]["dynamic_capacity"]:
                    # Load pre-calculated e_max_pu profiles
                    e_max_pu_data = xr.open_dataarray(ptes_e_max_pu_file)
                    e_max_pu = (
                        e_max_pu_data.sel(name=nodes)
                        .to_pandas()
                        .reindex(index=n.snapshots)
                    )
                else:
                    e_max_pu = 1

                n.add(
                    "Store",
                    nodes,
                    suffix=f" {heat_system} water pits",
                    bus=nodes + f" {heat_system} water pits",
                    e_cyclic=True,
                    e_nom_extendable=True,
                    e_max_pu=e_max_pu,
                    carrier=f"{heat_system} water pits",
                    standing_loss=costs.at[
                        "central water pit storage", "standing_losses"
                    ]
                    / 100,  # convert %/hour into unit/hour
                    capital_cost=costs.at["central water pit storage", "capital_cost"],
                    lifetime=costs.at["central water pit storage", "lifetime"],
                )

        if enable_ates and heat_system == HeatSystem.URBAN_CENTRAL:
            n.add("Carrier", f"{heat_system} aquifer thermal energy storage")

            n.add(
                "Bus",
                nodes + f" {heat_system} aquifer thermal energy storage",
                location=nodes,
                carrier=f"{heat_system} aquifer thermal energy storage",
                unit="MWh_th",
            )

            n.add(
                "Link",
                nodes + f" {heat_system} aquifer thermal energy storage charger",
                bus0=nodes + f" {heat_system} heat",
                bus1=nodes + f" {heat_system} aquifer thermal energy storage",
                efficiency=1.0,
                carrier=f"{heat_system} aquifer thermal energy storage charger",
                p_nom_extendable=True,
                lifetime=costs.at["central geothermal heat source", "lifetime"],
                marginal_cost=ates_marginal_cost_charger,
                capital_cost=costs.at["central geothermal heat source", "capital_cost"]
                * ates_capex_as_fraction_of_geothermal_heat_source
                / 2,
            )

            n.add(
                "Link",
                nodes + f" {heat_system} aquifer thermal energy storage discharger",
                bus1=nodes + f" {heat_system} heat",
                bus0=nodes + f" {heat_system} aquifer thermal energy storage",
                efficiency=1.0,
                carrier=f"{heat_system} aquifer thermal energy storage discharger",
                p_nom_extendable=True,
                lifetime=costs.at["central geothermal heat source", "lifetime"],
                capital_cost=costs.at["central geothermal heat source", "capital_cost"]
                * ates_capex_as_fraction_of_geothermal_heat_source
                / 2,
            )

            ates_e_nom_max = pd.read_csv(ates_e_nom_max, index_col=0)["ates_potential"]
            n.add(
                "Store",
                nodes,
                suffix=f" {heat_system} aquifer thermal energy storage",
                bus=nodes + f" {heat_system} aquifer thermal energy storage",
                e_cyclic=True,
                e_nom_extendable=True,
                e_nom_max=ates_e_nom_max[nodes],
                carrier=f"{heat_system} aquifer thermal energy storage",
                standing_loss=1 - ates_recovery_factor ** (1 / 8760),
                lifetime=costs.at["central geothermal heat source", "lifetime"],
            )

        ## Add heat pumps
        for heat_source in params.heat_pump_sources[heat_system.system_type.value]:
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
                else costs.loc[[costs_name_heat_pump], ["efficiency"]]
            )

            if heat_source in params.limited_heat_sources:
                # get potential
                p_max_source = pd.read_csv(
                    heat_source_profile_files[heat_source],
                    index_col=0,
                    parse_dates=True,
                ).squeeze()[nodes]

                # if only dimension is nodes, convert series to dataframe with columns as nodes and index as snapshots
                if p_max_source.ndim == 1:
                    p_max_source = pd.DataFrame(
                        [p_max_source] * len(n.snapshots),
                        index=n.snapshots,
                        columns=nodes,
                    )

                # add resource
                heat_carrier = f"{heat_system} {heat_source} heat"
                n.add("Carrier", heat_carrier)
                n.add(
                    "Bus",
                    nodes,
                    location=nodes,
                    suffix=f" {heat_carrier}",
                    carrier=heat_carrier,
                )

                # TODO: implement better handling of zero-cost heat sources
                try:
                    capital_cost = (
                        costs.at[
                            heat_system.heat_source_costs_name(heat_source),
                            "capital_cost",
                        ]
                        * overdim_factor
                    )
                    lifetime = costs.at[
                        heat_system.heat_source_costs_name(heat_source), "lifetime"
                    ]
                except KeyError:
                    logger.warning(
                        f"Heat source {heat_source} not found in cost data. Assuming zero cost and infinite lifetime."
                    )
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
                    p_nom_max=p_max_source.max(),
                    p_max_pu=p_max_source / p_max_source.max(),
                )
                # add heat pump converting source heat + electricity to urban central heat
                n.add(
                    "Link",
                    nodes,
                    suffix=f" {heat_system} {heat_source} heat pump",
                    bus0=nodes + f" {heat_system} heat",
                    bus1=nodes,
                    bus2=nodes + f" {heat_carrier}",
                    carrier=f"{heat_system} {heat_source} heat pump",
                    efficiency=(1 / cop_heat_pump.clip(lower=0.001)).squeeze(),
                    efficiency2=(1 - (1 / cop_heat_pump.clip(lower=0.001))).squeeze(),
                    capital_cost=costs.at[costs_name_heat_pump, "capital_cost"]
                    * overdim_factor,
                    p_nom_extendable=True,
                    p_min_pu=(
                        -cop_heat_pump / cop_heat_pump.clip(lower=0.001)
                    ).squeeze(),
                    p_max_pu=0,
                    lifetime=costs.at[costs_name_heat_pump, "lifetime"],
                )

                if heat_source in params.direct_utilisation_heat_sources:
                    # 1 if source temperature exceeds forward temperature, 0 otherwise:
                    efficiency_direct_utilisation = (
                        direct_heat_profile.sel(
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

            if (
                not options["district_heating"]["ptes"]["supplemental_heating"][
                    "enable"
                ]
                and options["district_heating"]["ptes"]["supplemental_heating"][
                    "booster_heat_pump"
                ]
            ):
                raise ValueError(
                    "'booster_heat_pump' is true, but 'enable' is false in 'supplemental_heating'."
                )

            if (
                heat_source in params.temperature_limited_stores
                and options["district_heating"]["ptes"]["supplemental_heating"][
                    "enable"
                ]
                and options["district_heating"]["ptes"]["supplemental_heating"][
                    "booster_heat_pump"
                ]
            ):
                n.add(
                    "Link",
                    nodes,
                    suffix=f" {heat_system} {heat_source} heat pump",
                    bus0=nodes + f" {heat_system} heat",
                    bus1=nodes,
                    bus2=nodes + f" {heat_system} water pits",
                    carrier=f"{heat_system} {heat_source} heat pump",
                    efficiency=(1 / (cop_heat_pump - 1).clip(lower=0.001)).squeeze(),
                    efficiency2=(1 - 1 / cop_heat_pump.clip(lower=0.001)).squeeze(),
                    capital_cost=costs.at[costs_name_heat_pump, "capital_cost"]
                    * overdim_factor,
                    p_nom_extendable=True,
                    p_min_pu=(
                        -cop_heat_pump / cop_heat_pump.clip(lower=0.001)
                    ).squeeze(),
                    p_max_pu=0,
                    lifetime=costs.at[costs_name_heat_pump, "lifetime"],
                )

            else:
                n.add(
                    "Link",
                    nodes,
                    suffix=f" {heat_system} {heat_source} heat pump",
                    bus0=nodes + f" {heat_system} heat",
                    bus1=nodes,
                    carrier=f"{heat_system} {heat_source} heat pump",
                    efficiency=(1 / cop_heat_pump.clip(lower=0.001)).squeeze(),
                    capital_cost=costs.at[costs_name_heat_pump, "capital_cost"]
                    * overdim_factor,
                    p_min_pu=(
                        -cop_heat_pump / cop_heat_pump.clip(lower=0.001)
                    ).squeeze(),
                    p_max_pu=0,
                    p_nom_extendable=True,
                    lifetime=costs.at[costs_name_heat_pump, "lifetime"],
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
                * costs.at[key, "capital_cost"]
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
                * costs.at[key, "capital_cost"]
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
                    heat_system.central_or_decentral + " solar thermal", "capital_cost"
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
                    carrier=f"urban central {fuel} CHP",
                    p_nom_extendable=True,
                    capital_cost=costs.at["central gas CHP", "capital_cost"]
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
                    carrier=f"urban central {fuel} CHP CC",
                    p_nom_extendable=True,
                    capital_cost=costs.at["central gas CHP", "capital_cost"]
                    * costs.at["central gas CHP", "efficiency"]
                    + costs.at["biomass CHP capture", "capital_cost"]
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
                capital_cost=costs.at["micro CHP", "capital_cost"],
                lifetime=costs.at["micro CHP", "lifetime"],
            )

    if options["retrofitting"]["retro_endogen"]:
        logger.info("Add retrofitting endogenously")

        # retrofitting data 'retro_data' with 'costs' [EUR/m^2] and heat
        # demand 'dE' [per unit of original heat demand] for each country and
        # different retrofitting strengths [additional insulation thickness in m]
        retro_data = pd.read_csv(
            retro_cost_file,
            index_col=[0, 1],
            skipinitialspace=True,
            header=[0, 1],
        )
        # heated floor area [10^6 * m^2] per country
        floor_area_file = pd.read_csv(floor_area_file, index_col=[0, 1])

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
            if "services" in name:
                sec = "services"
            elif "residential" in name:
                sec = "residential"
            elif "urban central" in name:
                sec = "tot"
            else:
                raise ValueError(f"Unknown sector in {name}")

            # get floor aread at node and region (urban/rural) in m^2
            floor_area_node = (
                pop_layout.loc[node].fraction * floor_area_file.loc[ct, "value"] * 10**6
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


def add_methanol(
    n: pypsa.Network,
    costs: pd.DataFrame,
    options: dict,
    spatial: SimpleNamespace,
    pop_layout: pd.DataFrame,
) -> None:
    """
    Add methanol-related components to the network.

    Adds methanol infrastructure including production, conversion, and power
    generation facilities based on specified options. Components can include
    biomass-to-methanol plants (with and without carbon capture), methanol
    power plants, and methanol reforming facilities.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object to be modified
    costs : pd.DataFrame
        Technology cost assumptions with technologies as index and cost parameters
        as columns
    options : dict
        Configuration options containing at least:
        - methanol: dict with boolean flags for different methanol technologies
        - biomass: bool indicating if biomass technologies are enabled
    spatial : SimpleNamespace
        Spatial resolution and location-specific parameters for component placement
    pop_layout : pd.DataFrame
        Population data per node used for methanol power plant placement

    Returns
    -------
    None
        Modifies the network object in-place by adding methanol-related components

    Notes
    -----
    The function checks the following methanol options:
    - biomass_to_methanol: Enables biomass to methanol conversion
    - biomass_to_methanol_cc: Enables biomass to methanol with carbon capture
    - methanol_to_power: Enables methanol power plants
    - methanol_reforming: Enables methanol reforming
    - methanol_reforming_cc: Enables methanol reforming with carbon capture
    """
    methanol_options = options["methanol"]
    if not any(
        v if isinstance(v, bool) else any(v.values()) for v in methanol_options.values()
    ):
        return

    logger.info("Add methanol")
    add_carrier_buses(
        n=n,
        carrier="methanol",
        costs=costs,
        spatial=spatial,
        options=options,
    )

    if options["biomass"]:
        if methanol_options["biomass_to_methanol"]:
            add_biomass_to_methanol(n=n, costs=costs)

        if methanol_options["biomass_to_methanol_cc"]:
            add_biomass_to_methanol_cc(n=n, costs=costs)

    if methanol_options["methanol_to_power"]:
        add_methanol_to_power(
            n=n,
            costs=costs,
            pop_layout=pop_layout,
            types=methanol_options["methanol_to_power"],
        )

    if methanol_options["methanol_reforming"]:
        add_methanol_reforming(n=n, costs=costs)

    if methanol_options["methanol_reforming_cc"]:
        add_methanol_reforming_cc(n=n, costs=costs)


def add_biomass(
    n,
    costs,
    options,
    spatial,
    cf_industry,
    pop_layout,
    biomass_potentials_file,
    biomass_transport_costs_file=None,
    nyears=1,
):
    """
    Add biomass-related components to the PyPSA network.

    This function adds various biomass-related components including biogas,
    solid biomass, municipal solid waste, biomass transport, and different
    biomass conversion technologies (CHP, boilers, BtL, BioSNG, etc.).

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        DataFrame containing technology cost assumptions
    options : dict
        Dictionary of configuration options including keys like:
        - gas_network : bool
        - biomass_transport : bool
        - biomass_spatial : bool
        - municipal_solid_waste : bool
        - biomass_to_liquid : bool
        - etc.
    spatial : object
        Object containing spatial information about different carriers (gas, biomass, etc.)
    cf_industry : dict
        Dictionary containing industrial sector configuration
    pop_layout : pd.DataFrame
        DataFrame containing population layout information
    biomass_potentials_file : str
        Path to CSV file containing biomass potentials data
    biomass_transport_costs_file : str, optional
        Path to CSV file containing biomass transport costs data.
        Required if biomass_transport or biomass_spatial options are True.
    nyears : float
        Number of years for which to scale the biomass potentials.

    Returns
    -------
    None
        The function modifies the network object in-place by adding
        biomass-related components.

    Notes
    -----
    The function adds various types of biomass-related components depending
    on the options provided, including:
    - Biogas and solid biomass generators
    - Municipal solid waste if enabled
    - Biomass transport infrastructure
    - Biomass conversion technologies (CHP, boilers, BtL, BioSNG)
    - Carbon capture options for different processes
    """
    logger.info("Add biomass")

    biomass_potentials = pd.read_csv(biomass_potentials_file, index_col=0) * nyears

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
        biomass_import_max_amount = (
            options["solid_biomass_import"]["max_amount"] * 1e6 * nyears
        )
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

        add_carrier_buses(
            n,
            carrier="oil",
            costs=costs,
            spatial=spatial,
            options=options,
            cf_industry=cf_industry,
        )

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

    if options["biogas_upgrading"]:
        n.add(
            "Link",
            spatial.gas.biogas_to_gas,
            bus0=spatial.gas.biogas,
            bus1=spatial.gas.nodes,
            bus2="co2 atmosphere",
            carrier="biogas to gas",
            capital_cost=costs.at["biogas", "capital_cost"]
            + costs.at["biogas upgrading", "capital_cost"],
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
            capital_cost=costs.at["biogas CC", "capital_cost"]
            + costs.at["biogas upgrading", "capital_cost"]
            + costs.at["biomass CHP capture", "capital_cost"]
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
        transport_costs = pd.read_csv(biomass_transport_costs_file, index_col=0)
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
        transport_costs = pd.read_csv(biomass_transport_costs_file, index_col=0)
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
            capital_cost=costs.at[key, "capital_cost"] * costs.at[key, "efficiency"],
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
            capital_cost=costs.at[key + " CC", "capital_cost"]
            * costs.at[key + " CC", "efficiency"]
            + costs.at["biomass CHP capture", "capital_cost"]
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
                * costs.at["biomass boiler", "capital_cost"]
                * options["overdimension_heat_generators"][
                    HeatSystem(name).central_or_decentral
                ],
                marginal_cost=costs.at["biomass boiler", "pelletizing cost"],
                lifetime=costs.at["biomass boiler", "lifetime"],
            )

    # Solid biomass to liquid fuel
    if options["biomass_to_liquid"]:
        add_carrier_buses(
            n,
            carrier="oil",
            costs=costs,
            spatial=spatial,
            options=options,
            cf_industry=cf_industry,
        )
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
            capital_cost=costs.at["BtL", "capital_cost"]
            * costs.at["BtL", "efficiency"],
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
            capital_cost=costs.at["BtL", "capital_cost"] * costs.at["BtL", "efficiency"]
            + costs.at["biomass CHP capture", "capital_cost"]
            * costs.at["BtL", "CO2 stored"],
            marginal_cost=costs.at["BtL", "VOM"] * costs.at["BtL", "efficiency"],
        )

    # Electrobiofuels (BtL with hydrogen addition to make more use of biogenic carbon).
    # Combination of efuels and biomass to liquid, both based on Fischer-Tropsch.
    # Experimental version - use with caution
    if options["electrobiofuels"]:
        add_carrier_buses(
            n,
            carrier="oil",
            costs=costs,
            spatial=spatial,
            options=options,
            cf_industry=cf_industry,
        )
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
            efficiency2=-costs.at["electrobiofuels", "efficiency-biomass"]
            / costs.at["electrobiofuels", "efficiency-hydrogen"],
            efficiency3=-costs.at["solid biomass", "CO2 intensity"]
            + costs.at["BtL", "CO2 stored"]
            * (1 - costs.at["Fischer-Tropsch", "capture rate"]),
            p_nom_extendable=True,
            capital_cost=costs.at["BtL", "capital_cost"] * costs.at["BtL", "efficiency"]
            + efuel_scale_factor
            * costs.at["Fischer-Tropsch", "capital_cost"]
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
            capital_cost=costs.at["BioSNG", "capital_cost"]
            * costs.at["BioSNG", "efficiency"],
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
            capital_cost=costs.at["BioSNG", "capital_cost"]
            * costs.at["BioSNG", "efficiency"]
            + costs.at["biomass CHP capture", "capital_cost"]
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
            capital_cost=costs.at["solid biomass to hydrogen", "capital_cost"]
            * costs.at["solid biomass to hydrogen", "efficiency"]
            + costs.at["biomass CHP capture", "capital_cost"]
            * costs.at["solid biomass", "CO2 intensity"],
            marginal_cost=0.0,
            lifetime=25,  # TODO: add value to technology-data
        )


def add_industry(
    n: pypsa.Network,
    costs: pd.DataFrame,
    industrial_demand_file: str,
    pop_layout: pd.DataFrame,
    pop_weighted_energy_totals: pd.DataFrame,
    options: dict,
    spatial: SimpleNamespace,
    cf_industry: dict,
    investment_year: int,
):
    """
    Add industry and their corresponding carrier buses to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        Costs data including carbon capture, fuel costs, etc.
    industrial_demand_file : str
        Path to CSV file containing industrial demand data
    pop_layout : pd.DataFrame
        Population layout data with index of nodes
    pop_weighted_energy_totals : pd.DataFrame
        Population-weighted energy totals including aviation and navigation data
    options : dict
        Dictionary of configuration options including:
        - biomass_spatial
        - biomass_transport
        - gas_network
        - methanol configuration
        - regional_oil_demand
        - shipping shares (hydrogen, methanol, oil)
        - and others
    spatial : object
        Object containing spatial configuration for different carriers
        (biomass, gas, oil, methanol, etc.)
    cf_industry : dict
        Industry-specific configuration parameters
    investment_year : int
        Year for which investment costs should be considered
    HeatSystem : Enum
        Enumeration defining different heat system types

    Returns
    -------
    None
        The function modifies the network object in-place by adding
        industry-related components.

    Notes
    -----
    This function adds multiple components to the network including:
    - Industrial demand for various carriers
    - Shipping and aviation infrastructure
    - Carbon capture facilities
    - Heat systems
    - Process emission handling
    """
    logger.info("Add industrial demand")
    # add oil buses for shipping, aviation and naptha for industry
    add_carrier_buses(
        n,
        carrier="oil",
        costs=costs,
        spatial=spatial,
        options=options,
        cf_industry=cf_industry,
    )
    add_carrier_buses(
        n,
        carrier="methanol",
        costs=costs,
        spatial=spatial,
        options=options,
        cf_industry=cf_industry,
    )

    nodes = pop_layout.index
    nhours = n.snapshot_weightings.generators.sum()
    nyears = nhours / 8760

    # 1e6 to convert TWh to MWh
    industrial_demand = pd.read_csv(industrial_demand_file, index_col=0) * 1e6 * nyears

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
        capital_cost=costs.at["cement capture", "capital_cost"]
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
        capital_cost=costs.at["cement capture", "capital_cost"]
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
        capital_cost=costs.at["methanolisation", "capital_cost"]
        * options["MWh_MeOH_per_MWh_H2"],  # EUR/MW_H2/a
        marginal_cost=options["MWh_MeOH_per_MWh_H2"]
        * costs.at["methanolisation", "VOM"],
        lifetime=costs.at["methanolisation", "lifetime"],
        efficiency=options["MWh_MeOH_per_MWh_H2"],
        efficiency2=-options["MWh_MeOH_per_MWh_H2"] / options["MWh_MeOH_per_MWh_e"],
        efficiency3=-options["MWh_MeOH_per_MWh_H2"] / options["MWh_MeOH_per_tCO2"],
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
                    * costs.at["decentral oil boiler", "capital_cost"]
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
        capital_cost=costs.at["Fischer-Tropsch", "capital_cost"]
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
        logger.warning(f"Changing HVC demand by {demand_factor * 100 - 100:+.2f}%.")

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
    # link to supply the naphtha for industry load
    n.add(
        "Link",
        spatial.oil.naphtha,
        bus0=spatial.oil.nodes,
        bus1=spatial.oil.naphtha,
        bus2=spatial.co2.process_emissions,
        carrier="naphtha for industry",
        p_nom_extendable=True,
        efficiency2=process_co2_per_naphtha,
    )

    non_sequestered = 1 - get(
        cf_industry["HVC_environment_sequestration_fraction"],
        investment_year,
    )
    # energetic efficiency from naphtha to HVC
    HVC_per_naphtha = (
        costs.at["oil", "CO2 intensity"] - process_co2_per_naphtha
    ) / costs.at["oil", "CO2 intensity"]

    # distribute HVC waste across population
    if len(spatial.oil.non_sequestered_hvc) == 1:
        HVC_potential = p_set_naphtha.sum() * nhours * non_sequestered * HVC_per_naphtha
    else:
        HVC_potential_sum = (
            p_set_naphtha.sum() * nhours * non_sequestered * HVC_per_naphtha
        )
        shares = pop_layout.total / pop_layout.total.sum()
        HVC_potential = shares.mul(HVC_potential_sum)
        HVC_potential.index = HVC_potential.index + " non-sequestered HVC"

    n.add("Carrier", "non-sequestered HVC")

    n.add(
        "Bus",
        spatial.oil.non_sequestered_hvc,
        location=spatial.oil.demand_locations,
        carrier="non-sequestered HVC",
        unit="MWh_LHV",
    )
    # add stores with population distributed potential - must be zero at the last step
    e_max_pu = pd.DataFrame(
        1, index=n.snapshots, columns=spatial.oil.non_sequestered_hvc
    )
    e_max_pu.iloc[-1, :] = 0

    n.add(
        "Store",
        spatial.oil.non_sequestered_hvc,
        bus=spatial.oil.non_sequestered_hvc,
        carrier="non-sequestered HVC",
        e_nom=HVC_potential,
        marginal_cost=0,
        e_initial=HVC_potential,
        e_max_pu=e_max_pu,
    )

    n.add(
        "Link",
        spatial.oil.demand_locations,
        suffix=" HVC to air",
        bus0=spatial.oil.non_sequestered_hvc,
        bus1="co2 atmosphere",
        carrier="HVC to air",
        p_nom_extendable=True,
        efficiency=costs.at["oil", "CO2 intensity"],
    )

    if cf_industry["waste_to_energy"] or cf_industry["waste_to_energy_cc"]:
        if options["biomass"] and options["municipal_solid_waste"]:
            n.add(
                "Link",
                spatial.msw.locations,
                bus0=spatial.msw.nodes,
                bus1=spatial.oil.non_sequestered_hvc,
                bus2="co2 atmosphere",
                carrier="municipal solid waste",
                p_nom_extendable=True,
                efficiency=1.0,
                efficiency2=-costs.at[
                    "oil", "CO2 intensity"
                ],  # because msw is co2 neutral and will be burned in waste CHP or decomposed as oil
            )

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
                bus0=spatial.oil.non_sequestered_hvc,
                bus1=spatial.nodes,
                bus2=urban_central_nodes,
                bus3="co2 atmosphere",
                carrier="waste CHP",
                p_nom_extendable=True,
                capital_cost=costs.at["waste CHP", "capital_cost"]
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
                bus0=spatial.oil.non_sequestered_hvc,
                bus1=spatial.nodes,
                bus2=urban_central_nodes,
                bus3="co2 atmosphere",
                bus4=spatial.co2.nodes,
                carrier="waste CHP CC",
                p_nom_extendable=True,
                capital_cost=costs.at["waste CHP CC", "capital_cost"]
                * costs.at["waste CHP CC", "efficiency"]
                + costs.at["biomass CHP capture", "capital_cost"]
                * costs.at["oil", "CO2 intensity"],
                marginal_cost=costs.at["waste CHP CC", "VOM"],
                efficiency=costs.at["waste CHP CC", "efficiency"],
                efficiency2=costs.at["waste CHP CC", "efficiency-heat"],
                efficiency3=costs.at["oil", "CO2 intensity"]
                * (1 - options["cc_fraction"]),
                efficiency4=costs.at["oil", "CO2 intensity"] * options["cc_fraction"],
                lifetime=costs.at["waste CHP CC", "lifetime"],
            )

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
        capital_cost=costs.at["cement capture", "capital_cost"],
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
        add_carrier_buses(
            n,
            carrier="coal",
            costs=costs,
            spatial=spatial,
            options=options,
            cf_industry=cf_industry,
        )

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


def add_aviation(
    n: pypsa.Network,
    costs: pd.DataFrame,
    pop_layout: pd.DataFrame,
    pop_weighted_energy_totals: pd.DataFrame,
    options: dict,
    spatial: SimpleNamespace,
) -> None:
    logger.info("Add aviation")

    nodes = pop_layout.index
    nhours = n.snapshot_weightings.generators.sum()

    demand_factor = options["aviation_demand_factor"]
    if demand_factor != 1:
        logger.warning(
            f"Changing aviation demand by {demand_factor * 100 - 100:+.2f}%."
        )

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
        tech = "methanol-to-kerosene"

        logger.info(f"Adding {tech}.")

        capital_cost = costs.at[tech, "capital_cost"] / costs.at[tech, "methanol-input"]

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
            efficiency3=costs.at["oil", "CO2 intensity"]
            / costs.at[tech, "methanol-input"],
            p_nom_extendable=True,
            lifetime=costs.at[tech, "lifetime"],
        )


def add_shipping(
    n: pypsa.Network,
    costs: pd.DataFrame,
    shipping_demand_file: str,
    pop_layout: pd.DataFrame,
    pop_weighted_energy_totals: pd.DataFrame,
    options: dict,
    spatial: SimpleNamespace,
    investment_year: int,
) -> None:
    logger.info("Add shipping")

    nodes = pop_layout.index
    nhours = n.snapshot_weightings.generators.sum()
    nyears = nhours / 8760

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
        pd.read_csv(shipping_demand_file, index_col=0).squeeze(axis=1) * nyears
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
                capital_cost=costs.at["H2 liquefaction", "capital_cost"],
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


def add_waste_heat(
    n: pypsa.Network,
    costs: pd.DataFrame,
    options: dict,
    cf_industry: dict,
) -> None:
    """
    Add industrial waste heat utilization capabilities to district heating systems.

    Modifies the network by adding waste heat outputs from various industrial processes
    to urban central heating systems. The amount of waste heat that can be utilized
    is controlled by efficiency parameters and option flags.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        DataFrame containing technology cost and efficiency parameters,
        particularly for methanolisation process
    options : dict
        Configuration dictionary containing boolean flags for different waste heat sources:
        - use_fischer_tropsch_waste_heat
        - use_methanation_waste_heat
        - use_haber_bosch_waste_heat
        - use_methanolisation_waste_heat
        - use_electrolysis_waste_heat
        - use_fuel_cell_waste_heat
    cf_industry : dict
        Dictionary containing conversion factors for industrial processes, including:
        - MWh_H2_per_tNH3_electrolysis
        - MWh_elec_per_tNH3_electrolysis
        - MWh_NH3_per_tNH3

    Returns
    -------
    None
        Modifies the network object in-place by adding waste heat connections

    Notes
    -----
    - Waste heat is only added to buses with carrier "urban central heat"
    - Default efficiency values (like 0.95 for Fischer-Tropsch) might need
      to be moved to configuration
    - The modification adds additional output buses (bus2, bus3, or bus4) to
      existing links representing industrial processes
    """
    logger.info("Add possibility to use industrial waste heat in district heating")

    # AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if not urban_central.empty:
        urban_central = urban_central.str[: -len(" urban central heat")]

        link_carriers = n.links.carrier.unique()

        # Fischer-Tropsch waste heat
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

        # Sabatier process waste heat
        if options["use_methanation_waste_heat"] and "Sabatier" in link_carriers:
            n.links.loc[urban_central + " Sabatier", "bus3"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " Sabatier", "efficiency3"] = (
                0.95 - n.links.loc[urban_central + " Sabatier", "efficiency"]
            ) * options["use_methanation_waste_heat"]

        # Haber-Bosch process waste heat
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

        # Methanolisation waste heat
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

        # Electrolysis waste heat
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

        # Fuel cell waste heat
        if options["use_fuel_cell_waste_heat"] and "H2 Fuel Cell" in link_carriers:
            n.links.loc[urban_central + " H2 Fuel Cell", "bus2"] = (
                urban_central + " urban central heat"
            )
            n.links.loc[urban_central + " H2 Fuel Cell", "efficiency2"] = (
                0.95 - n.links.loc[urban_central + " H2 Fuel Cell", "efficiency"]
            ) * options["use_fuel_cell_waste_heat"]


def add_agriculture(
    n: pypsa.Network,
    costs: pd.DataFrame,
    pop_layout: pd.DataFrame,
    pop_weighted_energy_totals: pd.DataFrame,
    investment_year: int,
    options: dict,
    spatial: SimpleNamespace,
) -> None:
    """
    Add agriculture, forestry and fishing sector loads to the network.

    Creates electrical loads, heat loads, and machinery-related components (both electric
    and oil-based) for the agricultural sector, distributed according to population weights.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object
    costs : pd.DataFrame
        DataFrame containing cost parameters, must include 'oil' index with 'CO2 intensity'
    pop_layout : pd.DataFrame
        Population distribution by node
    pop_weighted_energy_totals : pd.DataFrame
        Energy totals weighted by population, must include columns:
        ['total agriculture electricity', 'total agriculture heat',
         'total agriculture machinery']
    investment_year : int
        Year for which to get time-dependent parameters
    options : dict
        Configuration dictionary containing:
        - agriculture_machinery_electric_share: float or dict
        - agriculture_machinery_oil_share: float or dict
        - agriculture_machinery_fuel_efficiency: float
        - agriculture_machinery_electric_efficiency: float
        - regional_oil_demand: bool
    spatial : SimpleNamespace
        Namespace containing spatial definitions, must include:
        - oil.agriculture_machinery
        - oil.demand_locations
        - oil.nodes

    Returns
    -------
    None
        Modifies the network object in-place by adding loads and components

    Notes
    -----
    The function adds three types of loads:
    1. Agricultural electricity consumption
    2. Agricultural heat demand
    3. Agricultural machinery energy demand (split between electricity and oil)
    """
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
        """
        Define how attributes should be clustered.
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

        m = n.copy(snapshots=[])
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
    """Split bidirectional links into two unidirectional links to include transmission losses."""

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

    n.links["reversed"] = n.links.get("reversed", False)
    n.links = pd.concat([n.links, rev_links], sort=False)
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


def add_enhanced_geothermal(
    n,
    costs,
    costs_config,
    egs_potentials,
    egs_overlap,
    egs_config,
    egs_capacity_factors=None,
):
    """
    Add Enhanced Geothermal System (EGS) potential to the network model.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network container object.
    egs_potentials : str
        Path to CSV file containing EGS potential data.
    egs_overlap : str
        Path to CSV file defining overlap between gridded geothermal potential
        estimation and bus regions.
    costs : pd.DataFrame
        Technology cost assumptions including fields for lifetime, FOM, investment,
        and efficiency parameters.
    egs_config : dict
        Configuration for enhanced geothermal systems with keys:
        - var_cf : bool
            Whether to use time-varying capacity factors
        - flexible : bool
            Whether to add flexible operation using geothermal reservoir
        - max_hours : float
            Maximum hours of storage if flexible
        - max_boost : float
            Maximum power boost factor if flexible
    costs_config : dict
        General cost configuration containing:
        - fill_values : dict
            With key 'discount rate' for financial calculations
    egs_capacity_factors : str, optional
        Path to CSV file with time-varying capacity factors.
        Required if egs_config['var_cf'] is True.

    Returns
    -------
    None
        Modifies the network object in-place by adding EGS components.

    Notes
    -----
    Implements EGS with 2020 CAPEX from Aghahosseini et al 2021.
    The function adds various components to the network:
    - Geothermal heat generators and buses
    - Organic Rankine Cycle links for electricity generation
    - District heating links if urban central heat exists
    - Optional storage units for flexible operation
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

    assert (egs_potentials["capital_cost"] > 0).all(), (
        "Error in EGS cost, negative values found."
    )

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
        efficiency = pd.read_csv(egs_capacity_factors, parse_dates=True, index_col=0)
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
            efficiency=bus_eta.loc[n.snapshots],
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


def get_capacities_from_elec(n, carriers, component):
    """
    Gets capacities and efficiencies for {carrier} in n.{component} that were
    previously assigned in add_electricity.
    """
    component_list = ["generators", "storage_units", "links", "stores"]
    component_dict = {name: getattr(n, name) for name in component_list}
    e_nom_carriers = ["stores"]
    nom_col = {x: "e_nom" if x in e_nom_carriers else "p_nom" for x in component_list}
    eff_col = "efficiency"

    capacity_dict = {}
    efficiency_dict = {}
    for carrier in carriers:
        capacity_dict[carrier] = component_dict[component].query("carrier in @carrier")[
            nom_col[component]
        ]
        efficiency_dict[carrier] = component_dict[component].query(
            "carrier in @carrier"
        )[eff_col]

    return capacity_dict, efficiency_dict


def add_import_options(
    n: pypsa.Network,
    costs: pd.DataFrame,
    options: dict,
    gas_input_nodes: pd.DataFrame,
):
    """
    Add green energy import options.

    Parameters
    ----------
    n : pypsa.Network
    costs : pd.DataFrame
    options : dict
        Options from snakemake.params["sector"].
    gas_input_nodes : pd.DataFrame
        Locations of gas input nodes split by LNG and pipeline.
    """

    import_config = options["imports"]
    import_options = import_config["price"]
    logger.info(f"Adding import options:\n{pd.Series(import_options)}")

    if "methanol" in import_options:
        co2_intensity = costs.at["methanolisation", "carbondioxide-input"]

        n.add(
            "Link",
            spatial.methanol.nodes,
            suffix=" import",
            carrier="import methanol",
            bus0="co2 atmosphere",
            bus1=spatial.methanol.nodes,
            efficiency=1 / co2_intensity,
            marginal_cost=import_options["methanol"] / co2_intensity,
            p_nom=1e7,
        )

    if "oil" in import_options:
        co2_intensity = costs.at["oil", "CO2 intensity"]

        n.add(
            "Link",
            spatial.oil.nodes,
            suffix=" import",
            carrier="import oil",
            bus0="co2 atmosphere",
            bus1=spatial.oil.nodes,
            efficiency=1 / co2_intensity,
            marginal_cost=import_options["oil"] / co2_intensity,
            p_nom=1e7,
        )

    if "gas" in import_options:
        co2_intensity = costs.at["gas", "CO2 intensity"]

        p_nom = gas_input_nodes["lng"].dropna()
        p_nom.rename(lambda x: x + " gas", inplace=True)  #
        if len(spatial.gas.nodes) == 1:
            p_nom = p_nom.sum()
            nodes = spatial.gas.nodes
        else:
            nodes = p_nom.index

        n.add(
            "Link",
            nodes,
            suffix=" import",
            carrier="import gas",
            bus0="co2 atmosphere",
            bus1=nodes,
            efficiency=1 / co2_intensity,
            marginal_cost=import_options["gas"] / co2_intensity,
            p_nom=p_nom / co2_intensity,
        )

    if "NH3" in import_options:
        if options["ammonia"]:
            n.add(
                "Generator",
                spatial.ammonia.nodes,
                suffix=" import",
                bus=spatial.ammonia.nodes,
                carrier="import NH3",
                p_nom=1e7,
                marginal_cost=import_options["NH3"],
            )

        else:
            logger.warning(
                "Skipping specified ammonia imports because ammonia carrier is not present."
            )

    if "H2" in import_options:
        p_nom = gas_input_nodes["pipeline"].dropna()
        p_nom.rename(lambda x: x + " H2", inplace=True)

        n.add(
            "Generator",
            p_nom.index,
            suffix=" import",
            bus=p_nom.index,
            carrier="import H2",
            p_nom=p_nom,
            marginal_cost=import_options["H2"],
        )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_sector_network",
            opts="",
            clusters="10",
            sector_opts="",
            planning_horizons="2050",
        )

    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    options = snakemake.params.sector
    cf_industry = snakemake.params.industry

    investment_year = int(snakemake.wildcards.planning_horizons)

    n = pypsa.Network(snakemake.input.network)

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    nhours = n.snapshot_weightings.generators.sum()
    nyears = nhours / 8760

    costs = load_costs(snakemake.input.costs)

    pop_weighted_energy_totals = (
        pd.read_csv(snakemake.input.pop_weighted_energy_totals, index_col=0) * nyears
    )
    pop_weighted_heat_totals = (
        pd.read_csv(snakemake.input.pop_weighted_heat_totals, index_col=0) * nyears
    )
    pop_weighted_energy_totals.update(pop_weighted_heat_totals)

    fn = snakemake.input.gas_input_nodes_simplified
    gas_input_nodes = pd.read_csv(fn, index_col=0)

    if options.get("keep_existing_capacities", False):
        existing_capacities, existing_efficiencies = get_capacities_from_elec(
            n,
            carriers=options.get("conventional_generation").keys(),
            component="generators",
        )
    else:
        existing_capacities = existing_efficiencies = None

    carriers_to_keep = snakemake.params.pypsa_eur
    profiles = {
        key: snakemake.input[key]
        for key in snakemake.input.keys()
        if key.startswith("profile")
    }
    landfall_lengths = {
        tech: settings["landfall_length"]
        for tech, settings in snakemake.params.renewable.items()
        if "landfall_length" in settings.keys()
    }
    patch_electricity_network(n, costs, carriers_to_keep, profiles, landfall_lengths)

    fn = snakemake.input.heating_efficiencies
    year = int(snakemake.params["energy_totals_year"])
    heating_efficiencies = pd.read_csv(fn, index_col=[1, 0]).loc[year]

    spatial = define_spatial(pop_layout.index, options)

    if snakemake.params.foresight in ["overnight", "myopic", "perfect"]:
        add_lifetime_wind_solar(n, costs)

        conventional = snakemake.params.conventional_carriers
        for carrier in conventional:
            add_carrier_buses(
                n=n,
                carrier=carrier,
                costs=costs,
                spatial=spatial,
                options=options,
                cf_industry=cf_industry,
            )

    add_eu_bus(n)

    emission_prices = snakemake.params["emission_prices"]
    co2_price = (
        get(emission_prices["co2"], investment_year)
        if emission_prices["enable"]
        else 0.0
    )
    add_co2_tracking(
        n,
        costs,
        options,
        sequestration_potential_file=snakemake.input.sequestration_potential,
        co2_price=co2_price,
    )

    add_generation(
        n=n,
        costs=costs,
        pop_layout=pop_layout,
        conventionals=options["conventional_generation"],
        spatial=spatial,
        options=options,
        cf_industry=cf_industry,
        ext_carriers=snakemake.params.electricity.get("extendable_carriers", dict()),
        existing_capacities=existing_capacities,
        existing_efficiencies=existing_efficiencies,
    )

    add_storage_and_grids(
        n=n,
        costs=costs,
        pop_layout=pop_layout,
        h2_cavern_file=snakemake.input.h2_cavern,
        cavern_types=snakemake.params.sector["hydrogen_underground_storage_locations"],
        clustered_gas_network_file=snakemake.input.clustered_gas_network,
        gas_input_nodes=gas_input_nodes,
        spatial=spatial,
        options=options,
    )

    if options["transport"]:
        add_land_transport(
            n=n,
            costs=costs,
            transport_demand_file=snakemake.input.transport_demand,
            transport_data_file=snakemake.input.transport_data,
            avail_profile_file=snakemake.input.avail_profile,
            dsm_profile_file=snakemake.input.dsm_profile,
            temp_air_total_file=snakemake.input.temp_air_total,
            cf_industry=cf_industry,
            options=options,
            investment_year=investment_year,
            nodes=spatial.nodes,
        )

    if options["heating"]:
        add_heat(
            n=n,
            costs=costs,
            cop_profiles_file=snakemake.input.cop_profiles,
            direct_heat_source_utilisation_profile_file=snakemake.input.direct_heat_source_utilisation_profiles,
            hourly_heat_demand_total_file=snakemake.input.hourly_heat_demand_total,
            ptes_e_max_pu_file=snakemake.input.ptes_e_max_pu_profiles,
            ates_e_nom_max=snakemake.input.ates_potentials,
            ates_capex_as_fraction_of_geothermal_heat_source=snakemake.params.sector[
                "district_heating"
            ]["ates"]["capex_as_fraction_of_geothermal_heat_source"],
            ates_marginal_cost_charger=snakemake.params.sector["district_heating"][
                "ates"
            ]["marginal_cost_charger"],
            ates_recovery_factor=snakemake.params.sector["district_heating"]["ates"][
                "recovery_factor"
            ],
            enable_ates=snakemake.params.sector["district_heating"]["ates"]["enable"],
            ptes_direct_utilisation_profile=snakemake.input.ptes_direct_utilisation_profiles,
            district_heat_share_file=snakemake.input.district_heat_share,
            solar_thermal_total_file=snakemake.input.solar_thermal_total,
            retro_cost_file=snakemake.input.retro_cost,
            floor_area_file=snakemake.input.floor_area,
            heat_source_profile_files={
                source: snakemake.input[source]
                for source in snakemake.params.limited_heat_sources
                if source in snakemake.input.keys()
            },
            heat_dsm_profile_file=snakemake.input.heat_dsm_profile,
            params=snakemake.params,
            pop_weighted_energy_totals=pop_weighted_energy_totals,
            heating_efficiencies=heating_efficiencies,
            pop_layout=pop_layout,
            spatial=spatial,
            options=options,
            investment_year=investment_year,
        )

    if options["biomass"]:
        add_biomass(
            n=n,
            costs=costs,
            options=options,
            spatial=spatial,
            cf_industry=cf_industry,
            pop_layout=pop_layout,
            biomass_potentials_file=snakemake.input.biomass_potentials,
            biomass_transport_costs_file=snakemake.input.biomass_transport_costs,
            nyears=nyears,
        )

    if options["ammonia"]:
        add_ammonia(n, costs, pop_layout, spatial, cf_industry)

    if options["methanol"]:
        add_methanol(n, costs, options=options, spatial=spatial, pop_layout=pop_layout)

    if options["industry"]:
        add_industry(
            n=n,
            costs=costs,
            industrial_demand_file=snakemake.input.industrial_demand,
            pop_layout=pop_layout,
            pop_weighted_energy_totals=pop_weighted_energy_totals,
            options=options,
            spatial=spatial,
            cf_industry=cf_industry,
            investment_year=investment_year,
        )

    if options["shipping"]:
        add_shipping(
            n=n,
            costs=costs,
            shipping_demand_file=snakemake.input.shipping_demand,
            pop_layout=pop_layout,
            pop_weighted_energy_totals=pop_weighted_energy_totals,
            options=options,
            spatial=spatial,
            investment_year=investment_year,
        )

    if options["aviation"]:
        add_aviation(
            n=n,
            costs=costs,
            pop_layout=pop_layout,
            pop_weighted_energy_totals=pop_weighted_energy_totals,
            options=options,
            spatial=spatial,
        )

    if options["heating"]:
        add_waste_heat(n, costs, options, cf_industry)

    if options["agriculture"]:  # requires H and I
        add_agriculture(
            n,
            costs,
            pop_layout,
            pop_weighted_energy_totals,
            investment_year,
            options,
            spatial,
        )

    if options["dac"]:
        add_dac(n, costs)

    if not options["electricity_transmission_grid"]:
        decentral(n)

    if not options["H2_network"]:
        remove_h2_network(n)

    if options["co2_network"]:
        add_co2_network(
            n,
            costs,
            co2_network_cost_factor=snakemake.config["sector"][
                "co2_network_cost_factor"
            ],
        )

    if options["allam_cycle_gas"]:
        add_allam_gas(n, costs, pop_layout=pop_layout, spatial=spatial)

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
                snakemake.params.countries,
                snakemake.params.planning_horizons,
            )
        co2_cap = pd.read_csv(fn, index_col=0).squeeze()
        limit = co2_cap.loc[investment_year]
    else:
        limit = get(co2_budget, investment_year)
    add_co2limit(
        n,
        options,
        snakemake.input.co2_totals_name,
        snakemake.params.countries,
        nyears,
        limit,
    )

    maxext = snakemake.params["lines"]["max_extension"]
    if maxext is not None:
        limit_individual_line_extension(n, maxext)

    if options["electricity_distribution_grid"]:
        insert_electricity_distribution_grid(
            n, costs, options, pop_layout, snakemake.input.solar_rooftop_potentials
        )

    if options["enhanced_geothermal"].get("enable", False):
        logger.info("Adding Enhanced Geothermal Systems (EGS).")
        add_enhanced_geothermal(
            n,
            costs=costs,
            costs_config=snakemake.config["costs"],
            egs_potentials=snakemake.input["egs_potentials"],
            egs_overlap=snakemake.input["egs_overlap"],
            egs_config=snakemake.params["sector"]["enhanced_geothermal"],
            egs_capacity_factors="path/to/capacity_factors.csv",
        )

    if options["imports"]["enable"]:
        add_import_options(n, costs, options, gas_input_nodes)

    if options["gas_distribution_grid"]:
        insert_gas_distribution_costs(n, costs, options=options)

    if options["electricity_grid_connection"]:
        add_electricity_grid_connection(n, costs)

    for k, v in options["transmission_efficiency"].items():
        if k in options["transmission_efficiency"]["enable"]:
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
