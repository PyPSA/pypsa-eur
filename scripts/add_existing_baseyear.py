# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Adds existing power and heat generation capacities for initial planning
horizon.
"""

import logging
import re
from types import SimpleNamespace

import country_converter as coco
import geopandas as gpd
import numpy as np
import pandas as pd
import geopandas as gpd
import pycountry
import powerplantmatching as pm
import pypsa
import xarray as xr

from scripts._helpers import (
    configure_logging,
    sanitize_custom_columns,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.add_electricity import load_costs, sanitize_carriers
from scripts.build_energy_totals import cartesian
from scripts.definitions.heat_system import HeatSystem
from scripts.prepare_sector_network import cluster_heat_buses, define_spatial

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()
idx = pd.IndexSlice
spatial = SimpleNamespace()


def add_build_year_to_new_assets(n: pypsa.Network, baseyear: int) -> None:
    """
    Add build year to new assets in the network.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    baseyear : int
        Year in which optimized assets are built
    """
    # Give assets with lifetimes and no build year the build year baseyear
    for c in n.iterate_components(["Link", "Generator", "Store"]):
        assets = c.df.index[(c.df.lifetime != np.inf) & (c.df.build_year == 0)]
        c.df.loc[assets, "build_year"] = baseyear

        # add -baseyear to name
        rename = pd.Series(c.df.index, c.df.index)
        rename[assets] += f"-{str(baseyear)}"
        c.df.rename(index=rename, inplace=True)

        # rename time-dependent
        selection = n.component_attrs[c.name].type.str.contains(
            "series"
        ) & n.component_attrs[c.name].status.str.contains("Input")
        for attr in n.component_attrs[c.name].index[selection]:
            c.pnl[attr] = c.pnl[attr].rename(columns=rename)


def add_existing_renewables(
    n: pypsa.Network,
    costs: pd.DataFrame,
    df_agg: pd.DataFrame,
    countries: list[str],
) -> None:
    """
    Add existing renewable capacities to conventional power plant data.

    Parameters
    ----------
    df_agg : pd.DataFrame
        DataFrame containing conventional power plant data
    costs : pd.DataFrame
        Technology cost data with 'lifetime' column indexed by technology
    n : pypsa.Network
        Network containing topology and generator data
    countries : list
        List of country codes to consider

    Returns
    -------
    None
        Modifies df_agg in-place
    """
    tech_map = {"solar": "PV", "onwind": "Onshore", "offwind-ac": "Offshore"}

    irena = pm.data.IRENASTAT().powerplant.convert_country_to_alpha2()
    irena = irena.query("Country in @countries")
    irena = irena.groupby(["Technology", "Country", "Year"]).Capacity.sum()

    irena = irena.unstack().reset_index()

    for carrier, tech in tech_map.items():
        df = (
            irena[irena.Technology.str.contains(tech)]
            .drop(columns=["Technology"])
            .set_index("Country")
        )
        df.columns = df.columns.astype(int)

        # calculate yearly differences
        df.insert(loc=0, value=0.0, column="1999")
        df = df.diff(axis=1).drop("1999", axis=1).clip(lower=0)

        # distribute capacities among generators potential (p_nom_max)
        gen_i = n.generators.query("carrier == @carrier").index
        carrier_gens = n.generators.loc[gen_i]
        res_capacities = []
        for country, group in carrier_gens.groupby(
            carrier_gens.bus.map(n.buses.country)
        ):
            fraction = group.p_nom_max / group.p_nom_max.sum()
            res_capacities.append(cartesian(df.loc[country], fraction))
        res_capacities = pd.concat(res_capacities, axis=1).T

        for year in res_capacities.columns:
            for gen in res_capacities.index:
                bus_bin = re.sub(f" {carrier}.*", "", gen)
                bus, bin_id = bus_bin.rsplit(" ", maxsplit=1)
                name = f"{bus_bin} {carrier}-{year}"
                capacity = res_capacities.loc[gen, year]
                if capacity > 0.0:
                    cost_key = carrier.split("-", maxsplit=1)[0]
                    df_agg.at[name, "Fueltype"] = carrier
                    df_agg.at[name, "Capacity"] = capacity
                    df_agg.at[name, "DateIn"] = year
                    df_agg.at[name, "lifetime"] = costs.at[cost_key, "lifetime"]
                    df_agg.at[name, "DateOut"] = (
                        year + costs.at[cost_key, "lifetime"] - 1
                    )
                    df_agg.at[name, "bus"] = bus
                    df_agg.at[name, "resource_class"] = bin_id

    df_agg["resource_class"] = df_agg["resource_class"].fillna(0)


def add_power_capacities_installed_before_baseyear(
    n: pypsa.Network,
    costs: pd.DataFrame,
    grouping_years: list[int],
    baseyear: int,
    powerplants_file: str,
    countries: list[str],
    capacity_threshold: float,
    lifetime_values: dict[str, float],
) -> None:
    """
    Add power generation capacities installed before base year.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    costs : pd.DataFrame
        Technology costs
    grouping_years : list
        Intervals to group existing capacities
    baseyear : int
        Base year for analysis
    powerplants_file : str
        Path to powerplants CSV file
    countries : list
        List of countries to consider
    capacity_threshold : float
        Minimum capacity threshold
    lifetime_values : dict
        Default values for missing data
    """
    logger.debug(f"Adding power capacities installed before {baseyear}")

    df_agg = pd.read_csv(powerplants_file, index_col=0)

    rename_fuel = {
        "Hard Coal": "coal",
        "Lignite": "lignite",
        "Nuclear": "nuclear",
        "Oil": "oil",
        "OCGT": "OCGT",
        "CCGT": "CCGT",
        "Bioenergy": "urban central solid biomass CHP",
    }

    # Replace Fueltype "Natural Gas" with the respective technology (OCGT or CCGT)
    df_agg.loc[df_agg["Fueltype"] == "Natural Gas", "Fueltype"] = df_agg.loc[
        df_agg["Fueltype"] == "Natural Gas", "Technology"
    ]

    fueltype_to_drop = [
        "Hydro",
        "Wind",
        "Solar",
        "Geothermal",
        "Waste",
        "Other",
        "CCGT, Thermal",
    ]

    technology_to_drop = ["Pv", "Storage Technologies"]

    # drop unused fueltypes and technologies
    df_agg.drop(df_agg.index[df_agg.Fueltype.isin(fueltype_to_drop)], inplace=True)
    df_agg.drop(df_agg.index[df_agg.Technology.isin(technology_to_drop)], inplace=True)
    df_agg.Fueltype = df_agg.Fueltype.map(rename_fuel)

    # Intermediate fix for DateIn & DateOut
    # Fill missing DateIn
    biomass_i = df_agg.loc[df_agg.Fueltype == "urban central solid biomass CHP"].index
    mean = df_agg.loc[biomass_i, "DateIn"].mean()
    df_agg.loc[biomass_i, "DateIn"] = df_agg.loc[biomass_i, "DateIn"].fillna(int(mean))
    # Fill missing DateOut
    dateout = df_agg.loc[biomass_i, "DateIn"] + lifetime_values["lifetime"]
    df_agg.loc[biomass_i, "DateOut"] = df_agg.loc[biomass_i, "DateOut"].fillna(dateout)

    # include renewables in df_agg
    add_existing_renewables(
        df_agg=df_agg,
        costs=costs,
        n=n,
        countries=countries,
    )
    # drop assets which are already phased out / decommissioned
    phased_out = df_agg[df_agg["DateOut"] < baseyear].index
    df_agg.drop(phased_out, inplace=True)

    newer_assets = (df_agg.DateIn > max(grouping_years)).sum()
    if newer_assets:
        logger.warning(
            f"There are {newer_assets} assets with build year "
            f"after last power grouping year {max(grouping_years)}. "
            "These assets are dropped and not considered."
            "Consider to redefine the grouping years to keep them."
        )
        to_drop = df_agg[df_agg.DateIn > max(grouping_years)].index
        df_agg.drop(to_drop, inplace=True)

    df_agg["grouping_year"] = np.take(
        grouping_years, np.digitize(df_agg.DateIn, grouping_years, right=True)
    )

    # calculate (adjusted) remaining lifetime before phase-out (+1 because assuming
    # phase out date at the end of the year)
    df_agg["lifetime"] = df_agg.DateOut - df_agg["grouping_year"] + 1

    df = df_agg.pivot_table(
        index=["grouping_year", "Fueltype", "resource_class"],
        columns="bus",
        values="Capacity",
        aggfunc="sum",
    )

    lifetime = df_agg.pivot_table(
        index=["grouping_year", "Fueltype", "resource_class"],
        columns="bus",
        values="lifetime",
        aggfunc="mean",  # currently taken mean for clustering lifetimes
    )

    carrier = {
        "OCGT": "gas",
        "CCGT": "gas",
        "coal": "coal",
        "oil": "oil",
        "lignite": "lignite",
        "nuclear": "uranium",
        "urban central solid biomass CHP": "biomass",
    }

    for grouping_year, generator, resource_class in df.index:
        # capacity is the capacity in MW at each node for this
        capacity = df.loc[grouping_year, generator, resource_class]
        capacity = capacity[~capacity.isna()]
        capacity = capacity[capacity > capacity_threshold]
        suffix = "-ac" if generator == "offwind" else ""
        name_suffix = f" {generator}{suffix}-{grouping_year}"
        asset_i = capacity.index + name_suffix
        if generator in ["solar", "onwind", "offwind-ac"]:
            asset_i = capacity.index + " " + resource_class + name_suffix
            name_suffix = " " + resource_class + name_suffix
            cost_key = generator.split("-")[0]
            # to consider electricity grid connection costs or a split between
            # solar utility and rooftop as well, rather take cost assumptions
            # from existing network than from the cost database
            capital_cost = n.generators.loc[
                n.generators.carrier == generator + suffix, "capital_cost"
            ].mean()
            marginal_cost = n.generators.loc[
                n.generators.carrier == generator + suffix, "marginal_cost"
            ].mean()
            # check if assets are already in network (e.g. for 2020)
            already_build = n.generators.index.intersection(asset_i)
            new_build = asset_i.difference(n.generators.index)

            # this is for the year 2020
            if not already_build.empty:
                n.generators.loc[already_build, "p_nom"] = n.generators.loc[
                    already_build, "p_nom_min"
                ] = capacity.loc[already_build.str.replace(name_suffix, "")].values
            new_capacity = capacity.loc[new_build.str.replace(name_suffix, "")]

            name_suffix_by = f" {resource_class} {generator}{suffix}-{baseyear}"
            p_max_pu = n.generators_t.p_max_pu[capacity.index + name_suffix_by]

            if not new_build.empty:
                n.add(
                    "Generator",
                    new_capacity.index,
                    suffix=name_suffix,
                    bus=new_capacity.index,
                    carrier=generator,
                    p_nom=new_capacity,
                    marginal_cost=marginal_cost,
                    capital_cost=capital_cost,
                    efficiency=costs.at[cost_key, "efficiency"],
                    p_max_pu=p_max_pu.rename(columns=n.generators.bus),
                    build_year=grouping_year,
                    lifetime=costs.at[cost_key, "lifetime"],
                )

        else:
            bus0 = vars(spatial)[carrier[generator]].nodes
            if "EU" not in vars(spatial)[carrier[generator]].locations:
                bus0 = bus0.intersection(capacity.index + " " + carrier[generator])

            # check for missing bus
            missing_bus = pd.Index(bus0).difference(n.buses.index)
            if not missing_bus.empty:
                logger.info(f"add buses {bus0}")
                n.add(
                    "Bus",
                    bus0,
                    carrier=generator,
                    location=vars(spatial)[carrier[generator]].locations,
                    unit="MWh_el",
                )

            already_build = n.links.index.intersection(asset_i)
            new_build = asset_i.difference(n.links.index)
            lifetime_assets = lifetime.loc[
                grouping_year, generator, resource_class
            ].dropna()

            # this is for the year 2020
            if not already_build.empty:
                n.links.loc[already_build, "p_nom_min"] = capacity.loc[
                    already_build.str.replace(name_suffix, "")
                ].values

            if not new_build.empty:
                new_capacity = capacity.loc[new_build.str.replace(name_suffix, "")]

                if generator != "urban central solid biomass CHP":
                    n.add(
                        "Link",
                        new_capacity.index,
                        suffix=name_suffix,
                        bus0=bus0,
                        bus1=new_capacity.index,
                        bus2="co2 atmosphere",
                        carrier=generator,
                        marginal_cost=costs.at[generator, "efficiency"]
                        * costs.at[generator, "VOM"],  # NB: VOM is per MWel
                        capital_cost=costs.at[generator, "efficiency"]
                        * costs.at[
                            generator, "capital_cost"
                        ],  # NB: fixed cost is per MWel
                        p_nom=new_capacity / costs.at[generator, "efficiency"],
                        efficiency=costs.at[generator, "efficiency"],
                        efficiency2=costs.at[carrier[generator], "CO2 intensity"],
                        build_year=grouping_year,
                        lifetime=lifetime_assets.loc[new_capacity.index],
                    )
                else:
                    key = "central solid biomass CHP"
                    central_heat = n.buses.query(
                        "carrier == 'urban central heat'"
                    ).location.unique()
                    heat_buses = new_capacity.index.map(
                        lambda i: i + " urban central heat" if i in central_heat else ""
                    )

                    n.add(
                        "Link",
                        new_capacity.index,
                        suffix=name_suffix,
                        bus0=spatial.biomass.df.loc[new_capacity.index]["nodes"].values,
                        bus1=new_capacity.index,
                        bus2=heat_buses,
                        carrier=generator,
                        p_nom=new_capacity / costs.at[key, "efficiency"],
                        capital_cost=costs.at[key, "capital_cost"]
                        * costs.at[key, "efficiency"],
                        marginal_cost=costs.at[key, "VOM"],
                        efficiency=costs.at[key, "efficiency"],
                        build_year=grouping_year,
                        efficiency2=costs.at[key, "efficiency-heat"],
                        lifetime=lifetime_assets.loc[new_capacity.index],
                    )
        # check if existing capacities are larger than technical potential
        existing_large = n.generators[
            n.generators["p_nom_min"] > n.generators["p_nom_max"]
        ].index
        if len(existing_large):
            logger.warning(
                f"Existing capacities larger than technical potential for {existing_large},\
                           adjust technical potential to existing capacities"
            )
            n.generators.loc[existing_large, "p_nom_max"] = n.generators.loc[
                existing_large, "p_nom_min"
            ]


def get_efficiency(
    heat_system: HeatSystem,
    carrier: str,
    nodes: pd.Index,
    efficiencies: dict[str, float],
    costs: pd.DataFrame,
) -> pd.Series | float:
    """
    Computes the heating system efficiency based on the sector and carrier
    type.

    Parameters
    ----------
    heat_system : object
    carrier : str
        The type of fuel or energy carrier (e.g., 'gas', 'oil').
    nodes : pandas.Series
        A pandas Series containing node information used to match the heating efficiency data.
    efficiencies : dict
        A dictionary containing efficiency values for different carriers and sectors.
    costs : pandas.DataFrame
        A DataFrame containing boiler cost and efficiency data for different heating systems.

    Returns
    -------
    efficiency : pandas.Series or float
        A pandas Series mapping the efficiencies based on nodes for residential and services sectors, or a single
        efficiency value for other heating systems (e.g., urban central).

    Notes
    -----
    - For residential and services sectors, efficiency is mapped based on the nodes.
    - For other sectors, the default boiler efficiency is retrieved from the `costs` database.
    """

    if heat_system.value == "urban central":
        boiler_costs_name = getattr(heat_system, f"{carrier}_boiler_costs_name")
        efficiency = costs.at[boiler_costs_name, "efficiency"]
    elif heat_system.sector.value == "residential":
        key = f"{carrier} residential space efficiency"
        efficiency = nodes.str[:2].map(efficiencies[key])
    elif heat_system.sector.value == "services":
        key = f"{carrier} services space efficiency"
        efficiency = nodes.str[:2].map(efficiencies[key])
    else:
        raise ValueError(f"Heat system {heat_system} not defined.")

    return efficiency


def add_heating_capacities_installed_before_baseyear(
    n: pypsa.Network,
    costs: pd.DataFrame,
    baseyear: int,
    grouping_years: list[int],
    existing_capacities: pd.DataFrame,
    heat_pump_cop: xr.DataArray,
    heat_pump_source_types: dict[str, list[str]],
    efficiency_file: str,
    use_time_dependent_cop: bool,
    default_lifetime: int,
    energy_totals_year: int,
    capacity_threshold: float,
    use_electricity_distribution_grid: bool,
) -> None:
    """
    Add heating capacities installed before base year.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    costs : pd.DataFrame
        Technology costs
    baseyear : int
        Base year for analysis
    grouping_years : list
        Intervals to group capacities
    heat_pump_cop : xr.DataArray
        Heat pump coefficients of performance
    use_time_dependent_cop : bool
        Use time-dependent COPs
    heating_default_lifetime : int
        Default lifetime for heating systems
    existing_capacities : pd.DataFrame
        Existing heating capacity distribution
    heat_pump_source_types : dict
        Heat pump sources by system type
    efficiency_file : str
        Path to heating efficiencies file
    energy_totals_year : int
        Year for energy totals
    capacity_threshold : float
        Minimum capacity threshold
    use_electricity_distribution_grid : bool
        Whether to use electricity distribution grid
    """
    logger.debug(f"Adding heating capacities installed before {baseyear}")

    # Load heating efficiencies
    heating_efficiencies = pd.read_csv(efficiency_file, index_col=[1, 0]).loc[
        energy_totals_year
    ]

    ratios = []
    valid_grouping_years = []

    for heat_system in existing_capacities.columns.get_level_values(0).unique():
        heat_system = HeatSystem(heat_system)

        nodes = pd.Index(
            n.buses.location[n.buses.index.str.contains(f"{heat_system} heat")]
        )

        if (
            not heat_system == HeatSystem.URBAN_CENTRAL
        ) and use_electricity_distribution_grid:
            nodes_elec = nodes + " low voltage"
        else:
            nodes_elec = nodes

            too_large_grouping_years = [
                gy for gy in grouping_years if gy >= int(baseyear)
            ]
            if too_large_grouping_years:
                logger.warning(
                    f"Grouping years >= baseyear are ignored. Dropping {too_large_grouping_years}."
                )
            valid_grouping_years = pd.Series(
                [
                    int(grouping_year)
                    for grouping_year in grouping_years
                    if int(grouping_year) + default_lifetime > int(baseyear)
                    and int(grouping_year) < int(baseyear)
                ]
            )

            assert valid_grouping_years.is_monotonic_increasing

            if len(valid_grouping_years) == 0:
                logger.warning(
                    f"No valid grouping years found for {heat_system}. "
                    "No existing capacities will be added."
                )
                ratios = []
            else:
                # get number of years of each interval
                _years = valid_grouping_years.diff()
                # Fill NA from .diff() with value for the first interval
                _years[0] = valid_grouping_years[0] - baseyear + default_lifetime
                # Installation is assumed to be linear for the past
                ratios = _years / _years.sum()

        for ratio, grouping_year in zip(ratios, valid_grouping_years):
            # Add heat pumps
            for heat_source in heat_pump_source_types[heat_system.system_type.value]:
                costs_name = heat_system.heat_pump_costs_name(heat_source)

                efficiency = (
                    heat_pump_cop.sel(
                        heat_system=heat_system.system_type.value,
                        heat_source=heat_source,
                        name=nodes,
                    )
                    .to_pandas()
                    .reindex(index=n.snapshots)
                    if use_time_dependent_cop
                    else costs.at[costs_name, "efficiency"]
                )

                n.add(
                    "Link",
                    nodes,
                    suffix=f" {heat_system} {heat_source} heat pump-{grouping_year}",
                    bus0=nodes_elec,
                    bus1=nodes + " " + heat_system.value + " heat",
                    carrier=f"{heat_system} {heat_source} heat pump",
                    efficiency=efficiency,
                    capital_cost=costs.at[costs_name, "efficiency"]
                    * costs.at[costs_name, "capital_cost"],
                    p_nom=existing_capacities.loc[
                        nodes, (heat_system.value, f"{heat_source} heat pump")
                    ]
                    * ratio
                    / costs.at[costs_name, "efficiency"],
                    build_year=int(grouping_year),
                    lifetime=costs.at[costs_name, "lifetime"],
                )

            # add resistive heater, gas boilers and oil boilers
            n.add(
                "Link",
                nodes,
                suffix=f" {heat_system} resistive heater-{grouping_year}",
                bus0=nodes_elec,
                bus1=nodes + " " + heat_system.value + " heat",
                carrier=heat_system.value + " resistive heater",
                efficiency=costs.at[
                    heat_system.resistive_heater_costs_name, "efficiency"
                ],
                capital_cost=(
                    costs.at[heat_system.resistive_heater_costs_name, "efficiency"]
                    * costs.at[heat_system.resistive_heater_costs_name, "capital_cost"]
                ),
                p_nom=(
                    existing_capacities.loc[
                        nodes, (heat_system.value, "resistive heater")
                    ]
                    * ratio
                    / costs.at[heat_system.resistive_heater_costs_name, "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at[heat_system.resistive_heater_costs_name, "lifetime"],
            )

            efficiency = get_efficiency(
                heat_system, "gas", nodes, heating_efficiencies, costs
            )

            n.add(
                "Link",
                nodes,
                suffix=f" {heat_system} gas boiler-{grouping_year}",
                bus0="EU gas" if "EU gas" in spatial.gas.nodes else nodes + " gas",
                bus1=nodes + " " + heat_system.value + " heat",
                bus2="co2 atmosphere",
                carrier=heat_system.value + " gas boiler",
                efficiency=efficiency,
                efficiency2=costs.at["gas", "CO2 intensity"],
                capital_cost=(
                    costs.at[heat_system.gas_boiler_costs_name, "efficiency"]
                    * costs.at[heat_system.gas_boiler_costs_name, "capital_cost"]
                ),
                p_nom=(
                    existing_capacities.loc[nodes, (heat_system.value, "gas boiler")]
                    * ratio
                    / costs.at[heat_system.gas_boiler_costs_name, "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at[heat_system.gas_boiler_costs_name, "lifetime"],
            )

            efficiency = get_efficiency(
                heat_system, "oil", nodes, heating_efficiencies, costs
            )

            n.add(
                "Link",
                nodes,
                suffix=f" {heat_system} oil boiler-{grouping_year}",
                bus0=spatial.oil.nodes,
                bus1=nodes + " " + heat_system.value + " heat",
                bus2="co2 atmosphere",
                carrier=heat_system.value + " oil boiler",
                efficiency=efficiency,
                efficiency2=costs.at["oil", "CO2 intensity"],
                capital_cost=costs.at[heat_system.oil_boiler_costs_name, "efficiency"]
                * costs.at[heat_system.oil_boiler_costs_name, "capital_cost"],
                p_nom=(
                    existing_capacities.loc[nodes, (heat_system.value, "oil boiler")]
                    * ratio
                    / costs.at[heat_system.oil_boiler_costs_name, "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at[
                    f"{heat_system.central_or_decentral} gas boiler", "lifetime"
                ],
            )

            # delete links with p_nom=nan corresponding to extra nodes in country
            n.remove(
                "Link",
                [
                    index
                    for index in n.links.index.to_list()
                    if str(grouping_year) in index and np.isnan(n.links.p_nom[index])
                ],
            )

            # delete links with capacities below threshold
            n.remove(
                "Link",
                [
                    index
                    for index in n.links.index.to_list()
                    if str(grouping_year) in index
                    and n.links.p_nom[index] < capacity_threshold
                ],
            )


def prepare_plant_data(
    regions_fn: str,
    isi_database: str,
) -> tuple[pd.DataFrame, gpd.GeoDataFrame]:
    """
    Reads in the Fraunhofer ISI database with high resolution plant data and maps them to the bus regions.
    Returns the database as df as well as the regions as gdf.

    Parameters
    ----------
    regions_fn : str
        path to the onshore regions file
    isi_database: str
        path to the fraunhofer isi database
    """
    # add existing industry
    regions = gpd.read_file(regions_fn).set_index("name")

    isi_data = pd.read_excel(isi_database, sheet_name="Database", index_col=1)
    # assign bus region to each plant
    geometry = gpd.points_from_xy(isi_data["Longitude"], isi_data["Latitude"])
    plant_data = gpd.GeoDataFrame(isi_data, geometry=geometry, crs="EPSG:4326")
    plant_data = gpd.sjoin(plant_data, regions, how="inner", predicate="within")
    plant_data.rename(columns={"name": "bus"}, inplace=True)
    # filter for countries in model scope
    plant_data = plant_data[plant_data.Country.isin(snakemake.params.countries)]
    # replace UK with GB in Country column
    plant_data["Country"] = plant_data["Country"].replace("UK", "GB")
    # assign industry grouping year
    grouping_years = snakemake.params.existing_capacities["grouping_years_industry"]
    plant_data.loc[:, "Year of last modernisation"] = plant_data[
        "Year of last modernisation"
    ].replace("x", np.nan)
    plant_data["grouping_year"] = 0
    valid_mask = plant_data["Year of last modernisation"].notna()
    valid_years = plant_data.loc[valid_mask, "Year of last modernisation"]
    indices = np.searchsorted(grouping_years, valid_years, side="right")
    plant_data.loc[valid_years.index, "grouping_year"] = np.array(grouping_years)[
        indices
    ]

    return plant_data, regions


def country_to_code(name):
    try:
        return pycountry.countries.lookup(name).alpha_2
    except LookupError:
        return None


def prepare_gem_database():
    """
    Load GEM database of cement plants and map onto bus regions.
    """
    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    df = pd.read_excel(
        snakemake.input.gem_gspt,
        sheet_name="Plant Data",
        na_values=["N/A", "unknown", ">0"],
    )
    
    df["country_code"] = df["Country/Area"].apply(country_to_code)

    df = df[df.country_code.isin(snakemake.params.countries)]

    df["Start date"] = pd.to_numeric(
        df["Start date"].str.split("-").str[0], errors="coerce"
    )

    latlon = (
        df["Coordinates"]
        .str.split(", ", expand=True)
        .rename(columns={0: "lat", 1: "lon"})
    )
    geometry = gpd.points_from_xy(latlon["lon"], latlon["lat"])
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

    gdf = gpd.sjoin(gdf, regions, how="inner", predicate="within")

    gdf.rename(columns={"name": "bus"}, inplace=True)
    gdf["country"] = gdf.bus.str[:2]

    grouping_years = snakemake.params.existing_capacities["grouping_years_industry"]
    gdf = gdf[gdf["Start date"] < grouping_years[-1]]
    avg_age = gdf["Start date"].mean()
    gdf.fillna({"Start date": avg_age}, inplace=True)
    
    gdf["grouping_year"] = 0
    indices = np.searchsorted(grouping_years, gdf["Start date"], side="right")
    gdf.loc[:, "grouping_year"] = np.array(grouping_years)[indices]
    gdf["grouping_year"].replace(2025, 2020, inplace=True)

    return gdf


def prepare_plant_data(
    regions_fn: str,
    isi_database: str,
) -> tuple[pd.DataFrame, gpd.GeoDataFrame]:
    """
    Reads in the Fraunhofer ISI database with high resolution plant data and maps them to the bus regions.
    Returns the database as df as well as the regions as gdf.

    Parameters
    ----------
    regions_fn : str
        path to the onshore regions file
    isi_database: str
        path to the fraunhofer isi database
    """
    # add existing industry
    regions = gpd.read_file(regions_fn).set_index("name")

    isi_data = pd.read_excel(isi_database, sheet_name="Database", index_col=1)
    # assign bus region to each plant
    geometry = gpd.points_from_xy(isi_data["Longitude"], isi_data["Latitude"])
    plant_data = gpd.GeoDataFrame(isi_data, geometry=geometry, crs="EPSG:4326")
    plant_data = gpd.sjoin(plant_data, regions, how="inner", predicate="within")
    plant_data.rename(columns={"name": "bus"}, inplace=True)
    # filter for countries in model scope
    plant_data = plant_data[plant_data.Country.isin(snakemake.params.countries)]
    # replace UK with GB in Country column
    plant_data["Country"] = plant_data["Country"].replace("UK", "GB")
    # assign industry grouping year
    grouping_years = snakemake.params.existing_capacities["grouping_years_industry"]
    plant_data.loc[:, "Year of last modernisation"] = plant_data[
        "Year of last modernisation"
    ].replace("x", np.nan)
    plant_data["grouping_year"] = 0
    valid_mask = plant_data["Year of last modernisation"].notna()
    valid_years = plant_data.loc[valid_mask, "Year of last modernisation"]
    indices = np.searchsorted(grouping_years, valid_years, side="right")
    plant_data.loc[valid_years.index, "grouping_year"] = np.array(grouping_years)[
        indices
    ]

    return plant_data, regions


def add_existing_steel_plants(
    n: pypsa.Network,
) -> None:
    """
    Adds existing steel plants.
    The plants are running on natural gas only since the retrofitting to hydrogen would be associated with costs. Plants are not expected to run at a minimal part load to avoid forcing the use of natural gas in planning horizons with climate targets.
    Exhaust heat is not integrated since assuming that heat is integrated to make the current process more efficient.
    """
    logger.info("Adding existing steel plants.")

    plant_data, regions = prepare_plant_data(
        snakemake.input.regions_onshore,
        snakemake.input.isi_database,
    )
    steel_buses = n.loads[n.loads.carrier=="steel"].bus.to_list()
    steel_buses = [x.replace(" steel", "") for x in steel_buses]

    fh_drg = plant_data[plant_data["Process status qup"] == "Direct reduction NG"]
    drg = fh_drg.groupby(
        ["bus", "Country", "grouping_year", "Product"], as_index=False
    ).agg({
        "Production in tons (calibrated)": "sum",
        "Out": "mean",
        "Year of last modernisation": "mean",
    })

    fh_bof = plant_data[plant_data["Process status qup"] == "Blast furnace"]
    bof = fh_bof.groupby(
        ["bus", "Country", "grouping_year", "Product"], as_index=False
    ).agg({
        "Production in tons (calibrated)": "sum",
        "Out": "mean",
        "Year of last modernisation": "mean",
    })
    # fill Year of last modernisation NaNs with mean
    mean_year = int(bof["Year of last modernisation"].mean())
    bof.loc[:, "Year of last modernisation"] = bof["Year of last modernisation"].fillna(mean_year)
    drg.loc[:, "Year of last modernisation"] = drg["Year of last modernisation"].fillna(2025)
    # fill in lifetime
    bof.loc[:, "Lifetime"] = bof.loc[:, "Out"] - bof.loc[:, "Year of last modernisation"]
    bof.loc[:, "Lifetime"] = bof["Lifetime"].fillna(costs.at["blast furnace-basic oxygen furnace", "lifetime"])
    drg.loc[:, "Lifetime"] = drg.loc[:, "Out"] - drg.loc[:, "Year of last modernisation"]
    drg.loc[:, "Lifetime"] = drg["Lifetime"].fillna(costs.at["natural gas direct iron reduction furnace", "lifetime"])

    drg.index = (
        drg["bus"]
        + " gas DRI-"
        + drg["grouping_year"].astype(str)
    )
    drg = drg[drg.bus.isin(steel_buses)]
    bof.index = (
        bof["bus"]
        + " BOF-"
        + bof["grouping_year"].astype(str)
    )
    bof = bof[bof.bus.isin(steel_buses)]

    # add direct reduction with natural gas
    gas_input = costs.at["natural gas direct iron reduction furnace", "gas-input"]
    marginal_cost = (
            costs.at["iron ore DRI-ready", "commodity"]
            * costs.at["natural gas direct iron reduction furnace", "ore-input"]
            / gas_input
        )

    logger.info(f"Adding {len(drg)} gas DRI plant.")

    n.add(
        "Link",
        drg.index,
        bus0=[bus + " gas" for bus in drg.bus]
        if snakemake.params.sector["gas_network"]
        else "EU gas",
        bus1=[bus + " hbi" for bus in drg.bus] if not snakemake.params.sector["industry_relocation"] else "EU hbi",
        bus2=[bus + " gas DRI emission" for bus in drg.bus],
        p_nom=drg["Production in tons (calibrated)"]
        .mul(gas_input)
        .div(8760)
        .values,
        p_nom_extendable=False,
        carrier="gas DRI",
        efficiency=1/gas_input,
        efficiency2=costs.at["gas", "CO2 intensity"],
        capital_cost=costs.at["natural gas direct iron reduction furnace", "capital_cost"]
            / gas_input,
        marginal_cost=marginal_cost,
        build_year=drg["Year of last modernisation"],
        lifetime=drg["Lifetime"],
    )

    # add coal blast furnaces
    coal_input = costs.at["blast furnace-basic oxygen furnace", "coal-input"]
    marginal_cost = (
        costs.at["iron ore DRI-ready", "commodity"]
        * costs.at["blast furnace-basic oxygen furnace", "ore-input"] / coal_input
    )
    logger.info(f"Adding {len(bof)} BOF plant.")

    n.add(
        "Link",
        bof.index,
        bus0="EU coal",
        bus1=[bus + " steel" for bus in bof.bus],
        bus2=[bus + " BOF emission" for bus in bof.bus],
        p_nom=bof["Production in tons (calibrated)"]
        .mul(coal_input)
        .div(8760)
        .values,
        p_nom_extendable=False,
        carrier="BOF",
        efficiency=1/coal_input,
        efficiency2=costs.at["coal", "CO2 intensity"],
        marginal_cost=marginal_cost,
        capital_cost=costs.at["blast furnace-basic oxygen furnace", "capital_cost"] / coal_input,
        build_year=bof["Year of last modernisation"],
        lifetime=bof["Lifetime"],
    )


def add_existing_cement_plants(n):
    database = prepare_gem_database()
    # get plants that are still operating
    database = database[database["Operating status"] == "operating"]

    cement = database.groupby(
        ["bus", "country", "grouping_year"], as_index=False
    )['Cement Capacity (millions metric tonnes per annum)'].sum()

    cement.index = (cement["bus"]
        + " cement kiln-"
        + cement["grouping_year"].astype(str)
    )
    cement["clinker"] = (cement["bus"]
        + " cement production-"
        + cement["grouping_year"].astype(str)
    )
    # get rid of capacities without demand
    cement_buses = n.loads[n.loads.carrier=="cement"].bus.to_list()
    cement_buses = [x.replace(" cement", "") for x in cement_buses]
    cement = cement[cement.bus.isin(cement_buses)]

    logger.info(f"Adding {len(cement)} existing cement links.")

    # clinker production
    gas_input = costs.at["cement dry clinker", "gas-input"] + costs.at["cement dry clinker", "heat-input"]
    electricity_input = costs.at["cement dry clinker", "electricity-input"]
    co2_emission = 0.79/gas_input + costs.at["gas", "CO2 intensity"]

    n.add(
        "Link",
        cement.index,
        bus0=[bus + " gas" for bus in cement.bus]
        if snakemake.params.sector["gas_network"]
        else "EU gas",
        bus1=[bus + " clinker" for bus in cement.bus],
        bus2=cement.bus,
        bus3=[bus + " cement emission" for bus in cement.bus],
        carrier="cement kiln",
        p_nom_extendable=False,
        p_nom=cement['Cement Capacity (millions metric tonnes per annum)']
        .mul(gas_input)
        .mul(costs.at["cement finishing", "clinker-input"])
        .div(8760)
        .mul(1e6)
        .values,
        p_min_pu=0,
        capital_cost=costs.at["cement dry clinker", "capital_cost"] / gas_input,
        efficiency=1/gas_input,
        efficiency2=-electricity_input,
        efficiency3=co2_emission,
        build_year=cement.grouping_year,
        lifetime=costs.at["cement dry clinker", "lifetime"],
    )

    # cement finishing
    electricity_input = costs.at["cement finishing", "electricity-input"] / costs.at["cement finishing", "clinker-input"]
    clinker_input = costs.at["cement finishing", "clinker-input"]
    n.add(
        "Link",
        cement.clinker.to_list(),
        bus0=[bus + " clinker" for bus in cement.bus],
        bus1=[bus + " cement" for bus in cement.bus],
        bus2=cement.bus.to_list(),
        carrier="cement finishing",
        p_nom_extendable=False,
        p_nom=cement['Cement Capacity (millions metric tonnes per annum)']
        .mul(clinker_input)
        .div(8760)
        .mul(1e6)
        .values,
        p_min_pu=0,
        capital_cost=costs.at["cement finishing", "capital_cost"] / gas_input,
        efficiency=1/clinker_input,
        efficiency2=-electricity_input,
        build_year=cement.grouping_year.to_list(),
        lifetime=costs.at["cement finishing", "lifetime"],
    )


def add_existing_meoh_plants(n):

    logger.info("Adding existing methanol plants.")

    plant_data, regions = prepare_plant_data(
        snakemake.input.regions_onshore,
        snakemake.input.isi_database,
    )

    fh_meoh = plant_data[plant_data.Product == "Methanol"]
    fh_meoh["grouping_year"].replace(2025, 2020, inplace=True)
    fh_meoh = fh_meoh.groupby(
        ["bus", "Country", "grouping_year", "Product"], as_index=False
    )["Production in tons (calibrated)"].sum()

    fh_meoh.index = (
        fh_meoh["bus"]
        + " grey methanol-"
        + fh_meoh["grouping_year"].astype(str)
    )

    # grey methanol efficiency hard coded for now
    costs.at["grey methanol synthesis", "efficiency"] = 1/1.757469244
    capital_cost = (
        costs.at["SMR", "capital_cost"]
        + costs.at["methanolisation", "capital_cost"]
        * costs.at["grey methanol synthesis", "efficiency"]
    )
    co2_emissions = (
        costs.at["gas", "CO2 intensity"]
        - costs.at["grey methanol synthesis", "efficiency"]
        * costs.at["methanol", "CO2 intensity"]
    )
    n.add(
        "Link",
        fh_meoh.index,
        bus0=[bus + " gas" for bus in fh_meoh.bus]
        if snakemake.params.sector["gas_network"]
        else "EU gas",
        bus1="EU methanol",
        bus2="co2 atmosphere",
        p_nom_extendable=False,
        p_nom=fh_meoh["Production in tons (calibrated)"]
        .mul(snakemake.params.MWh_MeOH_per_tMeOH)
        .div(costs.at["grey methanol synthesis", "efficiency"]),
        carrier="grey methanol",
        efficiency=costs.at["grey methanol synthesis", "efficiency"],
        efficiency2=co2_emissions,
        capital_cost=capital_cost,
        build_year=fh_meoh.grouping_year,
        lifetime=costs.at["SMR", "lifetime"],
    )


def add_existing_ammonia_plants(
    n: pypsa.Network,
) -> None:
    """
    Adds existing Haber-Bosch plants.
    The plants are running on natural gas only since the retrofitting to hydrogen would be associated with costs. Plants are not expected to run at a minimal part load to avoid forcing the use of natural gas in planning horizons with climate targets.
    Exhaust heat is not integrated since assuming that heat is integrated to make the current process more efficient.
    """
    logger.info("Adding existing Haber-Bosch plants.")

    plant_data, regions = prepare_plant_data(
        snakemake.input.regions_onshore,
        snakemake.input.isi_database,
    )

    fh_ammonia = plant_data[plant_data.Product == "Ammonia"]

    fh_ammonia = fh_ammonia.groupby(
        ["bus", "Country", "grouping_year", "Product"], as_index=False
    )["Production in tons (calibrated)"].sum()

    fh_ammonia.index = (
        fh_ammonia["bus"]
        + " Haber-Bosch-SMR-"
        + fh_ammonia["grouping_year"].astype(str)
    )
    # add dataset for Non EU27 countries
    df = pd.read_csv(snakemake.input.ammonia, index_col=0)

    geometry = gpd.points_from_xy(df.Longitude, df.Latitude)
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

    gdf = gpd.sjoin(gdf, regions, how="inner", predicate="within")

    gdf.rename(columns={"name": "bus"}, inplace=True)
    gdf["Country"] = gdf.bus.str[:2]
    # filter for countries that are missing
    gdf = gdf[
        (~gdf.Country.isin(fh_ammonia.Country.unique()))
        & (gdf.Country.isin(snakemake.params.countries))
    ]
    # following approach from build_industrial_distribution_key.py
    for country in gdf.Country:
        facilities = gdf.query("Country == @country")
        production = facilities["Ammonia [kt/a]"]
        # assume 50% of the minimum production for missing values
        production = production.fillna(0.5 * facilities["Ammonia [kt/a]"].min())

    # missing data
    gdf.drop(gdf[gdf["Ammonia [kt/a]"].isna()].index, inplace=True)

    # get average plant age:
    avg_age = plant_data[plant_data.Product == "Ammonia"][
        "Year of last modernisation"
    ].mean()
    gdf["grouping_year"] = min(
        y
        for y in snakemake.params.existing_capacities["grouping_years_industry"]
        if y > avg_age
    )
    # match database
    gdf.index = (
        gdf["bus"] + " Haber-Bosch-SMR-" + gdf["grouping_year"].values.astype(str)
    )
    gdf.rename(
        columns={"Ammonia [kt/a]": "Production in tons (calibrated)"}, inplace=True
    )
    gdf["Production in tons (calibrated)"] *= 1e3

    ammonia_plants = pd.concat(
        [
            fh_ammonia,
            gdf[["bus", "Country", "grouping_year", "Production in tons (calibrated)"]],
        ]
    )

    # https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry.pdf
    # page 56: 1.83 t_CO2/t_NH3
    ch4_per_nh3 = (
        1.83 / costs.at["gas", "CO2 intensity"] / snakemake.params["MWh_NH3_per_tNH3"]
    )
    n.add(
        "Link",
        ammonia_plants.index,
        bus0=[bus + " gas" for bus in ammonia_plants.bus]
        if snakemake.params.sector["gas_network"]
        else "EU gas",
        bus1=[bus + " NH3" for bus in ammonia_plants.bus]
        if snakemake.params.sector["ammonia"]
        else "EU NH3",
        bus2=ammonia_plants.bus,
        bus3="co2 atmosphere",
        p_nom=ammonia_plants["Production in tons (calibrated)"]
        .mul(snakemake.params.MWh_NH3_per_tNH3)
        .div(ch4_per_nh3)
        .div(8760)
        .values,
        p_nom_extendable=False,
        carrier="Haber-Bosch",
        efficiency=1 / ch4_per_nh3,
        efficiency1=-costs.at["Haber-Bosch", "electricity-input"] / ch4_per_nh3,
        efficiency2=costs.at["gas", "CO2 intensity"],
        capital_cost=costs.at["Haber-Bosch", "capital_cost"]
        / costs.at["Haber-Bosch", "electricity-input"],
        marginal_cost=costs.at["Haber-Bosch", "VOM"]
        / costs.at["Haber-Bosch", "electricity-input"],
        build_year=ammonia_plants["grouping_year"],
        lifetime=costs.at["Haber-Bosch", "lifetime"],
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_existing_baseyear",
            configfiles="config/test/config.myopic.yaml",
            clusters="5",
            opts="",
            sector_opts="",
            planning_horizons=2030,
        )

    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    options = snakemake.params.sector

    baseyear = snakemake.params.baseyear

    n = pypsa.Network(snakemake.input.network)

    # define spatial resolution of carriers
    spatial = define_spatial(n.buses[n.buses.carrier == "AC"].index, options)
    add_build_year_to_new_assets(n, baseyear)

    Nyears = n.snapshot_weightings.generators.sum() / 8760.0
    costs = load_costs(
        snakemake.input.costs,
        snakemake.params.costs,
        nyears=Nyears,
    )

    grouping_years_power = snakemake.params.existing_capacities["grouping_years_power"]
    grouping_years_heat = snakemake.params.existing_capacities["grouping_years_heat"]
    add_power_capacities_installed_before_baseyear(
        n=n,
        costs=costs,
        grouping_years=grouping_years_power,
        baseyear=baseyear,
        powerplants_file=snakemake.input.powerplants,
        countries=snakemake.config["countries"],
        capacity_threshold=snakemake.params.existing_capacities["threshold_capacity"],
        lifetime_values=snakemake.params.costs["fill_values"],
    )

    if options["heating"]:
        # one could use baseyear here instead (but dangerous if no data)
        fn = snakemake.input.heating_efficiencies
        year = int(snakemake.params["energy_totals_year"])
        heating_efficiencies = pd.read_csv(fn, index_col=[1, 0]).loc[year]

        add_heating_capacities_installed_before_baseyear(
            n=n,
            costs=costs,
            baseyear=baseyear,
            grouping_years=grouping_years_heat,
            heat_pump_cop=xr.open_dataarray(snakemake.input.cop_profiles),
            use_time_dependent_cop=options["time_dep_hp_cop"],
            default_lifetime=snakemake.params.existing_capacities[
                "default_heating_lifetime"
            ],
            existing_capacities=pd.read_csv(
                snakemake.input.existing_heating_distribution,
                header=[0, 1],
                index_col=0,
            ),
            heat_pump_source_types=snakemake.params.heat_pump_sources,
            efficiency_file=snakemake.input.heating_efficiencies,
            energy_totals_year=snakemake.params["energy_totals_year"],
            capacity_threshold=snakemake.params.existing_capacities[
                "threshold_capacity"
            ],
            use_electricity_distribution_grid=options["electricity_distribution_grid"],
        )

    # Set defaults for missing missing values

    if options.get("cluster_heat_buses", False):
        cluster_heat_buses(n)

    # add existing industry plants
    if "steel" in snakemake.params.sector["endogenous_sectors"]:
        add_existing_steel_plants(n)
    if "cement" in snakemake.params.sector["endogenous_sectors"]:
        add_existing_cement_plants(n)
    if ("ammonia" in snakemake.params.sector["endogenous_sectors"]) and (snakemake.params.sector["ammonia"]):
        add_existing_ammonia_plants(n)
    if "methanol" in snakemake.params.sector["endogenous_sectors"]:
        add_existing_meoh_plants(n)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))

    sanitize_custom_columns(n)
    sanitize_carriers(n, snakemake.config)
    n.export_to_netcdf(snakemake.output[0])
