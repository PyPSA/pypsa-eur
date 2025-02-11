# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Adds existing power and heat generation capacities for initial planning
horizon.
"""

import logging
from collections import Counter
from types import SimpleNamespace

import country_converter as coco
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import xarray as xr
from _helpers import (
    configure_logging,
    sanitize_custom_columns,
    set_scenario_config,
    update_config_from_wildcards,
)
from add_electricity import sanitize_carriers, calculate_annuity
from definitions.heat_sector import HeatSector
from definitions.heat_system import HeatSystem
from prepare_sector_network import cluster_heat_buses, define_spatial, prepare_costs, calculate_steel_parameters

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

        # distribute capacities among nodes according to capacity factor
        # weighting with nodal_fraction
        elec_buses = n.buses.index[n.buses.carrier == "AC"].union(
            n.buses.index[n.buses.carrier == "DC"]
        )
        nodal_fraction = pd.Series(0.0, elec_buses)

        for country in n.buses.loc[elec_buses, "country"].unique():
            gens = n.generators.index[
                (n.generators.index.str[:2] == country)
                & (n.generators.carrier == carrier)
            ]
            cfs = n.generators_t.p_max_pu[gens].mean()
            cfs_key = cfs / cfs.sum()
            nodal_fraction.loc[n.generators.loc[gens, "bus"]] = cfs_key.groupby(
                n.generators.loc[gens, "bus"]
            ).sum()

        nodal_df = df.loc[n.buses.loc[elec_buses, "country"]]
        nodal_df.index = elec_buses
        nodal_df = nodal_df.multiply(nodal_fraction, axis=0)

        for year in nodal_df.columns:
            for node in nodal_df.index:
                name = f"{node}-{carrier}-{year}"
                capacity = nodal_df.loc[node, year]
                if capacity > 0.0:
                    cost_key = carrier.split("-")[0]
                    df_agg.at[name, "Fueltype"] = carrier
                    df_agg.at[name, "Capacity"] = capacity
                    df_agg.at[name, "DateIn"] = year
                    df_agg.at[name, "lifetime"] = costs.at[cost_key, "lifetime"]
                    df_agg.at[name, "DateOut"] = (
                        year + costs.at[cost_key, "lifetime"] - 1
                    )
                    df_agg.at[name, "bus"] = node


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
        index=["grouping_year", "Fueltype"],
        columns="bus",
        values="Capacity",
        aggfunc="sum",
    )

    lifetime = df_agg.pivot_table(
        index=["grouping_year", "Fueltype"],
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

    for grouping_year, generator in df.index:
        # capacity is the capacity in MW at each node for this
        capacity = df.loc[grouping_year, generator]
        capacity = capacity[~capacity.isna()]
        capacity = capacity[capacity > capacity_threshold]
        suffix = "-ac" if generator == "offwind" else ""
        name_suffix = f" {generator}{suffix}-{grouping_year}"
        name_suffix_by = f" {generator}{suffix}-{baseyear}"
        asset_i = capacity.index + name_suffix
        if generator in ["solar", "onwind", "offwind-ac"]:
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
            lifetime_assets = lifetime.loc[grouping_year, generator].dropna()

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
                        * costs.at[generator, "fixed"],  # NB: fixed cost is per MWel
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
                        capital_cost=costs.at[key, "fixed"]
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
        logger.warning(f"{heat_system} not defined.")

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
                    * costs.at[costs_name, "fixed"],
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
                    * costs.at[heat_system.resistive_heater_costs_name, "fixed"]
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
                    * costs.at[heat_system.gas_boiler_costs_name, "fixed"]
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
                * costs.at[heat_system.oil_boiler_costs_name, "fixed"],
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


def add_steel_industry_existing(n):

    # Steel capacities in Europe in kton of steel products per year
    capacities = pd.read_csv(snakemake.input.endoindustry_capacities, index_col=0)
    capacities = capacities[['EAF','DRI + EAF', 'Integrated steelworks']]
    start_dates = pd.read_csv(snakemake.input.endoindustry_start_dates, index_col=0)
    start_dates = start_dates[['EAF','DRI + EAF', 'Integrated steelworks']]
    keys = pd.read_csv(snakemake.input.industrial_distribution_key, index_col=0)

    capacities_bof = capacities["Integrated steelworks"]
    capacities_eaf = capacities["EAF"] + capacities["DRI + EAF"]
    capacities_bof.index = capacities.index
    capacities_eaf.index = capacities.index

    capacities_bof = capacities_bof * keys["Integrated steelworks"]

    capacities_eaf = capacities_eaf * keys["EAF"]

    start_dates_eaf = round((start_dates["EAF"] * capacities["EAF"] + start_dates["DRI + EAF"] * capacities["DRI + EAF"])/ capacities_eaf)
    start_dates_bof = round(start_dates["Integrated steelworks"])

    # Average age of assets in Iron and steel in Europe: 21-28 years, so I assume they are starting in 2000 in case https://www.energimyndigheten.se/4a9556/globalassets/energieffektivisering_/jag-ar-saljare-eller-tillverkare/dokument/produkter-med-krav/ugnar-industriella-och-laboratorie/annex-b_lifetime_energy.pdf
    start_dates_eaf = start_dates_eaf.where((start_dates_eaf >= 1000) & np.isfinite(start_dates_eaf), 2000)
    start_dates_bof = start_dates_bof.where((start_dates_bof >= 1000) & np.isfinite(start_dates_bof), 2000)

    nodes = pop_layout.index
    p_nom_bof = pd.DataFrame(index=nodes, columns=(["value"]))
    p_nom_eaf = pd.DataFrame(index=nodes, columns=(["value"]))

    p_nom_bof = capacities_bof / nhours  # get the hourly production capacity
    p_nom_eaf = capacities_eaf / nhours  # get the hourly production capacity

    # Should steel be produced at a constant rate during the year or not? 1 or 0
    prod_constantly = 0
    ramp_limit = 0

    # PARAMETERS
    nyears = n.snapshot_weightings.generators.sum() / 8760.0
    bof, eaf_ng, eaf_h2, tgr = calculate_steel_parameters(nyears)

    n.add(
        "Link",
        nodes,
        suffix=" BF-BOF-2020",
        bus0=spatial.iron.nodes,
        bus1=spatial.steel.nodes,
        bus2=spatial.coal.nodes,
        bus3=nodes,
        bus4=spatial.co2.bof,
        carrier="BF-BOF",
        p_nom=p_nom_bof * bof['iron input'],
        p_nom_extendable=False,
        marginal_cost=0,#opex_bof,
        efficiency=1 / bof['iron input'],
        efficiency2= -  bof['coal input'] /  bof['iron input'],  # MWhth coal per kt iron
        efficiency3= -  bof['elec input'] /  bof['iron input'],  # MWh electricity per kt iron
        efficiency4=  bof['emission factor'] /  bof['iron input'], # t CO2 per kt iron
        lifetime= bof['lifetime'],  # https://www.energimyndigheten.se/4a9556/globalassets/energieffektivisering_/jag-ar-saljare-eller-tillverkare/dokument/produkter-med-krav/ugnar-industriella-och-laboratorie/annex-b_lifetime_energy.pdf
        build_year=start_dates_bof,
    )

    n.add(
        "Link",
        nodes,
        suffix=" DRI-EAF-2020",
        bus0=spatial.iron.nodes,
        bus1=spatial.steel.nodes,
        bus2=spatial.syngas_dri.nodes,  # in this process is the reducing agent, it is not burnt
        bus3=nodes,
        carrier="DRI-EAF",
        p_nom=p_nom_eaf *  eaf_ng['iron input'],
        p_nom_extendable=False,
        marginal_cost=0,#opex_eaf,
        efficiency=1 / eaf_ng['iron input'],
        efficiency2= -1 / eaf_ng['iron input'], # one unit of dri gas per kt iron
        efficiency3= - eaf_ng['elec input'] / eaf_ng['iron input'], #MWh electricity per kt iron
        lifetime= eaf_ng['lifetime'],  # https://www.energimyndigheten.se/4a9556/globalassets/energieffektivisering_/jag-ar-saljare-eller-tillverkare/dokument/produkter-med-krav/ugnar-industriella-och-laboratorie/annex-b_lifetime_energy.pdf
        build_year=start_dates_eaf,
    )


def add_cement_industry_existing(n):

    # Cement capacities in Europe in kton of cement products per year
    capacities = pd.read_csv(snakemake.input.endoindustry_capacities, index_col=0)
    capacities = capacities['Cement']
    start_dates = pd.read_csv(snakemake.input.endoindustry_start_dates, index_col=0)
    start_dates = round(start_dates['Cement'])
    keys = pd.read_csv(snakemake.input.industrial_distribution_key, index_col=0)

    capacities = capacities * keys["Cement_SFI"]

    start_dates = start_dates.where(
        (start_dates >= 1000) & np.isfinite(start_dates), 2000
    )

    nodes = pop_layout.index
    p_nom = pd.DataFrame(index=nodes, columns=(["value"]))

    p_nom = capacities / nhours  # get the hourly production capacity

    # Should steel be produced at a constant rate during the year or not? 1 or 0
    prod_constantly = 0
    ramp_limit = 0

    ########### Add existing cement production capacities ############

    # Lifetimes
    lifetime_cement = 100

    # Capital costs
    discount_rate = 0.04
    capex_cement = 263000/nhours * calculate_annuity(lifetime_cement, discount_rate) # https://iea-etsap.org/E-TechDS/HIGHLIGHTS%20PDF/I03_cement_June%202010_GS-gct%201.pdf with CCS 558000 

    n.add(
        "Link",
        nodes,
        suffix=" Cement Plant-2020",
        bus0=spatial.limestone.nodes,
        bus1=spatial.cement.nodes,
        bus2=spatial.gas.nodes,
        bus3=spatial.co2.cement,
        carrier="cement plant",
        p_nom=p_nom,
        p_min_pu=prod_constantly,  # hot elements cannot be turned off easily
        p_nom_extendable=False,
        capital_cost=capex_cement,
        efficiency=1/1.28, # kt limestone/ kt clinker https://www.sciencedirect.com/science/article/pii/S2214157X22005974
        efficiency2= - 3420.1 / 3.6 * (1/1.28) / 0.5, # MWh/kt clinker https://www.sciencedirect.com/science/article/pii/S2214157X22005974
        efficiency3=500 * (1/1.28), #tCO2/kt cement
        lifetime=lifetime_cement, 
        build_year=start_dates,
    )


def add_chemicals_industry_existing(n, options):

    # Chemicals capacities in Europe in kton of cement products per year
    capacities = pd.read_csv(snakemake.input.endoindustry_capacities, index_col=0)
    start_dates = pd.read_csv(snakemake.input.endoindustry_start_dates, index_col=0)
    keys = pd.read_csv(snakemake.input.industrial_distribution_key, index_col=0)

    # Ammonia
    capacities_nh3 = capacities['Ammonia']
    start_dates_nh3 = start_dates['Ammonia']
    capacities_nh3 = capacities_nh3 * keys["Ammonia"]

    start_dates_nh3 = round(start_dates_nh3)
    start_dates_nh3 = start_dates_nh3.where((start_dates_nh3 >= 1000) & np.isfinite(start_dates_nh3), 2000)

    nodes = pop_layout.index
    p_nom_nh3 = pd.DataFrame(index=nodes, columns=(["value"]))

    p_nom_nh3 = capacities_nh3 / nhours  # get the hourly production capacity

    ########### Add existing ammonia production capacities ############

    n.add(
        "Link",
        nodes,
        suffix=" Haber-Bosch-2020",
        bus0=nodes,
        bus1=spatial.ammonia.nodes,
        bus2=nodes + " H2",
        p_nom_extendable=False,
        p_nom=p_nom_nh3,
        #p_min_pu=prod_constantly,  # hot elements cannot be turned off easily
        carrier="Haber-Bosch",
        efficiency=1 / costs.at["Haber-Bosch", "electricity-input"],
        efficiency2=-costs.at["Haber-Bosch", "hydrogen-input"]
        / costs.at["Haber-Bosch", "electricity-input"],
        capital_cost=costs.at["Haber-Bosch", "fixed"]
        / costs.at["Haber-Bosch", "electricity-input"],
        marginal_cost=costs.at["Haber-Bosch", "VOM"]
        / costs.at["Haber-Bosch", "electricity-input"],
        lifetime=costs.at["Haber-Bosch", "lifetime"],
        build_year=start_dates_nh3,
    )

    # Methanol

    capacities_meth = capacities['Methanol']
    start_dates_meth = start_dates['Methanol']
    capacities_meth = capacities_meth * keys["Chemical industry"] #ADB fix this with real methanol

    start_dates_meth = round(start_dates_meth)
    start_dates_meth = start_dates_meth.where((start_dates_meth >= 1000) & np.isfinite(start_dates_meth), 2000)

    p_nom_meth = pd.DataFrame(index=nodes, columns=(["value"]))

    p_nom_meth = capacities_meth / nhours  # get the hourly production capacity

    ########### Add existing methanol production capacities ############

    n.add(
        "Link",
        nodes,
        suffix = " methanolisation-2020",
        #spatial.h2.locations + " methanolisation-2020",
        bus0=spatial.h2.nodes,
        bus1=spatial.methanol.nodes,
        bus2=nodes,
        bus3=spatial.co2.nodes,
        carrier="methanolisation",
        p_nom_extendable=False,
        p_nom=p_nom_meth,
        p_min_pu=options["min_part_load_methanolisation"],
        capital_cost=costs.at["methanolisation", "fixed"]
        * options["MWh_MeOH_per_MWh_H2"],  # EUR/MW_H2/a
        marginal_cost=options["MWh_MeOH_per_MWh_H2"]
        * costs.at["methanolisation", "VOM"],
        lifetime=costs.at["methanolisation", "lifetime"],
        efficiency=options["MWh_MeOH_per_MWh_H2"],
        efficiency2=-options["MWh_MeOH_per_MWh_H2"] / options["MWh_MeOH_per_MWh_e"],
        efficiency3=-options["MWh_MeOH_per_MWh_H2"] / options["MWh_MeOH_per_tCO2"],
        build_year=start_dates_meth,
    )

    # HVC (Ethylene)

    if options["endo_industry"]["endo_hvc"]:
        capacities_hvc = capacities['Ethylene']
        start_dates_hvc = start_dates['Methanol']
        capacities_hvc = capacities_hvc * keys["Chemical industry"] #ADB fix this with real hvc

        start_dates_hvc = round(start_dates_hvc)
        start_dates_hvc = start_dates_meth.where((start_dates_meth >= 1000) & np.isfinite(start_dates_meth), 2000)

        p_nom_hvc = pd.DataFrame(index=nodes, columns=(["value"]))

        p_nom_hvc = capacities_hvc / nhours  # get the hourly production capacity

        ########### Add existing HVC production capacities ############

        naphtha_to_hvc = 2.31 * 12.47 * 1000 # kt oil / kt HVC * MWh/t oil * 1000 t / kt =   MWh oil / kt HVC

        n.add(
            "Link",
            spatial.hvc.locations,
            suffix = " naphtha steam cracker-2020",
            bus0=spatial.oil.nodes,
            bus1=spatial.hvc.nodes,
            bus2=nodes + " H2",
            bus3="co2 atmosphere",
            bus4=nodes,
            carrier="naphtha steam cracker",
            p_nom_extendable=False,
            capital_cost=725 * 1e3, #€/kt HVC
            efficiency=1/ naphtha_to_hvc, # MWh oil / kt HVC
            efficiency2= 21 * 33.3 / naphtha_to_hvc, # MWh H2 / kt HVC
            efficiency3= 819 / naphtha_to_hvc, # tCO2 / kt HVC
            efficiency4= - 135 / naphtha_to_hvc, # MWh electricity / kt HVC
            lifetime=30, 
            build_year=start_dates_hvc,
        )


def set_defaults(n):
    """
    Set default values for missing values in the network.

    Parameters
    ----------
        n (pypsa.Network): The network object.

    Returns
    -------
        None
    """
    if "Link" in n.components:
        if "reversed" in n.links.columns:
            # Replace NA values with default value False
            n.links.loc[n.links.reversed.isna(), "reversed"] = False
            n.links.reversed = n.links.reversed.astype(bool)


# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_existing_baseyear",
            configfiles="config/test/config.myopic.yaml",
            clusters="5",
            ll="v1.5",
            opts="",
            sector_opts="",
            planning_horizons=2030,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    options = snakemake.params.sector

    baseyear = snakemake.params.baseyear

    n = pypsa.Network(snakemake.input.network)

    # define spatial resolution of carriers
    spatial = define_spatial(n.buses[n.buses.carrier == "AC"].index, options)
    add_build_year_to_new_assets(n, baseyear)

    Nyears = n.snapshot_weightings.generators.sum() / 8760.0
    costs = prepare_costs(
        snakemake.input.costs,
        snakemake.params.costs,
        Nyears,
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

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    nhours = n.snapshot_weightings.generators.sum()
    if options["endo_industry"].get("enable"):
        add_steel_industry_existing(n)
        add_cement_industry_existing(n)
        if options["endo_industry"].get("endo_chemicals"):
            add_chemicals_industry_existing(n, options)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))

    sanitize_custom_columns(n)
    sanitize_carriers(n, snakemake.config)
    n.export_to_netcdf(snakemake.output[0])
