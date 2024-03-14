# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Adds existing power and heat generation capacities for initial planning
horizon.
"""

import logging
import os
import sys
from types import SimpleNamespace

import country_converter as coco
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from add_electricity import sanitize_carriers
from prepare_sector_network import cluster_heat_buses, define_spatial, prepare_costs

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()
idx = pd.IndexSlice
spatial = SimpleNamespace()

from build_powerplants import add_custom_powerplants


def add_build_year_to_new_assets(n, baseyear):
    """
    Parameters
    ----------
    n : pypsa.Network
    baseyear : int
        year in which optimized assets are built
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


def add_existing_renewables(df_agg):
    """
    Append existing renewables to the df_agg pd.DataFrame with the conventional
    power plants.
    """
    carriers = {"solar": "solar", "onwind": "onwind", "offwind": "offwind-ac"}

    for tech in ["solar", "onwind", "offwind"]:
        carrier = carriers[tech]

        df = pd.read_csv(snakemake.input[f"existing_{tech}"], index_col=0).fillna(0.0)
        df.columns = df.columns.astype(int)
        df.index = cc.convert(df.index, to="iso2")

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
                name = f"{node}-{tech}-{year}"
                capacity = nodal_df.loc[node, year]
                if capacity > 0.0:
                    df_agg.at[name, "Fueltype"] = tech
                    df_agg.at[name, "Capacity"] = capacity
                    df_agg.at[name, "DateIn"] = year
                    df_agg.at[name, "cluster_bus"] = node


def add_power_capacities_installed_before_baseyear(n, grouping_years, costs, baseyear):
    """
    Parameters
    ----------
    n : pypsa.Network
    grouping_years :
        intervals to group existing capacities
    costs :
        to read lifetime to estimate YearDecomissioning
    baseyear : int
    """
    logger.debug(
        f"Adding power capacities installed before {baseyear} from powerplants.csv"
    )

    df_agg = pd.read_csv(snakemake.input.powerplants, index_col=0)

    if snakemake.input.get("custom_powerplants"):
        df_agg = add_custom_powerplants(
            df_agg, snakemake.input.custom_powerplants, True
        )

    rename_fuel = {
        "Hard Coal": "coal",
        "Lignite": "lignite",
        "Nuclear": "nuclear",
        "Oil": "oil",
        "OCGT": "OCGT",
        "CCGT": "CCGT",
        "Bioenergy": "urban central solid biomass OP",
    }

    # If heat is considered, add CHPs in the add_heating_capacities function.
    # Assume that all oil power plants are not CHPs.
    if "H" in snakemake.wildcards.sector_opts.split("-"):
        df_agg = df_agg.query("Set != 'CHP'")
    elif "I" not in snakemake.wildcards.sector_opts.split("-") and "Industry" in df_agg.columns:
         df_agg.query("not Industry", inplace=True)

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

    # drop unused fueltyps and technologies
    df_agg.drop(df_agg.index[df_agg.Fueltype.isin(fueltype_to_drop)], inplace=True)
    df_agg.drop(df_agg.index[df_agg.Technology.isin(technology_to_drop)], inplace=True)
    df_agg.Fueltype = df_agg.Fueltype.map(rename_fuel)

    # Intermediate fix for DateIn & DateOut
    # Fill missing DateIn
    biomass_i = df_agg.loc[df_agg.Fueltype == "urban central solid biomass OP"].index
    mean = df_agg.loc[biomass_i, "DateIn"].mean()
    df_agg.loc[biomass_i, "DateIn"] = df_agg.loc[biomass_i, "DateIn"].fillna(int(mean))
    # Fill missing DateOut
    dateout = (
        df_agg.loc[biomass_i, "DateIn"]
        + snakemake.params.costs["fill_values"]["lifetime"]
    )
    df_agg.loc[biomass_i, "DateOut"] = df_agg.loc[biomass_i, "DateOut"].fillna(dateout)

    # drop assets which are already phased out / decommissioned
    phased_out = df_agg[df_agg["DateOut"] < baseyear].index
    df_agg.drop(phased_out, inplace=True)

    # assign clustered bus
    busmap_s = pd.read_csv(snakemake.input.busmap_s, index_col=0).squeeze()
    busmap = pd.read_csv(snakemake.input.busmap, index_col=0).squeeze()

    inv_busmap = {}
    for k, v in busmap.items():
        inv_busmap[v] = inv_busmap.get(v, []) + [k]

    clustermaps = busmap_s.map(busmap)
    clustermaps.index = clustermaps.index.astype(int)

    df_agg["cluster_bus"] = df_agg.bus.map(clustermaps)

    # include renewables in df_agg
    add_existing_renewables(df_agg)

    df_agg["grouping_year"] = np.take(
        grouping_years, np.digitize(df_agg.DateIn, grouping_years, right=True)
    )

    # calculate (adjusted) remaining lifetime before phase-out (+1 because assuming
    # phase out date at the end of the year)
    df_agg["lifetime"] = df_agg.DateOut - df_agg["grouping_year"] + 1

    df = df_agg.pivot_table(
        index=["grouping_year", "Fueltype"],
        columns="cluster_bus",
        values="Capacity",
        aggfunc="sum",
    )

    lifetime = df_agg.pivot_table(
        index=["grouping_year", "Fueltype"],
        columns="cluster_bus",
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
        "urban central solid biomass OP": "biomass",
    }

    for grouping_year, generator in df.index:
        # capacity is the capacity in MW at each node for this
        capacity = df.loc[grouping_year, generator]
        capacity = capacity[~capacity.isna()]
        capacity = capacity[
            capacity > snakemake.params.existing_capacities["threshold_capacity"]
        ]
        suffix = "-ac" if generator == "offwind" else ""
        name_suffix = f" {generator}{suffix}-{grouping_year}"
        asset_i = capacity.index + name_suffix
        if generator in ["solar", "onwind", "offwind"]:
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
                n.generators.loc[already_build, "p_nom_min"] = capacity.loc[
                    already_build.str.replace(name_suffix, "")
                ].values
            new_capacity = capacity.loc[new_build.str.replace(name_suffix, "")]

            if "m" in snakemake.wildcards.clusters:
                for ind in new_capacity.index:
                    # existing capacities are split evenly among regions in every country
                    inv_ind = list(inv_busmap[ind])

                    # for offshore the splitting only includes coastal regions
                    inv_ind = [
                        i for i in inv_ind if (i + name_suffix) in n.generators.index
                    ]

                    p_max_pu = n.generators_t.p_max_pu[
                        [i + name_suffix for i in inv_ind]
                    ]
                    p_max_pu.columns = [i + name_suffix for i in inv_ind]

                    n.madd(
                        "Generator",
                        [i + name_suffix for i in inv_ind],
                        bus=ind,
                        carrier=generator,
                        p_nom=new_capacity[ind]
                        / len(inv_ind),  # split among regions in a country
                        marginal_cost=marginal_cost,
                        capital_cost=capital_cost,
                        efficiency=costs.at[generator, "efficiency"],
                        p_max_pu=p_max_pu,
                        build_year=grouping_year,
                        lifetime=costs.at[generator, "lifetime"],
                    )

            else:
                p_max_pu = n.generators_t.p_max_pu[
                    capacity.index + f" {generator}{suffix}-{baseyear}"
                ]

                if not new_build.empty:
                    n.madd(
                        "Generator",
                        new_capacity.index,
                        suffix=" " + name_suffix,
                        bus=new_capacity.index,
                        carrier=generator,
                        p_nom=new_capacity,
                        marginal_cost=marginal_cost,
                        capital_cost=capital_cost,
                        efficiency=costs.at[generator, "efficiency"],
                        p_max_pu=p_max_pu.rename(columns=n.generators.bus),
                        build_year=grouping_year,
                        lifetime=costs.at[generator, "lifetime"],
                    )

        else:
            bus0 = vars(spatial)[carrier[generator]].nodes
            if "EU" not in vars(spatial)[carrier[generator]].locations:
                bus0 = bus0.intersection(capacity.index + " " + carrier[generator])

            # check for missing bus
            missing_bus = pd.Index(bus0).difference(n.buses.index)
            if not missing_bus.empty:
                logger.info(f"add buses {bus0}")
                n.madd(
                    "Bus",
                    bus0,
                    carrier=generator,
                    location=vars(spatial)[carrier[generator]].locations,
                    unit="MWh_el",
                )

            already_build = n.links.index.intersection(asset_i)
            new_build = asset_i.difference(n.links.index)
            if (grouping_year, generator) in lifetime.index:
                lifetime_assets = lifetime.loc[grouping_year, generator].dropna()
            else:
                lifetime_assets = costs.at[generator, "lifetime"]

            # this is for the year 2020
            if not already_build.empty:
                n.links.loc[already_build, "p_nom_min"] = capacity.loc[
                    already_build.str.replace(name_suffix, "")
                ].values

            if not new_build.empty:
                new_capacity = capacity.loc[new_build.str.replace(name_suffix, "")]

                if generator != "urban central solid biomass OP":
                    # missing lifetimes are filled with mean lifetime
                    # if mean cannot be built, lifetime is taken from costs.csv
                    if isinstance(lifetime_assets, pd.Series):
                        lifetime_assets = (
                            lifetime_assets.reindex(capacity.index)
                            .fillna(lifetime_assets.mean())
                            .fillna(costs.at[generator, "lifetime"])
                        )

                    n.madd(
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
                        lifetime=(
                            lifetime_assets.loc[new_capacity.index]
                            if isinstance(lifetime_assets, pd.Series)
                            else lifetime_assets
                        ),
                    )
                else:
                    # for the power only biomass plants, technology parameters of CHP plants are used
                    key = "central solid biomass CHP"
                    n.madd(
                        "Link",
                        new_capacity.index,
                        suffix=name_suffix,
                        bus0=spatial.biomass.df.loc[new_capacity.index]["nodes"].values,
                        bus1=new_capacity.index,
                        carrier=generator,
                        p_nom=new_capacity / costs.at[key, "efficiency"],
                        capital_cost=costs.at[key, "fixed"]
                        * costs.at[key, "efficiency"],
                        marginal_cost=costs.at[key, "VOM"],
                        efficiency=costs.at[key, "efficiency"],
                        build_year=grouping_year,
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


def add_heating_capacities_installed_before_baseyear(
    n,
    baseyear,
    grouping_years,
    ashp_cop,
    gshp_cop,
    time_dep_hp_cop,
    costs,
    default_lifetime,
):
    """
    Parameters
    ----------
    n : pypsa.Network
    baseyear : last year covered in the existing capacities database
    grouping_years : intervals to group existing capacities
        linear decommissioning of heating capacities from 2020 to 2045 is
        currently assumed heating capacities split between residential and
        services proportional to heating load in both 50% capacities
        in rural busess 50% in urban buses
    """
    logger.debug(f"Adding heating capacities installed before {baseyear}")

    existing_heating = pd.read_csv(
        snakemake.input.existing_heating_distribution, header=[0, 1], index_col=0
    )

    techs = existing_heating.columns.get_level_values(1).unique()

    for name in existing_heating.columns.get_level_values(0).unique():
        name_type = "central" if name == "urban central" else "decentral"

        nodes = pd.Index(n.buses.location[n.buses.index.str.contains(f"{name} heat")])

        if (name_type != "central") and options["electricity_distribution_grid"]:
            nodes_elec = nodes + " low voltage"
        else:
            nodes_elec = nodes

        heat_pump_type = "air" if "urban" in name else "ground"

        # Add heat pumps
        costs_name = f"decentral {heat_pump_type}-sourced heat pump"

        cop = {"air": ashp_cop, "ground": gshp_cop}

        if time_dep_hp_cop:
            efficiency = cop[heat_pump_type][nodes]
        else:
            efficiency = costs.at[costs_name, "efficiency"]

        for i, grouping_year in enumerate(grouping_years):
            if int(grouping_year) + default_lifetime <= int(baseyear):
                continue

            # installation is assumed to be linear for the past default_lifetime years
            ratio = (int(grouping_year) - int(grouping_years[i - 1])) / default_lifetime

            n.madd(
                "Link",
                nodes,
                suffix=f" {name} {heat_pump_type} heat pump-{grouping_year}",
                bus0=nodes_elec,
                bus1=nodes + " " + name + " heat",
                carrier=f"{name} {heat_pump_type} heat pump",
                efficiency=efficiency,
                capital_cost=costs.at[costs_name, "efficiency"]
                * costs.at[costs_name, "fixed"],
                p_nom=existing_heating.loc[nodes, (name, f"{heat_pump_type} heat pump")]
                * ratio
                / costs.at[costs_name, "efficiency"],
                build_year=int(grouping_year),
                lifetime=costs.at[costs_name, "lifetime"],
            )

            # add resistive heater, gas boilers and oil boilers
            n.madd(
                "Link",
                nodes,
                suffix=f" {name} resistive heater-{grouping_year}",
                bus0=nodes_elec,
                bus1=nodes + " " + name + " heat",
                carrier=name + " resistive heater",
                efficiency=costs.at[f"{name_type} resistive heater", "efficiency"],
                capital_cost=(
                    costs.at[f"{name_type} resistive heater", "efficiency"]
                    * costs.at[f"{name_type} resistive heater", "fixed"]
                ),
                p_nom=(
                    existing_heating.loc[nodes, (name, "resistive heater")]
                    * ratio
                    / costs.at[f"{name_type} resistive heater", "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at[f"{name_type} resistive heater", "lifetime"],
            )

            n.madd(
                "Link",
                nodes,
                suffix=f" {name} gas boiler-{grouping_year}",
                bus0="EU gas" if "EU gas" in spatial.gas.nodes else nodes + " gas",
                bus1=nodes + " " + name + " heat",
                bus2="co2 atmosphere",
                carrier=name + " gas boiler",
                efficiency=costs.at[f"{name_type} gas boiler", "efficiency"],
                efficiency2=costs.at["gas", "CO2 intensity"],
                capital_cost=(
                    costs.at[f"{name_type} gas boiler", "efficiency"]
                    * costs.at[f"{name_type} gas boiler", "fixed"]
                ),
                p_nom=(
                    existing_heating.loc[nodes, (name, "gas boiler")]
                    * ratio
                    / costs.at[f"{name_type} gas boiler", "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at[f"{name_type} gas boiler", "lifetime"],
            )

            n.madd(
                "Link",
                nodes,
                suffix=f" {name} oil boiler-{grouping_year}",
                bus0=spatial.oil.nodes,
                bus1=nodes + " " + name + " heat",
                bus2="co2 atmosphere",
                carrier=name + " oil boiler",
                efficiency=costs.at["decentral oil boiler", "efficiency"],
                efficiency2=costs.at["oil", "CO2 intensity"],
                capital_cost=costs.at["decentral oil boiler", "efficiency"]
                * costs.at["decentral oil boiler", "fixed"],
                p_nom=(
                    existing_heating.loc[nodes, (name, "oil boiler")]
                    * ratio
                    / costs.at["decentral oil boiler", "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at[f"{name_type} gas boiler", "lifetime"],
            )
            # add biomass boilers
            n.madd(
                "Link",
                nodes,
                suffix=f" {name} biomass boiler-{grouping_year}",
                bus0=spatial.biomass.nodes,
                bus1=nodes + " " + name + " heat",
                carrier=name + " biomass boiler",
                efficiency=costs.at["biomass boiler", "efficiency"],
                capital_cost=costs.at["biomass boiler", "efficiency"]
                * costs.at["biomass boiler", "fixed"],
                p_nom=(
                    existing_heating.loc[nodes, (name, "biomass boiler")]
                    * ratio
                    / costs.at["biomass boiler", "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at["biomass boiler", "lifetime"],
            )

            # delete links with p_nom=nan corresponding to extra nodes in country
            n.mremove(
                "Link",
                [
                    index
                    for index in n.links.index.to_list()
                    if str(grouping_year) in index and np.isnan(n.links.p_nom[index])
                ],
            )

            # delete links with capacities below threshold
            threshold = snakemake.params.existing_capacities["threshold_capacity"]
            n.mremove(
                "Link",
                [
                    index
                    for index in n.links.index.to_list()
                    if str(grouping_year) in index and n.links.p_nom[index] < threshold
                ],
            )

    # deal with CHP plants
    # rename fuel of CHPs - lignite not in DEA database
    rename_fuel = {
        "Hard Coal": "coal",
        "Coal": "coal",
        "Lignite": "lignite",
        "Natural Gas": "gas",
        "Bioenergy": "urban central solid biomass CHP",
        "Oil": "oil",
    }

    ppl = pd.read_csv(snakemake.input.powerplants, index_col=0)

    if snakemake.input.get("custom_powerplants"):
        ppl = add_custom_powerplants(ppl, snakemake.input.custom_powerplants, True)

    # drop assets which are already phased out / decommissioned
    # drop hydro, waste and oil fueltypes for CHP
    limit = np.max(grouping_years)
    drop_fueltypes = ["Hydro", "Other", "Waste", "nicht biogener Abfall"]
    chp = ppl.query(
        "Set == 'CHP' and (DateOut >= @baseyear or DateOut != DateOut) and (DateIn <= @limit or DateIn != DateIn) and Fueltype not in @drop_fueltypes"
    ).copy()
    # chp = ppl.query(
    #     "Set == 'CHP' and DateIn <= @limit and Fueltype not in @drop_fueltypes"
    # ).copy()
    # is_null_condition = ppl['DateOut'].isna()
    # chp = chp[is_null_condition | (chp.DateOut > baseyear)]

    # calculate remaining lifetime before phase-out (+1 because assuming
    # phase out date at the end of the year)
    chp["lifetime"] = chp.DateOut - chp.DateIn + 1
    chp.Fueltype = chp.Fueltype.map(rename_fuel)

    # assign clustered bus
    busmap_s = pd.read_csv(snakemake.input.busmap_s, index_col=0).squeeze()
    busmap = pd.read_csv(snakemake.input.busmap, index_col=0).squeeze()

    clustermaps = busmap_s.map(busmap)
    clustermaps.index = clustermaps.index.astype(int)
    chp["bus"] = chp["bus"].astype(int)

    chp["cluster_bus"] = chp.bus.map(clustermaps)

    chp["grouping_year"] = np.take(
        grouping_years, np.digitize(chp.DateIn, grouping_years, right=True)
    )

    # check if the CHPs were read in from MaStR for Germany
    if "Capacity_thermal" in chp.columns:
        if "I" not in snakemake.wildcards.sector_opts.split("-"):
            chp.query("Industry == False", inplace=True)

        thermal_capacity_b = ~chp.Capacity_thermal.isna()
        mastr_chp = chp[thermal_capacity_b]

        # CHPs without thermal capacity are handled later
        chp = chp[~thermal_capacity_b]

        # exclude small CHPs below 500 kW
        mastr_chp = mastr_chp.query("Capacity > 0.5 or Capacity_thermal > 0.5")

        # separate CHPs with substantial power output from those with little power output
        # ratio chosen for reasonable backpressure coefficients c_b
        mastr_chp_power = mastr_chp.query("Capacity > 0.5 * Capacity_thermal").copy()
        mastr_chp_heat = mastr_chp.query("Capacity <= 0.5 * Capacity_thermal").copy()

        mastr_chp_power["p_nom"] = mastr_chp_power.eval("Capacity / Efficiency")
        mastr_chp_power["c_b"] = mastr_chp_power.eval("Capacity / Capacity_thermal")
        mastr_chp_power["c_b"] = mastr_chp_power["c_b"].clip(
            upper=costs.at["CCGT", "c_b"]
        )  # exclude outliers
        mastr_chp_power["efficiency-heat"] = mastr_chp_power.eval("Efficiency / c_b")

        # these CHPs are mainly biomass CHPs
        mastr_chp_heat["efficiency-heat"] = costs.at[
            "central solid biomass CHP", "efficiency-heat"
        ]
        mastr_chp_heat["p_nom"] = (
            mastr_chp_heat.Capacity_thermal / mastr_chp_heat["efficiency-heat"]
        )
        mastr_chp_heat["Efficiency"] = mastr_chp_heat.eval("Capacity / p_nom")
        eff_total_max = costs.loc[
            "central solid biomass CHP", ["efficiency-heat", "efficiency"]
        ].sum()
        eff_heat = mastr_chp_heat["efficiency-heat"]
        mastr_chp_heat["Efficiency"] = mastr_chp_heat["Efficiency"].clip(
            upper=eff_total_max - eff_heat
        )

        mastr_chp = pd.concat([mastr_chp_power, mastr_chp_heat])

        mastr_chp_efficiency_power = mastr_chp.pivot_table(
            index=["grouping_year", "Fueltype"],
            columns="cluster_bus",
            values="Efficiency",
            aggfunc=lambda x: np.average(x, weights=mastr_chp.loc[x.index, "p_nom"]),
        )

        mastr_chp_efficiency_heat = mastr_chp.pivot_table(
            index=["grouping_year", "Fueltype"],
            columns="cluster_bus",
            values="efficiency-heat",
            aggfunc=lambda x: np.average(x, weights=mastr_chp.loc[x.index, "p_nom"]),
        )

        mastr_chp_p_nom = mastr_chp.pivot_table(
            index=["grouping_year", "Fueltype"],
            columns="cluster_bus",
            values="p_nom",
            aggfunc="sum",
        )

        keys = {
            "coal": "central coal CHP",
            "gas": "central gas CHP",
            "waste": "waste CHP",
            "oil": "central gas CHP",
            "lignite": "central coal CHP",
        }
        # add everything as Link
        for grouping_year, generator in mastr_chp_p_nom.index:
            # capacity is the capacity in MW at each node for this
            p_nom = mastr_chp_p_nom.loc[grouping_year, generator].dropna()
            threshold = snakemake.params.existing_capacities["threshold_capacity"]
            p_nom = p_nom[p_nom > threshold]

            # e.g. DE1 0 coal CHP-1980
            asset_i = p_nom.index + f" {generator} CHP-{grouping_year}"

            efficiency_power = mastr_chp_efficiency_power.loc[grouping_year, generator]
            efficiency_heat = mastr_chp_efficiency_heat.loc[grouping_year, generator]

            if generator != "urban central solid biomass CHP":
                # lignite CHPs are not in DEA database - use coal CHP parameters
                key = keys[generator]
                n.madd(
                    "Link",
                    asset_i,
                    bus0=vars(spatial)[generator].nodes,  # EU gas/coal/lignite EU
                    bus1=p_nom.index,
                    bus2=p_nom.index + " urban central heat",
                    bus3="co2 atmosphere",
                    carrier=f"urban central {generator} CHP",
                    p_nom=p_nom.values,
                    capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"],
                    marginal_cost=costs.at[key, "VOM"],
                    efficiency=efficiency_power.dropna(),
                    efficiency2=efficiency_heat.dropna(),
                    efficiency3=costs.at[generator, "CO2 intensity"],
                    build_year=grouping_year,
                    lifetime=costs.at[key, "lifetime"],
                )
            else:
                key = "central solid biomass CHP"
                n.madd(
                    "Link",
                    p_nom.index,
                    suffix=key,
                    bus0=spatial.biomass.df.loc[p_nom.index]["nodes"].values,
                    bus1=p_nom.index,
                    bus2=p_nom.index + " urban central heat",
                    carrier=generator,
                    p_nom=p_nom.values,
                    capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"],
                    marginal_cost=costs.at[key, "VOM"],
                    efficiency=efficiency_power,
                    efficiency2=efficiency_heat,
                    build_year=grouping_year,
                    lifetime=costs.at[key, "lifetime"],
                )

    # CHPs that are not from MaStR
    chp_nodal_p_nom = chp.pivot_table(
        index=["grouping_year", "Fueltype"],
        columns="cluster_bus",
        values="Capacity",
        aggfunc="sum",
    )
    for grouping_year, generator in chp_nodal_p_nom.index:
        p_nom = chp_nodal_p_nom.loc[grouping_year, generator].dropna()
        threshold = snakemake.params.existing_capacities["threshold_capacity"]
        p_nom = p_nom[p_nom > threshold]

        # e.g. DE1 0 coal CHP-1980
        asset_i = p_nom.index + f" {generator} CHP-{grouping_year}"

        if generator != "urban central solid biomass CHP":
            # lignite CHPs are not in DEA database - use coal CHP parameters
            key = keys[generator]

            n.madd(
                "Link",
                asset_i,
                bus0=vars(spatial)[generator].nodes,  # EU gas/coal/lignite EU
                bus1=p_nom.index,
                bus2=p_nom.index + " urban central heat",
                bus3="co2 atmosphere",
                carrier=f"urban central {generator} CHP",
                p_nom=p_nom.values / costs.at[key, "efficiency"],
                capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"],
                marginal_cost=costs.at[key, "VOM"],
                efficiency=costs.at[key, "efficiency"],
                efficiency2=costs.at[key, "efficiency"] / costs.at[key, "c_b"],
                efficiency3=costs.at[generator, "CO2 intensity"],
                build_year=grouping_year,
                lifetime=costs.at[key, "lifetime"],
            )
        else:
            key = "central solid biomass CHP"
            n.madd(
                "Link",
                p_nom.index,
                suffix=" " + generator + "-" + str(grouping_year),
                bus0=spatial.biomass.df.loc[p_nom.index]["nodes"].values,
                bus1=p_nom.index,
                bus2=p_nom.index + " urban central heat",
                carrier=generator,
                p_nom=p_nom.values / costs.at[key, "efficiency"],
                capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"],
                marginal_cost=costs.at[key, "VOM"],
                efficiency=costs.at[key, "efficiency"],
                efficiency2=costs.at[key, "efficiency-heat"],
                build_year=grouping_year,
                lifetime=costs.at[key, "lifetime"],
            )


# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_existing_baseyear",
            # configfiles="config/test/config.myopic.yaml",
            simpl="",
            clusters="22",
            ll="v1.2",
            opts="",
            sector_opts="8760-T-H-B-I-A-dist1",
            planning_horizons=2020,
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
        n, grouping_years_power, costs, baseyear
    )

    if options["heating"]:
        time_dep_hp_cop = options["time_dep_hp_cop"]
        ashp_cop = (
            xr.open_dataarray(snakemake.input.cop_air_total)
            .to_pandas()
            .reindex(index=n.snapshots)
        )
        gshp_cop = (
            xr.open_dataarray(snakemake.input.cop_soil_total)
            .to_pandas()
            .reindex(index=n.snapshots)
        )
        default_lifetime = snakemake.params.existing_capacities[
            "default_heating_lifetime"
        ]

        add_heating_capacities_installed_before_baseyear(
            n,
            baseyear,
            grouping_years_heat,
            ashp_cop,
            gshp_cop,
            time_dep_hp_cop,
            costs,
            default_lifetime,
        )

    if options.get("cluster_heat_buses", False):
        cluster_heat_buses(n)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))

    sanitize_carriers(n, snakemake.config)

    n.export_to_netcdf(snakemake.output[0])
