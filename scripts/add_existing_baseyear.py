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
import powerplantmatching as pm
import pypsa
import xarray as xr
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from add_electricity import sanitize_carriers
from prepare_sector_network import cluster_heat_buses, define_spatial, prepare_costs
from pypsa.io import import_components_from_dataframe

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


def add_existing_renewables(df_agg, costs):
    """
    Append existing renewables to the df_agg pd.DataFrame with the conventional
    power plants.
    """
    tech_map = {"solar": "PV", "onwind": "Onshore", "offwind-ac": "Offshore"}

    countries = snakemake.config["countries"]  # noqa: F841
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
        "Bioenergy": "solid biomass",
    }

    # If heat is considered, add CHPs in the add_heating_capacities function.
    # Assume that all oil power plants are not CHPs.
    if "H" in snakemake.wildcards.sector_opts.split("-"):
        df_agg = df_agg.query("Set != 'CHP'")
    elif (
        "I" not in snakemake.wildcards.sector_opts.split("-")
        and "Industry" in df_agg.columns
    ):
        df_agg["Industry"].fillna(False, inplace=True)
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

    # drop unused fueltypes and technologies
    df_agg.drop(df_agg.index[df_agg.Fueltype.isin(fueltype_to_drop)], inplace=True)
    df_agg.drop(df_agg.index[df_agg.Technology.isin(technology_to_drop)], inplace=True)
    df_agg.Fueltype = df_agg.Fueltype.map(rename_fuel)

    # Intermediate fix for DateIn & DateOut
    # Fill missing DateIn
    biomass_i = df_agg.loc[df_agg.Fueltype == "solid biomass"].index
    mean = df_agg.loc[biomass_i, "DateIn"].mean()
    df_agg.loc[biomass_i, "DateIn"] = df_agg.loc[biomass_i, "DateIn"].fillna(int(mean))
    # Fill missing DateOut
    dateout = (
        df_agg.loc[biomass_i, "DateIn"]
        + snakemake.params.costs["fill_values"]["lifetime"]
    )
    df_agg.loc[biomass_i, "DateOut"] = df_agg.loc[biomass_i, "DateOut"].fillna(dateout)

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
    add_existing_renewables(df_agg, costs)

    # add chp plants
    add_chp_plants(n, grouping_years, costs, baseyear, clustermaps)

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
        "solid biomass": "biomass",
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

            if "m" in snakemake.wildcards.clusters:
                for ind in new_capacity.index:
                    # existing capacities are split evenly among regions in every country
                    inv_ind = list(inv_busmap[ind])

                    # for offshore the splitting only includes coastal regions
                    inv_ind = [
                        i for i in inv_ind if (i + name_suffix_by) in n.generators.index
                    ]

                    p_max_pu = n.generators_t.p_max_pu[
                        [i + name_suffix_by for i in inv_ind]
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
                        efficiency=costs.at[cost_key, "efficiency"],
                        p_max_pu=p_max_pu,
                        build_year=grouping_year,
                        lifetime=costs.at[cost_key, "lifetime"],
                    )

            else:
                p_max_pu = n.generators_t.p_max_pu[capacity.index + name_suffix_by]

                if not new_build.empty:
                    n.madd(
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

                if generator != "solid biomass":
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
                        bus0=spatial.biomass.df.loc[new_capacity.index, "nodes"].values,
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


def add_chp_plants(n, grouping_years, costs, baseyear, clustermaps):
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

    # calculate remaining lifetime before phase-out (+1 because assuming
    # phase out date at the end of the year)
    chp.Fueltype = chp.Fueltype.map(rename_fuel)

    # assign clustered bus
    chp["bus"] = chp["bus"].astype(int)
    chp["cluster_bus"] = chp.bus.map(clustermaps)

    chp["grouping_year"] = np.take(
        grouping_years, np.digitize(chp.DateIn, grouping_years, right=True)
    )
    chp["lifetime"] = chp.DateOut - chp["grouping_year"] + 1

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

            efficiency_power = mastr_chp_efficiency_power.loc[grouping_year, generator]
            efficiency_heat = mastr_chp_efficiency_heat.loc[grouping_year, generator]

            if generator != "urban central solid biomass CHP":
                # lignite CHPs are not in DEA database - use coal CHP parameters
                key = keys[generator]
                if "EU" in vars(spatial)[generator].locations:
                    bus0 = vars(spatial)[generator].nodes
                else:
                    bus0 = vars(spatial)[generator].df.loc[p_nom.index, "nodes"]
                n.madd(
                    "Link",
                    p_nom.index,
                    suffix=f" urban central {generator} CHP-{grouping_year}",
                    bus0=bus0,
                    bus1=p_nom.index,
                    bus2=p_nom.index + " urban central heat",
                    bus3="co2 atmosphere",
                    carrier=f"urban central {generator} CHP",
                    p_nom=p_nom,
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
                    suffix=f" urban {key}-{grouping_year}",
                    bus0=spatial.biomass.df.loc[p_nom.index]["nodes"],
                    bus1=p_nom.index,
                    bus2=p_nom.index + " urban central heat",
                    carrier=generator,
                    p_nom=p_nom,
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

        if generator != "urban central solid biomass CHP":
            # lignite CHPs are not in DEA database - use coal CHP parameters
            key = keys[generator]
            if "EU" in vars(spatial)[generator].locations:
                bus0 = vars(spatial)[generator].nodes
            else:
                bus0 = vars(spatial)[generator].df.loc[p_nom.index, "nodes"]
            n.madd(
                "Link",
                p_nom.index,
                suffix=f" urban central {generator} CHP-{grouping_year}",
                bus0=bus0,
                bus1=p_nom.index,
                bus2=p_nom.index + " urban central heat",
                bus3="co2 atmosphere",
                carrier=f"urban central {generator} CHP",
                p_nom=p_nom / costs.at[key, "efficiency"],
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
                suffix=f" urban {key}-{grouping_year}",
                bus0=spatial.biomass.df.loc[p_nom.index]["nodes"],
                bus1=p_nom.index,
                bus2=p_nom.index + " urban central heat",
                carrier=generator,
                p_nom=p_nom / costs.at[key, "efficiency"],
                capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"],
                marginal_cost=costs.at[key, "VOM"],
                efficiency=costs.at[key, "efficiency"],
                efficiency2=costs.at[key, "efficiency-heat"],
                build_year=grouping_year,
                lifetime=costs.at[key, "lifetime"],
            )


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
        in rural buses 50% in urban buses
    """
    logger.debug(f"Adding heating capacities installed before {baseyear}")

    existing_heating = pd.read_csv(
        snakemake.input.existing_heating_distribution, header=[0, 1], index_col=0
    )

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

        too_large_grouping_years = [gy for gy in grouping_years if gy >= int(baseyear)]
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
                bus0=spatial.biomass.df.loc[nodes, "nodes"].values,
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


def add_h2_retro(n, baseyear, params):
    """
    Function to add H2 retrofitting of existing plants.
    """
    logger.info("Add H2 retrofitting.")
    plant_types = [
        ("OCGT", "retrofitted H2 OCGT"),
        ("CCGT", "retrofitted H2 CCGT"),
        ("urban central gas CHP", "urban central retrofitted H2 CHP"),
    ]
    start = params.retrofit_start

    for original_carrier, new_carrier in plant_types:
        # Query to filter the DataFrame
        plant_i = n.links.query(f"carrier == '{original_carrier}'")
        # Further filtering based on build_year excluding the current planning horizon
        plant_i = plant_i.loc[(plant_i.build_year >= start) & (plant_i.build_year != baseyear)].index
        # Set plants to extendable for constraint in solve_network()
        n.links.loc[plant_i, "p_nom_extendable"] = True

        # get all retrofitted plants and set p_nom_extendable to True
        h2_counterpart = n.links.query(f"carrier == '{new_carrier}'").index
        n.links.loc[h2_counterpart, "p_nom_extendable"] = True

        # check if the plants are not having a h2 counterpart yet
        h2_counterpart = set(
            index.replace(new_carrier, original_carrier) for index in h2_counterpart
        )
        plant_i = [index for index in plant_i if index not in h2_counterpart]

        if not plant_i:
            logger.info(f"No more {original_carrier} retrofitting potential.")
            continue

        n.links.loc[plant_i, "p_nom_max"] = n.links.loc[plant_i, "p_nom"]

        df = n.links.loc[plant_i].copy()
        # Adjust bus 0
        df["bus0"] = df.bus0.map(n.buses.location) + " H2"
        # Rename carrier and index
        df["carrier"] = df.carrier.apply(
            lambda x: x.replace(original_carrier, new_carrier)
        )
        df.rename(
            index=lambda x: x.replace(original_carrier, new_carrier),
            inplace=True,
        )
        df.loc[:, "capital_cost"] = params.retrofit_cost
        df.loc[:, "efficiency"] = params.retrofit_efficiency
        # Set p_nom_max to gas plant p_nom and existing capacity to zero
        df.loc[:, "p_nom"] = 0
        # Set CO2 emissions to 0: CHP plants have co2 emissions at bus3, gas plants at bus2
        if original_carrier != "urban central gas CHP":
            df.loc[:, "efficiency2"] = 0.0
        else:
            df.loc[:, "efficiency3"] = 0.0
        # Build_year and lifetime will stay the same as decommissioning of gas plant
        # Add retrofitted plant to network
        import_components_from_dataframe(n, df, "Link")


# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_existing_baseyear",
            configfiles="config/scenarios.automated.yaml",
            simpl="",
            clusters="22",
            ll="vopt",
            opts="",
            sector_opts="none",
            planning_horizons=2020,
            run="KN2045_Bal_v4",
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

    if snakemake.params.H2_retrofit_plants:
        # allow retrofitting of gas power plants to H2
        add_h2_retro(
            n,
            baseyear,
            snakemake.params,
        )

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))

    sanitize_carriers(n, snakemake.config)

    n.export_to_netcdf(snakemake.output[0])

# %%
