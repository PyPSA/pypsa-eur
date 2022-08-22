#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Concats pypsa networks of single investment periods to one network.

Created on Tue Aug 16 10:40:41 2022

@author: lisa
"""
import pypsa
import pandas as pd
from helper import override_component_attrs, update_config_with_sector_opts
from pypsa.io import import_components_from_dataframe
from add_existing_baseyear import add_build_year_to_new_assets
from six import iterkeys
from pypsa.descriptors import expand_series
import re
import logging
logger = logging.getLogger(__name__)

# helper functions ---------------------------------------------------
def get_missing(df, n, c):
    """Get in network n missing assets of df for component c.

        Input:
            df: pandas DataFrame, static values of pypsa components
            n : pypsa Network to which new assets should be added
            c : string, pypsa component.list_name (e.g. "generators")
        Return:
            pd.DataFrame with static values of missing assets
        """
    df_final = getattr(n, c)
    missing_i = df.index.difference(df_final.index)
    return df.loc[missing_i]


def get_social_discount(t, r=0.01):
    """Calculate for a given time t the social discount."""
    return 1 / (1 + r) ** t


def get_investment_weighting(time_weighting, r=0.01):
    """Define cost weighting.

    Returns cost weightings depending on the the time_weighting (pd.Series)
    and the social discountrate r
    """
    end = time_weighting.cumsum()
    start = time_weighting.cumsum().shift().fillna(0)
    return pd.concat([start, end], axis=1).apply(
        lambda x: sum([get_social_discount(t, r) for t in range(int(x[0]), int(x[1]))]),
        axis=1,
    )

# --------------------------------------------------------------------
def concat_networks(years):
    """Concat given pypsa networks and adds build_year.

        Return:
            n : pypsa.Network for the whole planning horizon

        """

    # input paths of sector coupling networks
    network_paths = [snakemake.input.brownfield_network] + snakemake.input.network[1:]
    # final concatenated network
    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network( override_component_attrs=overrides)


    # iterate over single year networks and concat to perfect foresight network
    for i, network_path in enumerate(network_paths):
        year = years[i]
        network = pypsa.Network(network_path, override_component_attrs=overrides)
        network.lines["carrier"] = "AC"
        add_build_year_to_new_assets(network, year)

        # static ----------------------------------
        # (1) add buses and carriers
        for component in network.iterate_components(["Bus", "Carrier"]):
            df_year = component.df
            # get missing assets
            missing = get_missing(df_year, n, component.list_name)
            import_components_from_dataframe(n, missing, component.name)
        # (2) add generators, links, stores and loads
        for component in network.iterate_components(
            ["Generator", "Link", "Store", "Load", "Line", "StorageUnit"]
        ):

            df_year = component.df.copy()
            missing = get_missing(df_year, n, component.list_name)

            import_components_from_dataframe(n, missing, component.name)

        # time variant --------------------------------------------------
        network_sns = pd.MultiIndex.from_product([[year], network.snapshots])
        snapshots = n.snapshots.drop("now", errors="ignore").union(network_sns)
        n.set_snapshots(snapshots)

        for component in network.iterate_components():
            pnl = getattr(n, component.list_name + "_t")
            for k in iterkeys(component.pnl):
                pnl_year = component.pnl[k].copy().reindex(snapshots, level=1)
                if pnl_year.empty and ~(component.name=="Load" and k=="p_set"): continue
                if component.name == "Load":
                    static_load = network.loads.loc[network.loads.p_set != 0]
                    static_load_t = expand_series(
                        static_load.p_set, network_sns
                    ).T
                    pnl_year = pd.concat([pnl_year.reindex(network_sns),
                                          static_load_t], axis=1)
                    columns = (pnl[k].columns.union(pnl_year.columns)).unique()
                    pnl[k] = pnl[k].reindex(columns=columns)
                    pnl[k].loc[pnl_year.index, pnl_year.columns] = pnl_year

                else:
                    # this is to avoid adding multiple times assets with infinit lifetime as ror
                    cols = pnl_year.columns.difference(pnl[k].columns)
                    pnl[k] = pd.concat([pnl[k], pnl_year[cols]], axis=1)


        n.snapshot_weightings.loc[year,:] = network.snapshot_weightings.values
    # set investment periods
    n.investment_periods = n.snapshots.levels[0]
    # weighting of the investment period -> assuming last period same weighting as the period before
    time_w = n.investment_periods.to_series().diff().shift(-1).fillna(method="ffill")
    n.investment_period_weightings["years"] = time_w
    # set objective weightings
    objective_w = get_investment_weighting(n.investment_period_weightings["years"],
                                           social_discountrate)
    n.investment_period_weightings["objective"] = objective_w
    # all former static loads are now time-dependent -> set static = 0
    n.loads["p_set"] = 0

    return n

def adjust_stores(n):
    # cylclic constraint
    cyclic_i = n.stores[n.stores.e_cyclic].index
    n.stores.loc[cyclic_i, "e_cyclic_per_period"] = True
    n.stores.loc[cyclic_i, "e_cyclic"] = False
    # co2 store assumptions
    co2_i = n.stores[n.stores.carrier.isin(["co2", "co2 stored"])].index
    n.stores.loc[co2_i, "e_cyclic_per_period"] = False

    return n

def set_phase_out(n, carrier, ct, phase_out_year):
    df = n.links[(n.links.carrier.isin(carrier))& (n.links.bus1.str[:2]==ct)]
    assets_i = df[df[["build_year", "lifetime"]].sum(axis=1) > phase_out_year].index
    n.links.loc[assets_i, "lifetime"] = (phase_out_year - n.links.loc[assets_i, "build_year"]).astype(float)

def set_all_phase_outs(n):
    planned= [(["nuclear"], "DE", 2022),
              (["nuclear"], "BE", 2025),
              (["nuclear"], "ES", 2027),
              (["coal", "lignite"], "DE", 2038),
              (["coal", "lignite"], "ES", 2027),
              (["coal", "lignite"], "FR", 2022),
              (["coal", "lignite"], "GB", 2024),
              (["coal", "lignite"], "IT", 2025)]
    for carrier, ct, phase_out_year in planned:
        set_phase_out(n, carrier,ct, phase_out_year)
    # remove assets which are already phased out
    remove_i = n.links[n.links[["build_year", "lifetime"]].sum(axis=1)<years[0]].index
    n.mremove("Link", remove_i)


def set_carbon_constraints(n, opts):
    """Add global constraints for carbon emissions."""
    budget = (
        snakemake.config["co2_budget"]["1p7"] * 1e9
    )  # budget for + 1.7 Celsius for Europe
    for o in opts:
        # temporal clustering
        m = re.match(r"^\d+p\d$", o, re.IGNORECASE)
        if m is not None:
            budget = snakemake.config["co2_budget"][m.group(0)] * 1e9
    logger.info("add carbon budget of {}".format(budget))
    n.add(
        "GlobalConstraint",
        "Budget",
        type="Co2constraint",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=budget,
    )

    if not "noco2neutral" in opts:
        logger.info("Add carbon neutrality constraint.")
        n.add(
            "GlobalConstraint",
            "Co2neutral",
            type="Co2constraint",
            carrier_attribute="co2_emissions",
            investment_period=n.snapshots.levels[0][-1],
            sense="<=",
            constant=0,
        )

    if "co2min" in opts:
        emissions_1990 = 4.53693
        emissions_2019 = 3.344096
        target_2030 = 0.45*emissions_1990
        annual_reduction = (emissions_2019-target_2030)/11
        first_year = n.snapshots.levels[0][0]
        time_weightings = n.investment_period_weightings.loc[first_year, "years"]
        co2min = emissions_2019-((first_year-2019)*annual_reduction)
        logger.info("add minimum emissions for {} of {} t CO2/a".format(first_year, co2min))
        n.add(
            "GlobalConstraint",
            "Co2min",
            type="Co2constraint",
            carrier_attribute="co2_emissions",
            sense=">=",
            investment_period=first_year,
            constant=co2min*1e9*time_weightings,
        )

    return n
#%%
if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'prepare_perfect_foresight',
            simpl='',
            opts="",
            clusters="45",
            lv=1.0,
            sector_opts='365H-T-H-B-I-A-solar+p3-dist1-co2min',
        )

    update_config_with_sector_opts(snakemake.config, snakemake.wildcards.sector_opts)
    # parameters -----------------------------------------------------------
    years = snakemake.config["scenario"]["planning_horizons"]
    social_discountrate = snakemake.config["costs"]["social_discountrate"]
    logger.info("Concat networks of investment period {} with social discount rate of {}%"
                .format(years, social_discountrate*100))

    # concat prenetworks of planning horizon to single network ------------
    n = concat_networks(years)

    # set phase outs
    set_all_phase_outs(n)
    # adjust stores to multi period investment
    n = adjust_stores(n)

    # set carbon constraints
    opts = snakemake.wildcards.sector_opts.split('-')
    n = set_carbon_constraints(n, opts)

    # export network
    n.export_to_netcdf(snakemake.output[0])
