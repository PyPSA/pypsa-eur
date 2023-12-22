# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Concats pypsa networks of single investment periods to one network.
"""

import logging
import re

import numpy as np
import pandas as pd
import pypsa
from _helpers import update_config_with_sector_opts
from add_existing_baseyear import add_build_year_to_new_assets
from pypsa.descriptors import expand_series
from pypsa.io import import_components_from_dataframe
from six import iterkeys

logger = logging.getLogger(__name__)


# helper functions ---------------------------------------------------
def get_missing(df, n, c):
    """
    Get in network n missing assets of df for component c.

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
    """
    Calculate for a given time t and social discount rate r [per unit] the
    social discount.
    """
    return 1 / (1 + r) ** t


def get_investment_weighting(time_weighting, r=0.01):
    """
    Define cost weighting.

    Returns cost weightings depending on the the time_weighting
    (pd.Series) and the social discountrate r
    """
    end = time_weighting.cumsum()
    start = time_weighting.cumsum().shift().fillna(0)
    return pd.concat([start, end], axis=1).apply(
        lambda x: sum(get_social_discount(t, r) for t in range(int(x[0]), int(x[1]))),
        axis=1,
    )


def add_year_to_constraints(n, baseyear):
    """
    Add investment period to global constraints and rename index.

    Parameters
    ----------
    n : pypsa.Network
    baseyear : int
        year in which optimized assets are built
    """

    for c in n.iterate_components(["GlobalConstraint"]):
        c.df["investment_period"] = baseyear
        c.df.rename(index=lambda x: x + "-" + str(baseyear), inplace=True)


def hvdc_transport_model(n):
    """
    Convert AC lines to DC links for multi-decade optimisation with line
    expansion.

    Losses of DC links are assumed to be 3% per 1000km
    """

    logger.info("Convert AC lines to DC links to perform multi-decade optimisation.")

    n.madd(
        "Link",
        n.lines.index,
        bus0=n.lines.bus0,
        bus1=n.lines.bus1,
        p_nom_extendable=True,
        p_nom=n.lines.s_nom,
        p_nom_min=n.lines.s_nom,
        p_min_pu=-1,
        efficiency=1 - 0.03 * n.lines.length / 1000,
        marginal_cost=0,
        carrier="DC",
        length=n.lines.length,
        capital_cost=n.lines.capital_cost,
    )

    # Remove AC lines
    logger.info("Removing AC lines")
    lines_rm = n.lines.index
    n.mremove("Line", lines_rm)

    # Set efficiency of all DC links to include losses depending on length
    n.links.loc[n.links.carrier == "DC", "efficiency"] = (
        1 - 0.03 * n.links.loc[n.links.carrier == "DC", "length"] / 1000
    )


def adjust_electricity_grid(n, year, years):
    """
    Add carrier to lines. Replace AC lines with DC links in case of line
    expansion. Add lifetime to DC links in case of line expansion.

    Parameters
    ----------
    n    : pypsa.Network
    year : int
           year in which optimized assets are built
    years: list
           investment periods
    """
    n.lines["carrier"] = "AC"
    links_i = n.links[n.links.carrier == "DC"].index
    if n.lines.s_nom_extendable.any() or n.links.loc[links_i, "p_nom_extendable"].any():
        hvdc_transport_model(n)
        links_i = n.links[n.links.carrier == "DC"].index
        n.links.loc[links_i, "lifetime"] = 100
        if year != years[0]:
            n.links.loc[links_i, "p_nom_min"] = 0
            n.links.loc[links_i, "p_nom"] = 0


# --------------------------------------------------------------------
def concat_networks(years):
    """
    Concat given pypsa networks and adds build_year.

    Return:
        n : pypsa.Network for the whole planning horizon
    """

    # input paths of sector coupling networks
    network_paths = [snakemake.input.brownfield_network] + [
        snakemake.input[f"network_{year}"] for year in years[1:]
    ]
    # final concatenated network
    n = pypsa.Network()

    # iterate over single year networks and concat to perfect foresight network
    for i, network_path in enumerate(network_paths):
        year = years[i]
        network = pypsa.Network(network_path)
        adjust_electricity_grid(network, year, years)
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
                if pnl_year.empty and ~(component.name == "Load" and k == "p_set"):
                    continue
                if component.name == "Load":
                    static_load = network.loads.loc[network.loads.p_set != 0]
                    static_load_t = expand_series(static_load.p_set, network_sns).T
                    pnl_year = pd.concat(
                        [pnl_year.reindex(network_sns), static_load_t], axis=1
                    )
                    columns = (pnl[k].columns.union(pnl_year.columns)).unique()
                    pnl[k] = pnl[k].reindex(columns=columns)
                    pnl[k].loc[pnl_year.index, pnl_year.columns] = pnl_year

                else:
                    # this is to avoid adding multiple times assets with
                    # infinite lifetime as ror
                    cols = pnl_year.columns.difference(pnl[k].columns)
                    pnl[k] = pd.concat([pnl[k], pnl_year[cols]], axis=1)

        n.snapshot_weightings.loc[year, :] = network.snapshot_weightings.values

        # (3) global constraints
        for component in network.iterate_components(["GlobalConstraint"]):
            add_year_to_constraints(network, year)
            import_components_from_dataframe(n, component.df, component.name)

    # set investment periods
    n.investment_periods = n.snapshots.levels[0]
    # weighting of the investment period -> assuming last period same weighting as the period before
    time_w = n.investment_periods.to_series().diff().shift(-1).fillna(method="ffill")
    n.investment_period_weightings["years"] = time_w
    # set objective weightings
    objective_w = get_investment_weighting(
        n.investment_period_weightings["years"], social_discountrate
    )
    n.investment_period_weightings["objective"] = objective_w
    # all former static loads are now time-dependent -> set static = 0
    n.loads["p_set"] = 0
    n.loads_t.p_set.fillna(0, inplace=True)

    return n


def adjust_stores(n):
    """
    Make sure that stores still behave cyclic over one year and not whole
    modelling horizon.
    """
    # cyclic constraint
    cyclic_i = n.stores[n.stores.e_cyclic].index
    n.stores.loc[cyclic_i, "e_cyclic_per_period"] = True
    n.stores.loc[cyclic_i, "e_cyclic"] = False
    # non cyclic store assumptions
    non_cyclic_store = ["co2", "co2 stored", "solid biomass", "biogas", "Li ion"]
    co2_i = n.stores[n.stores.carrier.isin(non_cyclic_store)].index
    n.stores.loc[co2_i, "e_cyclic_per_period"] = False
    n.stores.loc[co2_i, "e_cyclic"] = False
    # e_initial at beginning of each investment period
    e_initial_store = ["solid biomass", "biogas"]
    co2_i = n.stores[n.stores.carrier.isin(e_initial_store)].index
    n.stores.loc[co2_i, "e_initial_per_period"] = True
    # n.stores.loc[co2_i, "e_initial"] *= 10
    # n.stores.loc[co2_i, "e_nom"] *= 10
    e_initial_store = ["co2 stored"]
    co2_i = n.stores[n.stores.carrier.isin(e_initial_store)].index
    n.stores.loc[co2_i, "e_initial_per_period"] = True

    return n


def set_phase_out(n, carrier, ct, phase_out_year):
    """
    Set planned phase outs for given carrier,country (ct) and planned year of
    phase out (phase_out_year).
    """
    df = n.links[(n.links.carrier.isin(carrier)) & (n.links.bus1.str[:2] == ct)]
    # assets which are going to be phased out before end of their lifetime
    assets_i = df[df[["build_year", "lifetime"]].sum(axis=1) > phase_out_year].index
    build_year = n.links.loc[assets_i, "build_year"]
    # adjust lifetime
    n.links.loc[assets_i, "lifetime"] = (phase_out_year - build_year).astype(float)


def set_all_phase_outs(n):
    # TODO move this to a csv or to the config
    planned = [
        (["nuclear"], "DE", 2022),
        (["nuclear"], "BE", 2025),
        (["nuclear"], "ES", 2027),
        (["coal", "lignite"], "DE", 2030),
        (["coal", "lignite"], "ES", 2027),
        (["coal", "lignite"], "FR", 2022),
        (["coal", "lignite"], "GB", 2024),
        (["coal", "lignite"], "IT", 2025),
        (["coal", "lignite"], "DK", 2030),
        (["coal", "lignite"], "FI", 2030),
        (["coal", "lignite"], "HU", 2030),
        (["coal", "lignite"], "SK", 2030),
        (["coal", "lignite"], "GR", 2030),
        (["coal", "lignite"], "IE", 2030),
        (["coal", "lignite"], "NL", 2030),
        (["coal", "lignite"], "RS", 2030),
    ]
    for carrier, ct, phase_out_year in planned:
        set_phase_out(n, carrier, ct, phase_out_year)
    # remove assets which are already phased out
    remove_i = n.links[n.links[["build_year", "lifetime"]].sum(axis=1) < years[0]].index
    n.mremove("Link", remove_i)


def set_carbon_constraints(n, opts):
    """
    Add global constraints for carbon emissions.
    """
    budget = None
    for o in opts:
        # other budgets
        m = re.match(r"^\d+p\d$", o, re.IGNORECASE)
        if m is not None:
            budget = snakemake.config["co2_budget"][m.group(0)] * 1e9
    if budget != None:
        logger.info(f"add carbon budget of {budget}")
        n.add(
            "GlobalConstraint",
            "Budget",
            type="Co2Budget",
            carrier_attribute="co2_emissions",
            sense="<=",
            constant=budget,
            investment_period=n.investment_periods[-1],
        )

        # drop other CO2 limits
        drop_i = n.global_constraints[n.global_constraints.type == "co2_limit"].index
        n.mremove("GlobalConstraint", drop_i)

        n.add(
            "GlobalConstraint",
            "carbon_neutral",
            type="co2_limit",
            carrier_attribute="co2_emissions",
            sense="<=",
            constant=0,
            investment_period=n.investment_periods[-1],
        )

    # set minimum CO2 emission constraint to avoid too fast reduction
    if "co2min" in opts:
        emissions_1990 = 4.53693
        emissions_2019 = 3.344096
        target_2030 = 0.45 * emissions_1990
        annual_reduction = (emissions_2019 - target_2030) / 11
        first_year = n.snapshots.levels[0][0]
        time_weightings = n.investment_period_weightings.loc[first_year, "years"]
        co2min = emissions_2019 - ((first_year - 2019) * annual_reduction)
        logger.info(f"add minimum emissions for {first_year} of {co2min} t CO2/a")
        n.add(
            "GlobalConstraint",
            f"Co2Min-{first_year}",
            type="Co2min",
            carrier_attribute="co2_emissions",
            sense=">=",
            investment_period=first_year,
            constant=co2min * 1e9 * time_weightings,
        )

    return n


def adjust_lvlimit(n):
    """
    Convert global constraints for single investment period to one uniform if
    all attributes stay the same.
    """
    c = "GlobalConstraint"
    cols = ["carrier_attribute", "sense", "constant", "type"]
    glc_type = "transmission_volume_expansion_limit"
    if (n.df(c)[n.df(c).type == glc_type][cols].nunique() == 1).all():
        glc = n.df(c)[n.df(c).type == glc_type][cols].iloc[[0]]
        glc.index = pd.Index(["lv_limit"])
        remove_i = n.df(c)[n.df(c).type == glc_type].index
        n.mremove(c, remove_i)
        import_components_from_dataframe(n, glc, c)

    return n


def adjust_CO2_glc(n):
    c = "GlobalConstraint"
    glc_name = "CO2Limit"
    glc_type = "primary_energy"
    mask = (n.df(c).index.str.contains(glc_name)) & (n.df(c).type == glc_type)
    n.df(c).loc[mask, "type"] = "co2_limit"

    return n


def add_H2_boilers(n):
    """
    Gas boilers can be retrofitted to run with H2.

    Add H2 boilers for heating for all existing gas boilers.
    """
    c = "Link"
    logger.info("Add H2 boilers.")
    # existing gas boilers
    mask = n.links.carrier.str.contains("gas boiler") & ~n.links.p_nom_extendable
    gas_i = n.links[mask].index
    df = n.links.loc[gas_i]
    # adjust bus 0
    df["bus0"] = df.bus1.map(n.buses.location) + " H2"
    # rename carrier and index
    df["carrier"] = df.carrier.apply(
        lambda x: x.replace("gas boiler", "retrofitted H2 boiler")
    )
    df.rename(
        index=lambda x: x.replace("gas boiler", "retrofitted H2 boiler"), inplace=True
    )
    # todo, costs for retrofitting
    df["capital_costs"] = 100
    # set existing capacity to zero
    df["p_nom"] = 0
    df["p_nom_extendable"] = True
    # add H2 boilers to network
    import_components_from_dataframe(n, df, c)


def apply_time_segmentation_perfect(
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
    raw = raw.dropna(axis=1)
    sn_weightings = {}

    for year in raw.index.levels[0]:
        logger.info(f"Find representative snapshots for {year}.")
        raw_t = raw.loc[year]
        # normalise all time-dependent data
        annual_max = raw_t.max().replace(0, 1)
        raw_t = raw_t.div(annual_max, level=0)
        # get representative segments
        agg = tsam.TimeSeriesAggregation(
            raw_t,
            hoursPerPeriod=len(raw_t),
            noTypicalPeriods=1,
            noSegments=int(segments),
            segmentation=True,
            solver=solver_name,
        )
        segmented = agg.createTypicalPeriods()

        weightings = segmented.index.get_level_values("Segment Duration")
        offsets = np.insert(np.cumsum(weightings[:-1]), 0, 0)
        timesteps = [raw_t.index[0] + pd.Timedelta(f"{offset}h") for offset in offsets]
        snapshots = pd.DatetimeIndex(timesteps)
        sn_weightings[year] = pd.Series(
            weightings, index=snapshots, name="weightings", dtype="float64"
        )

    sn_weightings = pd.concat(sn_weightings)
    n.set_snapshots(sn_weightings.index)
    n.snapshot_weightings = n.snapshot_weightings.mul(sn_weightings, axis=0)

    return n


def set_temporal_aggregation_SEG(n, opts, solver_name):
    """
    Aggregate network temporally with tsam.
    """
    for o in opts:
        # segments with package tsam
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            segments = int(m[1])
            logger.info(f"Use temporal segmentation with {segments} segments")
            n = apply_time_segmentation_perfect(n, segments, solver_name=solver_name)
            break
    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_perfect_foresight",
            simpl="",
            opts="",
            clusters="37",
            ll="v1.5",
            sector_opts="1p7-4380H-T-H-B-I-A-solar+p3-dist1",
        )

    update_config_with_sector_opts(snakemake.config, snakemake.wildcards.sector_opts)
    # parameters -----------------------------------------------------------
    years = snakemake.config["scenario"]["planning_horizons"]
    opts = snakemake.wildcards.sector_opts.split("-")
    social_discountrate = snakemake.config["costs"]["social_discountrate"]
    for o in opts:
        if "sdr" in o:
            social_discountrate = float(o.replace("sdr", "")) / 100

    logger.info(
        f"Concat networks of investment period {years} with social discount rate of {social_discountrate * 100}%"
    )

    # concat prenetworks of planning horizon to single network ------------
    n = concat_networks(years)

    # temporal aggregate
    opts = snakemake.wildcards.sector_opts.split("-")
    solver_name = snakemake.config["solving"]["solver"]["name"]
    n = set_temporal_aggregation_SEG(n, opts, solver_name)

    # adjust global constraints lv limit if the same for all years
    n = adjust_lvlimit(n)
    # adjust global constraints CO2 limit
    n = adjust_CO2_glc(n)
    # adjust stores to multi period investment
    n = adjust_stores(n)

    # set phase outs
    set_all_phase_outs(n)

    # add H2 boiler
    add_H2_boilers(n)

    # set carbon constraints
    opts = snakemake.wildcards.sector_opts.split("-")
    n = set_carbon_constraints(n, opts)

    # export network
    n.export_to_netcdf(snakemake.output[0])
