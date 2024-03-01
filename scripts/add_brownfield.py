# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Prepares brownfield data from previous planning horizon.
"""

import logging

logger = logging.getLogger(__name__)

import pandas as pd

idx = pd.IndexSlice

import numpy as np
import pypsa
from _helpers import update_config_with_sector_opts
from add_existing_baseyear import add_build_year_to_new_assets


def add_brownfield(n, n_p, year):
    logger.info(f"Preparing brownfield for the year {year}")

    # electric transmission grid set optimised capacities of previous as minimum
    n.lines.s_nom_min = n_p.lines.s_nom_opt
    dc_i = n.links[n.links.carrier == "DC"].index
    n.links.loc[dc_i, "p_nom_min"] = n_p.links.loc[dc_i, "p_nom_opt"]

    for c in n_p.iterate_components(["Link", "Generator", "Store"]):
        attr = "e" if c.name == "Store" else "p"

        # first, remove generators, links and stores that track
        # CO2 or global EU values since these are already in n
        n_p.mremove(c.name, c.df.index[c.df.lifetime == np.inf])

        # remove assets whose build_year + lifetime < year
        n_p.mremove(c.name, c.df.index[c.df.build_year + c.df.lifetime < year])

        # remove assets if their optimized nominal capacity is lower than a threshold
        # since CHP heat Link is proportional to CHP electric Link, make sure threshold is compatible
        chp_heat = c.df.index[
            (c.df[f"{attr}_nom_extendable"] & c.df.index.str.contains("urban central"))
            & c.df.index.str.contains("CHP")
            & c.df.index.str.contains("heat")
        ]

        threshold = snakemake.params.threshold_capacity

        if not chp_heat.empty:
            threshold_chp_heat = (
                threshold
                * c.df.efficiency[chp_heat.str.replace("heat", "electric")].values
                * c.df.p_nom_ratio[chp_heat.str.replace("heat", "electric")].values
                / c.df.efficiency[chp_heat].values
            )
            n_p.mremove(
                c.name,
                chp_heat[c.df.loc[chp_heat, f"{attr}_nom_opt"] < threshold_chp_heat],
            )

        n_p.mremove(
            c.name,
            c.df.index[
                (c.df[f"{attr}_nom_extendable"] & ~c.df.index.isin(chp_heat))
                & (c.df[f"{attr}_nom_opt"] < threshold)
            ],
        )
        check_transport = ( (snakemake.config["sector"]["land_transport_electric_share"][year] is not None ) 
        or (snakemake.config["sector"]["land_transport_fuel_cell_share"][year] is not None)
        or (snakemake.config["sector"]["land_transport_ice_share"][year] is not None))
        
        if not snakemake.config["sector"]["endogenous_transport"] or check_transport:
            n_p.mremove(
                c.name,
                c.df.index[c.df.carrier.str.contains("land transport" or "V2G" or "EV battery storage" or "BEV charger")
                ],
            )

        # copy over assets but fix their capacity
        c.df[f"{attr}_nom"] = c.df[f"{attr}_nom_opt"]
        c.df[f"{attr}_nom_extendable"] = False

        n.import_components_from_dataframe(c.df, c.name)

        # copy time-dependent
        selection = n.component_attrs[c.name].type.str.contains(
            "series"
        ) & n.component_attrs[c.name].status.str.contains("Input")
        for tattr in n.component_attrs[c.name].index[selection]:
            n.import_series_from_dataframe(c.pnl[tattr], c.name, tattr)

        # deal with gas network
        pipe_carrier = ["gas pipeline"]
        if snakemake.params.H2_retrofit:
            # drop capacities of previous year to avoid duplicating
            to_drop = n.links.carrier.isin(pipe_carrier) & (n.links.build_year != year)
            n.mremove("Link", n.links.loc[to_drop].index)

            # subtract the already retrofitted from today's gas grid capacity
            h2_retrofitted_fixed_i = n.links[
                (n.links.carrier == "H2 pipeline retrofitted")
                & (n.links.build_year != year)
            ].index
            gas_pipes_i = n.links[n.links.carrier.isin(pipe_carrier)].index
            CH4_per_H2 = 1 / snakemake.params.H2_retrofit_capacity_per_CH4
            fr = "H2 pipeline retrofitted"
            to = "gas pipeline"
            # today's pipe capacity
            pipe_capacity = n.links.loc[gas_pipes_i, "p_nom"]
            # already retrofitted capacity from gas -> H2
            already_retrofitted = (
                n.links.loc[h2_retrofitted_fixed_i, "p_nom"]
                .rename(lambda x: x.split("-2")[0].replace(fr, to))
                .groupby(level=0)
                .sum()
            )
            remaining_capacity = (
                pipe_capacity
                - CH4_per_H2
                * already_retrofitted.reindex(index=pipe_capacity.index).fillna(0)
            )
            n.links.loc[gas_pipes_i, "p_nom"] = remaining_capacity
        else:
            new_pipes = n.links.carrier.isin(pipe_carrier) & (
                n.links.build_year == year
            )
            n.links.loc[new_pipes, "p_nom"] = 0.0
            n.links.loc[new_pipes, "p_nom_min"] = 0.0

def adjust_EVs(n, n_p, year):
    # set p_min_pu and p_max_pu for solved network, so that only the EV link for the current time horizon
    # is not constraint (constraining all land transport link with p_min_pu/p_max_pu while p_nom_extentable=True leads to infeasible by gurobi)
    lifetime_EV = n.links.lifetime[n.links[(n.links.carrier=="land transport EV") ].index[0]]
    i = 0
    while (year-lifetime_EV+i) < year:
        if not n.links_t.efficiency[(n.links.filter(like="land transport EV-"+str(int(year-lifetime_EV+i)),axis=0)).index].empty:
            p_set =  n_p.loads_t.p[n.loads[n.loads.carrier.str.contains('land transport demand')].index]
            eff = n.links_t.efficiency[(n.links.filter(like="land transport EV-"+str(int(year-lifetime_EV+i)),axis=0)).index]
            p_set = p_set.add_suffix(' EV-'+str(int(year-lifetime_EV+i)))
            p_set = p_set.drop([col for col in p_set.columns if col in p_set.columns and col not in eff.columns], axis=1)
            pnom = (p_set.divide(eff)).max()
            pu=p_set.divide(eff)/pnom
            n.links_t.p_min_pu[(n.links.filter(like="land transport EV-"+str(int(year-lifetime_EV+i)),axis=0)).index] = pu.values
            n.links_t.p_max_pu[(n.links.filter(like="land transport EV-"+str(int(year-lifetime_EV+i)),axis=0)).index] = pu.values
        i = i+1   

def disable_grid_expansion_if_LV_limit_hit(n):
    if not "lv_limit" in n.global_constraints.index:
        return

    total_expansion = (
        n.lines.eval("s_nom_min * length").sum()
        + n.links.query("carrier == 'DC'").eval("p_nom_min * length").sum()
    ).sum()

    lv_limit = n.global_constraints.at["lv_limit", "constant"]

    # allow small numerical differences
    if lv_limit - total_expansion < 1:
        logger.info(f"LV is already reached, disabling expansion and LV limit")
        extendable_acs = n.lines.query("s_nom_extendable").index
        n.lines.loc[extendable_acs, "s_nom_extendable"] = False
        n.lines.loc[extendable_acs, "s_nom"] = n.lines.loc[extendable_acs, "s_nom_min"]

        extendable_dcs = n.links.query("carrier == 'DC' and p_nom_extendable").index
        n.links.loc[extendable_dcs, "p_nom_extendable"] = False
        n.links.loc[extendable_dcs, "p_nom"] = n.links.loc[extendable_dcs, "p_nom_min"]

        n.global_constraints.drop("lv_limit", inplace=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_brownfield",
            simpl="",
            clusters="37",
            opts="",
            ll="v1.0",
            sector_opts="168H-T-H-B-I-solar+p3-dist1",
            planning_horizons=2030,
        )

    logging.basicConfig(level=snakemake.config["logging"]["level"])

    update_config_with_sector_opts(snakemake.config, snakemake.wildcards.sector_opts)

    logger.info(f"Preparing brownfield from the file {snakemake.input.network_p}")

    year = int(snakemake.wildcards.planning_horizons)

    n = pypsa.Network(snakemake.input.network)

    add_build_year_to_new_assets(n, year)

    n_p = pypsa.Network(snakemake.input.network_p)

    add_brownfield(n, n_p, year)

    disable_grid_expansion_if_LV_limit_hit(n)

    opts = snakemake.wildcards.sector_opts.split("-")
    if "T" in opts and snakemake.config["sector"]["endogenous_transport"]:
        adjust_EVs(n, n_p, year)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
