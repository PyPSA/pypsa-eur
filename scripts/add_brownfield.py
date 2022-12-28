# coding: utf-8

import logging
logger = logging.getLogger(__name__)

import pandas as pd
idx = pd.IndexSlice

import pypsa
import yaml
import numpy as np

from add_existing_baseyear import add_build_year_to_new_assets
from helper import override_component_attrs, update_config_with_sector_opts
from solve_network import basename


def add_brownfield(n, n_p, year):

    print("adding brownfield")

    # electric transmission grid set optimised capacities of previous as minimum
    n.lines.s_nom_min = n_p.lines.s_nom_opt
    dc_i = n.links[n.links.carrier=="DC"].index
    n.links.loc[dc_i, "p_nom_min"] = n_p.links.loc[dc_i, "p_nom_opt"]

    for c in n_p.iterate_components(["Link", "Generator", "Store"]):

        attr = "e" if c.name == "Store" else "p"

        # first, remove generators, links and stores that track
        # CO2 or global EU values since these are already in n
        n_p.mremove(
            c.name,
            c.df.index[c.df.lifetime==np.inf]
        )

        # remove assets whose build_year + lifetime < year
        n_p.mremove(
            c.name,
            c.df.index[c.df.build_year + c.df.lifetime < year]
        )

        # remove assets if their optimized nominal capacity is lower than a threshold
        # since CHP heat Link is proportional to CHP electric Link, make sure threshold is compatible
        chp_heat = c.df.index[(
            c.df[attr + "_nom_extendable"]
            & c.df.index.str.contains("urban central")
            & c.df.index.str.contains("CHP")
            & c.df.index.str.contains("heat")
        )]

        threshold = snakemake.config['existing_capacities']['threshold_capacity']

        if not chp_heat.empty:
            threshold_chp_heat = (threshold
                * c.df.efficiency[chp_heat.str.replace("heat", "electric")].values
                * c.df.p_nom_ratio[chp_heat.str.replace("heat", "electric")].values
                / c.df.efficiency[chp_heat].values
            )
            n_p.mremove(
                c.name,
                chp_heat[c.df.loc[chp_heat, attr + "_nom_opt"] < threshold_chp_heat]
            )

        n_p.mremove(
            c.name,
            c.df.index[c.df[attr + "_nom_extendable"] & ~c.df.index.isin(chp_heat) & (c.df[attr + "_nom_opt"] < threshold)]
        )

        # copy over assets but fix their capacity
        c.df[attr + "_nom"] = c.df[attr + "_nom_opt"]
        c.df[attr + "_nom_extendable"] = False

        n.import_components_from_dataframe(c.df, c.name)

        # copy time-dependent
        selection = (
            n.component_attrs[c.name].type.str.contains("series")
            & n.component_attrs[c.name].status.str.contains("Input")
        )
        for tattr in n.component_attrs[c.name].index[selection]:
            n.import_series_from_dataframe(c.pnl[tattr], c.name, tattr)

        # deal with gas network
        pipe_carrier = ['gas pipeline']
        if snakemake.config["sector"]['H2_retrofit']:
            # drop capacities of previous year to avoid duplicating
            to_drop = n.links.carrier.isin(pipe_carrier) & (n.links.build_year!=year)
            n.mremove("Link", n.links.loc[to_drop].index)

            # subtract the already retrofitted from today's gas grid capacity
            h2_retrofitted_fixed_i = n.links[(n.links.carrier=='H2 pipeline retrofitted') & (n.links.build_year!=year)].index
            gas_pipes_i =  n.links[n.links.carrier.isin(pipe_carrier)].index
            CH4_per_H2 = 1 / snakemake.config["sector"]["H2_retrofit_capacity_per_CH4"]
            fr = "H2 pipeline retrofitted"
            to = "gas pipeline"
            # today's pipe capacity
            pipe_capacity = n.links.loc[gas_pipes_i, 'p_nom']
            # already retrofitted capacity from gas -> H2
            already_retrofitted = (n.links.loc[h2_retrofitted_fixed_i, 'p_nom']
                                   .rename(lambda x: basename(x).replace(fr, to)).groupby(level=0).sum())
            remaining_capacity = pipe_capacity - CH4_per_H2 * already_retrofitted.reindex(index=pipe_capacity.index).fillna(0)
            n.links.loc[gas_pipes_i, "p_nom"] = remaining_capacity
        else:
            new_pipes = n.links.carrier.isin(pipe_carrier) & (n.links.build_year==year)
            n.links.loc[new_pipes, "p_nom"] = 0.
            n.links.loc[new_pipes, "p_nom_min"] = 0.



#%%
if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'add_brownfield',
            simpl='',
            clusters="37",
            opts="",
            lv=1.0,
            sector_opts='168H-T-H-B-I-solar+p3-dist1',
            planning_horizons=2030,
        )

    update_config_with_sector_opts(snakemake.config, snakemake.wildcards.sector_opts)

    print(snakemake.input.network_p)
    logging.basicConfig(level=snakemake.config['logging_level'])

    year = int(snakemake.wildcards.planning_horizons)

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    add_build_year_to_new_assets(n, year)

    n_p = pypsa.Network(snakemake.input.network_p, override_component_attrs=overrides)

    add_brownfield(n, n_p, year)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
