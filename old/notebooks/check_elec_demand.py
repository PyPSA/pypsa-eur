# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:40:10 2024

@author: alice
"""

import pypsa

# Retrieve the new network

n = pypsa.Network("../results/baseline/postnetworks/base_s_39_lvopt___2050.nc")
timestep = n.snapshot_weightings.iloc[0, 0]  # HOURS TO MULTIPLY FOR THE ANNUAL VALUES

capex = n.statistics.capex()
installed_capex = n.statistics.installed_capacity()
expanded_capex = n.statistics.expanded_capacity()

revenue = n.statistics.revenue()

supply = n.statistics.supply()
withdrawal = n.statistics.withdrawal()

# %%

# Match supply and withdrawal of electricity

withd = (
    n.statistics.withdrawal()
)  # every bus which is positive  in a link (entering the link)
supply = (
    n.statistics.supply()
)  # every bus which is negative in a link (exiting the link)
nuke = n.links[n.links.index.str.contains("nuclear")]

# Exogenous electrical demand
with_load = n.statistics.withdrawal(comps="Load")
with_load = (
    with_load[with_load.index.str.contains("electricity|EV", case=False)].sum() / 1e6
)  # TWh

# Endogenous electrical demand for HEATING -> they only have a negative output whi
with_th_links = n.statistics.withdrawal(comps="Link")
with_th_links = (
    with_th_links[with_th_links.index.str.contains("thermal|heat", case=False)].sum()
    / 1e6
)  # TWh

with_th_links = (
    n.links_t.p0.loc[
        :, n.links[(n.links.index.str.contains("thermal|heat", case=False))].index
    ]
    .sum()
    .sum()
    * timestep
    / 1e6
)  # TWh

# Endogenous electrical demand for METHANOLISATIONS
with_meth_links = (
    n.links_t.p2.loc[
        :,
        n.links[
            (n.links.bus2.str.contains(r"\d$", regex=True)) & (n.links.efficiency2 < 0)
        ].index,
    ]
    .sum()
    .sum()
    * timestep
    / 1e6
)  # TWh

# Endogenous electrical demand for DAC
with_dac_links = (
    n.links_t.p1.loc[:, n.links[(n.links.carrier.str.contains("DAC"))].index]
    .sum()
    .sum()
    * timestep
    / 1e6
)  # TWh

tot_dem = with_load + with_dac_links + with_meth_links + with_th_links

# Power supply
supply_gen = n.statistics.supply(comps="Generator")
supply_gen = (
    supply_gen[~supply_gen.index.str.contains("thermal|oil|gas", case=False)].sum()
    / 1e6
)  # TWh

supply_links = n.statistics.supply(comps="Link")
supply_links = (
    supply_links[
        supply_links.index.str.contains(
            "Cycle|Open|nuclear|coal| oil|lignite", case=False
        )
    ].sum()
    / 1e6
)  # TWh

supply_chp_links = (
    -n.links_t.p1.loc[
        :, n.links[(n.links.carrier.str.contains("CHP", case=False))].index
    ]
    .sum()
    .sum()
    * timestep
    / 1e6
)  # TWh

tot_sup = supply_gen + supply_links + supply_chp_links

curt = tot_sup - tot_dem
