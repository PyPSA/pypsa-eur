# SPDX-FileCopyrightText: : 2023- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


def custom_extra_functionality(n, snapshots, snakemake):
    """
    Add custom extra functionality constraints.
    """
    devices = snakemake.config["load"]["demand_response"]["devices"]
    for dev, cfg in devices.items(): 
        S = cfg["shift_time"]
        DR_i = n.stores.index[n.stores.carrier == f"DR {dev}"]
        if DR_i.empty:
            continue

        up_i = n.links.index[n.links.carrier == f"DR {dev} up"]
        down_i = n.links.index[n.links.carrier == f"DR {dev} down"]

        L_up = n.model["Link-p"].sel(name=up_i)
        L_up = L_up.assign_coords(name=n.links.loc[up_i, "bus0"].values)
        L_down = n.model["Link-p"].sel(name=down_i)
        L_down = L_down.assign_coords(name=n.links.loc[down_i, "bus1"].values)

        DR_act = n.model["Store-e"].sel(name=DR_i)
        DR_act = DR_act.assign_coords(name=n.stores.loc[DR_i, "bus"].map(n.buses.location).values)

        sum_ups = sum(L_up.roll(snapshot = -k) for k in range (1, S+1))
        n.model.add_constraints(sum_ups >= -DR_act, name=f"DR-{dev}-delayed")

        sum_downs = sum(L_down.roll(snapshot = -k) for k in range (1, S+1))
        n.model.add_constraints(sum_downs >= DR_act, name=f"DR-{dev}-advanced")


