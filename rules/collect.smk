# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


localrules:
    all,
    cluster_networks,
    prepare_elec_networks,
    prepare_sector_networks,
    solve_elec_networks,
    solve_sector_networks,


rule cluster_networks:
    input:
        expand(
            resources("networks/clustered.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule compose_networks:
    input:
        expand(
            resources("networks/composed.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_networks:
    input:
        expand(
            RESULTS + "networks/solved.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule plot_balance_maps:
    input:
        lambda w: expand(
            (RESULTS + "maps/balance_map_{carrier}-{planning_horizon}.pdf"),
            **config["scenario"],
            run=config["run"]["name"],
            carrier=config_provider("plotting", "balance_map", "bus_carriers")(w),
        ),


rule plot_power_networks_clustered:
    input:
        expand(
            resources("maps/power-network.pdf"),
            **config["scenario"],
            run=config["run"]["name"],
        ),
