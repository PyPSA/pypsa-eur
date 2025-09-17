# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


localrules:
    all,
    cluster_networks,


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
            resources("networks/composed_{horizon}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
            horizon=config["planning_horizons"],
        ),


rule solve_networks:
    input:
        expand(
            RESULTS + "networks/solved_{horizon}.nc",
            **config["scenario"],
            run=config["run"]["name"],
            horizon=config["planning_horizons"],
        ),


rule plot_balance_maps:
    input:
        lambda w: expand(
            (RESULTS + "maps/{carrier}_balance_map_{horizon}.pdf"),
            **config["scenario"],
            run=config["run"]["name"],
            horizon=config["planning_horizons"],
            carrier=config_provider("plotting", "balance_map", "bus_carriers")(w),
        ),


rule plot_power_networks:
    input:
        expand(
            resources("maps/power-network.pdf"),
            **config["scenario"],
            run=config["run"]["name"],
        ),
