# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from numpy import atleast_1d


localrules:
    all,
    cluster_networks,


rule process_costs:
    input:
        lambda w: (
            expand(
                resources(
                    f"costs_{config_provider('costs', 'year')(w)}_processed.csv"
                ),
                run=config["run"]["name"],
            )
            if config_provider("foresight")(w) == "overnight"
            else expand(
                resources("costs_{horizon}_processed.csv"),
                **config["scenario"],
                run=config["run"]["name"],
            )
        ),


rule cluster_networks:
    input:
        expand(
            resources("networks/clustered.nc"),
            run=config["run"]["name"],
        ),


rule compose_networks:
    input:
        expand(
            resources("networks/composed_{horizon}.nc"),
            run=config["run"]["name"],
            horizon=config["planning_horizons"],
        ),


rule solve_networks:
    input:
        expand(
            RESULTS + "networks/solved_{horizon}.nc",
            run=config["run"]["name"],
            horizon=atleast_1d(config["planning_horizons"])[-1],
        ),


rule plot_balance_maps:
    input:
        lambda w: (
            expand(
                (RESULTS + "maps/{carrier}_balance_map_{horizon}.pdf"),
                run=config["run"]["name"],
                horizon=config["planning_horizons"],
                carrier=config_provider("plotting", "balance_map", "bus_carriers")(w),
            )
            if config["foresight"] != "perfect"
            else []
        ),


rule plot_power_networks:
    input:
        (
            expand(
                resources("maps/clustered_network.pdf"),
                run=config["run"]["name"],
            )
            if config["foresight"] != "perfect"
            else []
        ),
