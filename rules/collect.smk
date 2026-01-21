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
                run=config["run"]["name"],
                horizon=config["planning_horizons"],
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


def balance_map_paths(kind, w):
    """
    kind = "static" or "interactive"
    """
    cfg_key = "balance_map" if kind == "static" else "balance_map_interactive"
    ext = "pdf" if kind == "static" else "html"

    if config["foresight"] == "perfect":
        return []

    return expand(
        RESULTS + f"maps/{kind}/{{carrier}}_balance_map_{{horizon}}.{ext}",
        run=config["run"]["name"],
        horizon=config["planning_horizons"],
        carrier=config_provider("plotting", cfg_key, "bus_carriers")(w),
    )


rule plot_balance_maps:
    input:
        static=lambda w: balance_map_paths("static", w),
        interactive=lambda w: balance_map_paths("interactive", w),


rule plot_balance_maps_static:
    input:
        lambda w: balance_map_paths("static", w),


rule plot_balance_maps_interactive:
    input:
        lambda w: balance_map_paths("interactive", w),


rule plot_power_networks:
    input:
        (
            expand(
                resources("maps/clustered_power_network.pdf"),
                run=config["run"]["name"],
            )
            if config["foresight"] != "perfect"
            else []
        ),
