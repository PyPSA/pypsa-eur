# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


localrules:
    all,
    cluster_networks,


rule process_costs:
    """Collects the processed technology cost tables for all runs and planning horizons."""
    input:
        expand(
            resources("costs_{horizon}_processed.csv"),
            run=config["run"]["name"],
            horizon=config["planning_horizons"],
        ),


rule cluster_networks:
    """Collects the clustered networks and bus maps for all runs."""
    input:
        expand(
            resources("networks/clustered.nc"),
            run=config["run"]["name"],
        ),
        expand(
            resources("busmap.csv"),
            run=config["run"]["name"],
        ),


rule compose_networks:
    """Collects the composed networks for all runs and planning horizons."""
    input:
        expand(
            resources("networks/composed_{horizon}.nc"),
            run=config["run"]["name"],
            horizon=config["planning_horizons"],
        ),


rule solve_networks:
    """Collects the solved networks for all runs at the final planning horizon."""
    input:
        expand(
            RESULTS + "networks/solved_{horizon}.nc",
            run=config["run"]["name"],
            horizon=config["planning_horizons"][-1],
        ),


rule solve_operations_networks:
    """Collects the operational dispatch networks for all runs and planning horizons."""
    input:
        expand(
            RESULTS + "networks/operations_{horizon}.nc",
            run=config["run"]["name"],
            horizon=(
                config["planning_horizons"][-1]
                if config["foresight"] == "perfect"
                else config["planning_horizons"]
            ),
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
        RESULTS + f"maps/{kind}/balance_map_{{carrier}}_{{horizon}}.{ext}",
        run=config["run"]["name"],
        horizon=config["planning_horizons"],
        carrier=config_provider("plotting", cfg_key, "bus_carriers")(w),
    )


rule plot_balance_maps:
    """Collects the static and interactive balance maps for all runs, horizons and carriers."""
    input:
        static=lambda w: balance_map_paths("static", w),
        interactive=lambda w: balance_map_paths("interactive", w),


rule plot_balance_maps_static:
    """Collects the static balance maps for all runs, horizons and carriers."""
    input:
        lambda w: balance_map_paths("static", w),


rule plot_balance_maps_interactive:
    """Collects the interactive balance maps for all runs, horizons and carriers."""
    input:
        lambda w: balance_map_paths("interactive", w),


rule plot_power_networks:
    """Collects the clustered power network maps for all runs."""
    input:
        (
            expand(
                resources("maps/clustered_network.pdf"),
                run=config["run"]["name"],
            )
            if config["foresight"] != "perfect"
            else []
        ),
