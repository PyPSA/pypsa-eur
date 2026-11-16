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


def balance_map_paths(w, kind=None):
    """Balance map targets; `kind` is "static", "interactive" or None for both."""
    if config["foresight"] == "perfect":
        return []
    kinds = ["static", "interactive"] if kind is None else [kind]
    paths = []
    for k in kinds:
        cfg_key = "balance_map" if k == "static" else "balance_map_interactive"
        ext = "pdf" if k == "static" else "html"
        paths.extend(
            expand(
                RESULTS + f"maps/{k}/balance_map_{{carrier}}_{{horizon}}.{ext}",
                run=config["run"]["name"],
                horizon=config["planning_horizons"],
                carrier=config_provider("plotting", cfg_key, "bus_carriers")(w),
            )
        )
    return paths


def sector_network_plot_paths(w):
    """Hydrogen and methane network map targets if the sector model builds them."""
    if config["foresight"] == "perfect" or not config_provider("sector", "enabled")(w):
        return []
    networks = {"H2_network": "h2", "gas_network": "ch4"}
    return [
        path
        for key, name in networks.items()
        if config_provider("sector", key)(w)
        for path in expand(
            RESULTS + f"maps/static/{name}_network_{{horizon}}.pdf",
            horizon=config["planning_horizons"],
            run=config["run"]["name"],
        )
    ]


rule plot_balance_maps:
    """Collects the static and interactive balance maps for all runs, horizons and carriers."""
    input:
        static=lambda w: balance_map_paths(w, "static"),
        interactive=lambda w: balance_map_paths(w, "interactive"),


rule plot_balance_maps_static:
    """Collects the static balance maps for all runs, horizons and carriers."""
    input:
        lambda w: balance_map_paths(w, "static"),


rule plot_balance_maps_interactive:
    """Collects the interactive balance maps for all runs, horizons and carriers."""
    input:
        lambda w: balance_map_paths(w, "interactive"),


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
