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
            resources("networks/base_s_{clusters}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule prepare_elec_networks:
    input:
        expand(
            resources("networks/base_s_{clusters}_elec_{opts}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule prepare_sector_networks:
    input:
        expand(
            resources(
                "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_elec_networks:
    input:
        expand(
            RESULTS + "networks/base_s_{clusters}_elec_{opts}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks:
    input:
        expand(
            RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks_perfect:
    input:
        expand(
            RESULTS
            + "maps/base_s_{clusters}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule validate_elec_networks:
    input:
        expand(
            RESULTS + "figures/.statistics_plots_base_s_{clusters}_elec_{opts}",
            **config["scenario"],
            run=config["run"]["name"],
        ),
        expand(
            RESULTS + "figures/.validation_{kind}_plots_base_s_{clusters}_elec_{opts}",
            **config["scenario"],
            run=config["run"]["name"],
            kind=["production", "prices", "cross_border"],
        ),


rule plot_statistics:
    input:
        lambda w: expand(
            (
                RESULTS
                + "statistics/csv/base_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}/country_{country}/.statistics_{carrier}_csv"
            ),
            **config["scenario"],
            run=config["run"]["name"],
            country=config_provider("plotting", "statistics")(w).get(
                "countries", "all"
            ),
            carrier=config_provider("plotting", "statistics")(w).get(
                "carriers", "all"
            ),
            allow_missing=True,
        ),
        expand(
            RESULTS
            + "statistics/figures/single/base_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}/country_{country}/.statistics_{carrier}_plots",
            **config["scenario"],
            country=config_provider("plotting", "statistics")(run).get(
                "countries", "all"
            ),
            carrier=config_provider("plotting", "statistics")(run).get(
                "carriers", "all"
            ),
            run=config["run"]["name"],
        ),
        expand(
            RESULTS
            + "statistics/figures/comparison/country_{country}/.statistics_{carrier}_plots",
            country=config_provider("plotting", "statistics")(run).get(
                "countries", "all"
            ),
            carrier=config_provider("plotting", "statistics")(run).get(
                "carriers", "all"
            ),
            run=config["run"]["name"],
        ),
        expand(
            "results/statistics/"
            + config_provider("plotting", "statistics")(run).get(
                "comparison_folder", "''"
            )
            + "/"
            + "figures/country_{country}/.statistics_{carrier}_plots",
            country=config_provider("plotting", "statistics")(run).get(
                "countries", "all"
            ),
            carrier=config_provider("plotting", "statistics")(run).get(
                "carriers", "all"
            ),
            run=config["run"]["name"],
        )
        if config_provider("plotting", "statistics")(run).get("comparison_folder", "")
        != ""
        else [],
