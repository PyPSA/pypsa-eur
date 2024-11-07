# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
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
            resources("networks/base_s_{clusters}_elec_l{ll}_{opts}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule prepare_sector_networks:
    input:
        expand(
            RESULTS
            + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_elec_networks:
    input:
        expand(
            RESULTS + "networks/base_s_{clusters}_elec_l{ll}_{opts}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks:
    input:
        expand(
            RESULTS
            + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks_perfect:
    input:
        expand(
            RESULTS
            + "maps/base_s_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule validate_elec_networks:
    input:
        expand(
            RESULTS + "figures/.statistics_plots_base_s_{clusters}_elec_l{ll}_{opts}",
            **config["scenario"],
            run=config["run"]["name"],
        ),
        expand(
            RESULTS
            + "figures/.validation_{kind}_plots_base_s_{clusters}_elec_l{ll}_{opts}",
            **config["scenario"],
            run=config["run"]["name"],
            kind=["production", "prices", "cross_border"],
        ),


rule plot_statistics:
    input:
        [
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
                RESULTS
                + "statistics/figures/single/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}/country_{country}/.statistics_{carrier}_plots",
                **config["scenario"],
                country=config_provider("plotting", "statistics")(run).get(
                    "countries", "all"
                ),
                carrier=config_provider("plotting", "statistics")(run).get(
                    "carriers", "all"
                ),
                run=config["run"]["name"],
            ),
            (
                expand(
                    "results/statistics/"
                    + config_provider("plotting", "statistics")(run).get(
                        "comparison_folder", "''"
                    )
                    + "/"
                    + "figures/country_{country}/.statistics_{carrier}_plots",
                    **config["scenario"],
                    country=config_provider("plotting", "statistics")(run).get(
                        "countries", "all"
                    ),
                    carrier=config_provider("plotting", "statistics")(run).get(
                        "carriers", "all"
                    ),
                    run=config["run"]["name"],
                )
                if config_provider("plotting", "statistics")(run).get(
                    "comparison_folder", ""
                )
                != ""
                else []
            ),
        ],
