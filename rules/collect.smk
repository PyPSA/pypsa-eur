# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


localrules:
    all,
    cluster_networks,
    extra_components_networks,
    prepare_elec_networks,
    prepare_sector_networks,
    solve_elec_networks,
    solve_sector_networks,


rule cluster_networks:
    input:
        expand(RESOURCES + "networks/elec_s{simpl}_{clusters}.nc", **config["scenario"]),


rule extra_components_networks:
    input:
        expand(
            RESOURCES + "networks/elec_s{simpl}_{clusters}_ec.nc", **config["scenario"]
        ),


rule prepare_elec_networks:
    input:
        expand(
            RESOURCES + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule prepare_sector_networks:
    input:
        expand(
            RESULTS
            + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule solve_elec_networks:
    input:
        expand(
            RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule solve_sector_networks:
    input:
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule solve_sector_networks_perfect:
    input:
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years.nc",
            **config["scenario"]
        ),


rule validate_elec_networks:
    input:
        expand(
            RESULTS
            + "figures/.statistics_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
            **config["scenario"]
        ),
        expand(
            RESULTS
            + "figures/.validation_{kind}_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
            **config["scenario"],
            kind=["production", "prices", "cross_border"]
        ),
