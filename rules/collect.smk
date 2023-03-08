# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

localrules: all, cluster_all_networks, extra_components_all_networks, prepare_all_networks, prepare_sector_networks, solve_all_elec_networks, solve_all_networks, plot_all_networks


rule all:
    input:
        RESULTS + "graphs/costs.pdf",
    default_target: True


rule cluster_all_networks:
    input:
        expand(RESOURCES + "networks/elec_s{simpl}_{clusters}.nc", **config["scenario"]),


rule extra_components_all_networks:
    input:
        expand(
            RESOURCES + "networks/elec_s{simpl}_{clusters}_ec.nc", **config["scenario"]
        ),


rule prepare_all_networks:
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


rule solve_all_elec_networks:
    input:
        expand(
            RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule solve_all_networks:
    input:
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule plot_all_networks:
    input:
        expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"]
        ),