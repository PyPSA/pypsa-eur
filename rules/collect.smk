# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
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
        expand(
            resources("networks/elec_s{simpl}_{clusters}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule extra_components_networks:
    input:
        expand(
            resources("networks/elec_s{simpl}_{clusters}_ec.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule prepare_elec_networks:
    input:
        expand(
            resources("networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule prepare_sector_networks:
    input:
        expand(
            RESULTS
            + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_elec_networks:
    input:
        expand(
            RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks:
    input:
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks_perfect:
    input:
        expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule validate_elec_networks:
    input:
        expand(
            RESULTS
            + "figures/.statistics_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
            **config["scenario"],
            run=config["run"]["name"],
        ),
        expand(
            RESULTS
            + "figures/.validation_{kind}_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
            **config["scenario"],
            run=config["run"]["name"],
            kind=["production", "prices", "cross_border"],
        ),


# rule plot_resources:
#     input:
#         resources("graphics/power-network-unclustered.pdf"),
#         resources("graphics/gas-network-unclustered.pdf"),
#         resources("graphics/wind-energy-density.pdf"),
#         resources("graphics/weather-map-irradiation.pdf"),
#         resources("graphics/industrial-sites.pdf"),
#         resources("graphics/powerplants.pdf"),
#         resources("graphics/salt-caverns.pdf"),
#         expand(
#             resources("graphics/power-network-{clusters}.pdf"), **config["scenario"], run=config["run"]["name"],
#         ),
#         expand(
#             resources("graphics/salt-caverns-{clusters}-nearshore.pdf"),
#             **config["scenario"], run=config["run"]["name"],
#         ),
#         expand(
#             resources("graphics/biomass-potentials-{clusters}-biogas.pdf"),
#             **config["scenario"], run=config["run"]["name"],
#         ),
# rule plot_statistics:
#     input:
#         [
#             expand(
#                 RESULTS
#                 + "statistics/figures/comparison/country_{country}/.statistics_{carrier}_plots",
#                 country=config["plotting"].get("countries", "all"),
#                 carrier=config["plotting"].get("carriers", ["all"]),
#                 run=config["run"]["name"],
#             ),
#             expand(
#                 RESULTS
#                 + "statistics/figures/single/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}/country_{country}/.statistics_{carrier}_plots",
#                 **config["scenario"],
#                 country=config["plotting"].get("countries", "all"),
#                 carrier=config["plotting"].get("carriers", ["all"]),
#                 run=config["run"]["name"],
#             ),
#         ],
