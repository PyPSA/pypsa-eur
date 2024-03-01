# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule solve_sector_network:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=RESULTS
        + "prenetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS + "config.yaml",
    output:
        RESULTS
        + "postnetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    shadow:
        "shallow"
    log:
        solver=RESULTS
        + "logs/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem"),
        walltime=config_provider("solving", "walltime", default="12:00:00"),
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"


rule solve_operations_network_other_year:
    input:
        pre=RDIR
        + "/prenetworks/elec{weather_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
        post=RDIR
        + "/postnetworks/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
    output:
        RDIR
        + "/operations/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_{weather_year}.nc",
    shadow:
        "shallow"
    log:
        solver=RDIR
        + "/logs/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_{weather_year}_solver.log",
        python=RDIR
        + "/logs/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_{weather_year}_python.log",
    threads: 4
    resources:
        mem_mb=10000,
    benchmark:
        (
            RDIR
            + "/benchmarks/solve_operations_network/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_{weather_year}"
        )
    script:
        "../scripts/solve_operations_network_other_year.py"


# rule solve_operations_network_other_year_myopic:
#     input:
#         overrides="data/override_component_attrs",
#         pre=RDIR
#         + "/prenetworks/elec{weather_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
#         post=RDIR
#         + "/postnetworks/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
#         previous=solved_previous_year,
#     output:
#         RDIR
#         + "/operations/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_{weather_year}_myopic.nc",
#     shadow:
#         "shallow"
#     log:
#         solver=RDIR
#         + "/logs/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_{weather_year}_solver.log",
#         python=RDIR
#         + "/logs/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_{weather_year}_python.log",
#     threads: 4
#     resources:
#         mem_mb=10000,
#     benchmark:
#         (
#             RDIR
#             + "/benchmarks/solve_operations_network_myopic/elec{capacity_year}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_{weather_year}"
#         )
#     script:
#         "../scripts/solve_operations_network_other_year_myopic.py"
