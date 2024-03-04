# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule solve_network:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=resources("networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"),
        config=RESULTS + "config.yaml",
    output:
        network=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        solver=normpath(
            RESULTS
            + "logs/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
        ),
        python=RESULTS
        + "logs/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
    benchmark:
        RESULTS + "benchmarks/solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    threads: solver_threads
    resources:
        mem_mb=memory,
        walltime=config_provider("solving", "walltime", default="12:00:00"),
    shadow:
        "minimal"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"


rule solve_operations_network:
    params:
        options=config_provider("solving", "options"),
    input:
        network=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        network=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op.nc",
    log:
        solver=normpath(
            RESULTS
            + "logs/solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_solver.log"
        ),
        python=RESULTS
        + "logs/solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_python.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    threads: 4
    resources:
        mem_mb=(lambda w: 10000 + 372 * int(w.clusters)),
        walltime=config_provider("solving", "walltime", default="12:00:00"),
    shadow:
        "minimal"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_operations_network.py"
