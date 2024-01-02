# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule solve_network:
    params:
        solving=config["solving"],
        foresight=config["foresight"],
        planning_horizons=config["scenario"]["planning_horizons"],
        co2_sequestration_potential=config["sector"].get(
            "co2_sequestration_potential", 200
        ),
    input:
        network=RESOURCES + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        config=RESULTS + "config.yaml",
    output:
        network=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        solver=normpath(
            LOGS + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
        ),
        python=LOGS
        + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
    benchmark:
        BENCHMARKS + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    threads: 4
    resources:
        mem_mb=memory,
        walltime=config["solving"].get("walltime", "12:00:00"),
    shadow:
        "minimal"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"


rule solve_operations_network:
    params:
        options=config["solving"]["options"],
    input:
        network=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        network=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op.nc",
    log:
        solver=normpath(
            LOGS
            + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_solver.log"
        ),
        python=LOGS
        + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_python.log",
    benchmark:
        (
            BENCHMARKS
            + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    threads: 4
    resources:
        mem_mb=(lambda w: 10000 + 372 * int(w.clusters)),
        walltime=config["solving"].get("walltime", "12:00:00"),
    shadow:
        "minimal"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_operations_network.py"
