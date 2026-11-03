# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule solve_network:
    input:
        network=lambda w: (
            resources("networks/base_s_stoch_{clusters}_elec_{opts}.nc").format(**w)
            if config["stochastic_scenarios"]["enable"]
            else resources("networks/base_s_{clusters}_elec_{opts}.nc").format(**w)
        ),
    output:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}_{solver}.nc",
        config=RESULTS + "configs/config.base_s_{clusters}_elec_{opts}_{solver}.yaml",
        model=(
            RESULTS + "models/base_s_{clusters}_elec_{opts}_{solver}.nc"
            if config["solving"]["options"]["store_model"]
            else []
        ),
    log:
        solver=normpath(
            RESULTS + "logs/solve_network/base_s_{clusters}_elec_{opts}_{solver}_solver.log"
        ),
        memory=RESULTS
        + "logs/solve_network/base_s_{clusters}_elec_{opts}_{solver}_memory.log",
        python=RESULTS
        + "logs/solve_network/base_s_{clusters}_elec_{opts}_{solver}_python.log",
        failure_benchmark=RESULTS
        + "benchmarks/solve_network/base_s_{clusters}_elec_{opts}_{solver}.failed",
    benchmark:
        (RESULTS + "benchmarks/solve_network/base_s_{clusters}_elec_{opts}_{solver}")
    shadow:
        shadow_config
    threads: solver_threads
    resources:
        mem_mb=memory,
        runtime=config_provider("solving", "runtime", default="6h"),
        solver_slots=1,
    params:
        solving=solving_for_solver,
        smspp_benchmark=RESULTS
        + "benchmarks/solve_network/base_s_{clusters}_elec_{opts}_{solver}.smspp",
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    message:
        "Solving electricity network optimization with {wildcards.solver} for {wildcards.clusters} clusters and {wildcards.opts} electric options"
    script:
        scripts("solve_network_supervisor.py")


rule solve_operations_network:
    input:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}_{solver}.nc",
    output:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}_{solver}_op.nc",
    log:
        solver=normpath(
            RESULTS
            + "logs/solve_operations_network/base_s_{clusters}_elec_{opts}_{solver}_op_solver.log"
        ),
        python=RESULTS
        + "logs/solve_operations_network/base_s_{clusters}_elec_{opts}_{solver}_op_python.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_operations_network/base_s_{clusters}_elec_{opts}_{solver}"
        )
    shadow:
        shadow_config
    threads: 4
    resources:
        mem_mb=memory,
        runtime=config_provider("solving", "runtime", default="6h"),
    params:
        options=config_provider("solving", "options"),
        solving=solving_for_solver,
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    message:
        "Solving electricity network operations optimization with {wildcards.solver} for {wildcards.clusters} clusters and {wildcards.opts} electric options"
    script:
        scripts("solve_operations_network.py")
