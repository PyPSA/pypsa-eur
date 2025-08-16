# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule solve_network:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        sector=config_provider("sector"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=resources("networks/composed_{horizon}.nc"),
    output:
        network=RESULTS + "networks/solved_{horizon}.nc",
    log:
        solver=normpath(RESULTS + "logs/solve_network/solver_{horizon}.log"),
        memory=RESULTS + "logs/solve_network/memory_{horizon}.log",
        python=RESULTS + "logs/solve_network/python_{horizon}.log",
    benchmark:
        (RESULTS + "benchmarks/solve_network_{horizon}.log")
    threads: solver_threads
    resources:
        mem_mb=memory,
        runtime=config_provider("solving", "runtime", default="6h"),
    shadow:
        shadow_config
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"
