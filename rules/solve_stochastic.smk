# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from pathlib import Path
import re
import yaml


rule build_stochastic_network:
    input:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
    output:
        network=resources(
            "networks/base_s_stoch_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        config=RESULTS
        + "configs/config.base_s_stoch_{clusters}_{opts}_{sector_opts}_{planning_horizons}.yaml",
    log:
        python=RESULTS
        + "logs/build_stochastic_network/base_s_stoch_{clusters}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/build_stochastic_network/base_s_stoch_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    shadow:
        shadow_config
    threads: 1
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="1h"),
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        stochastic_scenarios=config_provider("stochastic_scenarios"),
    message:
        "Building stochastic sector-coupled network"
    script:
        scripts("build_stochastic_network.py")


rule solve_sector_network:
    input:
        network=resources(
            "networks/base_s_stoch_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
    output:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS
        + "configs/config.base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.yaml",
        model=(
            RESULTS
            + "models/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
            if config["solving"]["options"]["store_model"]
            else []
        ),
    log:
        solver=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    shadow:
        shadow_config
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    message:
        "Solving stochastic sector-coupled network with overnight investment optimization for {wildcards.clusters} clusters, {wildcards.planning_horizons} planning horizons, {wildcards.opts} electric options and {wildcards.sector_opts} sector options"
    script:
        scripts("solve_network.py")


if config["stochastic_scenarios"]["export"]["average"]:

    rule export_stochastic_average:
        input:
            network=RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        output:
            network=RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__avg.nc",
        log:
            python=RESULTS
            + "logs/export_stochastic_views/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__avg.log",
        benchmark:
            (
                RESULTS
                + "benchmarks/export_stochastic_views/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__avg"
            )
        threads: 1
        resources:
            mem_mb=8000,
        params:
            scenarios_file=config["stochastic_scenarios"]["file"],
            mode="average",
        message:
            "Exporting deterministic average view from stochastic solution"
        script:
            scripts("export_stochastic_views.py")


if config["stochastic_scenarios"]["export"]["scenarios"]:

    rule export_stochastic_scenario:
        input:
            network=RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        output:
            network=RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
        log:
            python=RESULTS
            + "logs/export_stochastic_views/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.log",
        benchmark:
            (
                RESULTS
                + "benchmarks/export_stochastic_views/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}"
            )
        wildcard_constraints:
            stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
        threads: 1
        resources:
            mem_mb=8000,
        params:
            scenarios_file=config["stochastic_scenarios"]["file"],
            mode="scenario",
            scenario=lambda w: w.stoch_scenario,
        message:
            "Exporting deterministic scenario view from stochastic solution for {wildcards.stoch_scenario}"
        script:
            scripts("export_stochastic_views.py")
