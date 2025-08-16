# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

def ptes_operation_profiles(w):
    """
    Return a dict of only the PTES profiles that are enabled in config,
    keyed by the same names you’d have used in `input:`
    """
    profiles = {}
    if config_provider(
        "sector", "district_heating", "ptes", "discharge_boosting_required"
    )(w):
        profiles["boost_per_discharge_profile"] = resources(
            "boost_per_discharge_profile_base_s_{clusters}_{planning_horizons}.nc"
        )
        profiles["cop_profiles"] = resources(
            "cop_profiles_base_s_{clusters}_{planning_horizons}.nc"
        )

    if config_provider(
        "sector", "district_heating", "ptes", "charge_boosting_required"
    )(w):
        profiles["boost_per_charge_profile"] = resources(
            "boost_per_charge_profile_base_s_{clusters}_{planning_horizons}.nc"
        )

    return profiles


rule solve_sector_network:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        unpack(ptes_operation_profiles),
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
    output:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS
        + "configs/config.base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.yaml",
    shadow:
        shadow_config
    log:
        solver=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"
