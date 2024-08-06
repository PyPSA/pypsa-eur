# SPDX-FileCopyrightText: : 2023-4 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule add_existing_baseyear:
    params:
        baseyear=config_provider("scenario", "planning_horizons", 0),
        sector=config_provider("sector"),
        existing_capacities=config_provider("existing_capacities"),
        costs=config_provider("costs"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
    input:
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        powerplants=resources("powerplants.csv"),
        busmap_s=resources("busmap_elec_s{simpl}.csv"),
        busmap=resources("busmap_elec_s{simpl}_{clusters}.csv"),
        clustered_pop_layout=resources("pop_layout_elec_s{simpl}_{clusters}.csv"),
        costs=lambda w: resources(
            "costs_{}.csv".format(
                config_provider("scenario", "planning_horizons", 0)(w)
            )
        ),
        cop_profiles=resources("cop_profiles_elec_s{simpl}_{clusters}.nc"),
        existing_heating_distribution=resources(
            "existing_heating_distribution_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
        ),
    output:
        RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    wildcard_constraints:
        # TODO: The first planning_horizon needs to be aligned across scenarios
        # snakemake does not support passing functions to wildcard_constraints
        # reference: https://github.com/snakemake/snakemake/issues/2703
        planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
    threads: 1
    resources:
        mem_mb=2000,
    log:
        RESULTS
        + "logs/add_existing_baseyear_elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/add_existing_baseyear/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_existing_baseyear.py"


def input_profile_tech_brownfield(w):
    return {
        f"profile_{tech}": resources(f"profile_{tech}.nc")
        for tech in config_provider("electricity", "renewable_carriers")(w)
        if tech != "hydro"
    }


rule add_brownfield:
    params:
        H2_retrofit=config_provider("sector", "H2_retrofit"),
        H2_retrofit_capacity_per_CH4=config_provider(
            "sector", "H2_retrofit_capacity_per_CH4"
        ),
        threshold_capacity=config_provider("existing_capacities", "threshold_capacity"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        carriers=config_provider("electricity", "renewable_carriers"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
    input:
        unpack(input_profile_tech_brownfield),
        simplify_busmap=resources("busmap_elec_s{simpl}.csv"),
        cluster_busmap=resources("busmap_elec_s{simpl}_{clusters}.csv"),
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        network_p=solved_previous_horizon,  #solved network at previous time step
        costs=resources("costs_{planning_horizons}.csv"),
        cop_profiles=resources("cop_profiles_elec_s{simpl}_{clusters}.nc"),
    output:
        RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    threads: 4
    resources:
        mem_mb=10000,
    log:
        RESULTS
        + "logs/add_brownfield_elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/add_brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_brownfield.py"


ruleorder: add_existing_baseyear > add_brownfield


rule solve_sector_network_myopic:
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
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        costs=resources("costs_{planning_horizons}.csv"),
    output:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS
        + "configs/config.elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.yaml",
    shadow:
        "shallow"
    log:
        solver=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"
