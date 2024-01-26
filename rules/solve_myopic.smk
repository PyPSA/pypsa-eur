# SPDX-FileCopyrightText: : 2023-4 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule build_existing_heating_distribution:
    params:
        baseyear=config["scenario"]["planning_horizons"][0],
        sector=config["sector"],
        existing_capacities=config["existing_capacities"],
    input:
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
        clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}_{clusters}.csv",
        clustered_pop_energy_layout=RESOURCES
        + "pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
        district_heat_share=RESOURCES
        + "district_heat_share_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    output:
        existing_heating_distribution=RESOURCES
        + "existing_heating_distribution_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    wildcard_constraints:
        planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
    threads: 1
    resources:
        mem_mb=2000,
    log:
        LOGS
        + "build_existing_heating_distribution_elec_s{simpl}_{clusters}_{planning_horizons}.log",
    benchmark:
        (
            BENCHMARKS
            + "build_existing_heating_distribution/elec_s{simpl}_{clusters}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_existing_heating_distribution.py"


rule add_existing_baseyear:
    params:
        baseyear=config["scenario"]["planning_horizons"][0],
        sector=config["sector"],
        existing_capacities=config["existing_capacities"],
        costs=config["costs"],
    input:
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        powerplants=RESOURCES + "powerplants.csv",
        busmap_s=RESOURCES + "busmap_elec_s{simpl}.csv",
        busmap=RESOURCES + "busmap_elec_s{simpl}_{clusters}.csv",
        clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}_{clusters}.csv",
        costs="data/costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
        cop_soil_total=RESOURCES + "cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_air_total=RESOURCES + "cop_air_total_elec_s{simpl}_{clusters}.nc",
        existing_heating_distribution=RESOURCES
        + "existing_heating_distribution_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        existing_solar="data/existing_infrastructure/solar_capacity_IRENA.csv",
        existing_onwind="data/existing_infrastructure/onwind_capacity_IRENA.csv",
        existing_offwind="data/existing_infrastructure/offwind_capacity_IRENA.csv",
    output:
        RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    wildcard_constraints:
        planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
    threads: 1
    resources:
        mem_mb=2000,
    log:
        LOGS
        + "add_existing_baseyear_elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            BENCHMARKS
            + "add_existing_baseyear/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_existing_baseyear.py"


rule add_brownfield:
    params:
        H2_retrofit=config["sector"]["H2_retrofit"],
        H2_retrofit_capacity_per_CH4=config["sector"]["H2_retrofit_capacity_per_CH4"],
        threshold_capacity=config["existing_capacities"]["threshold_capacity"],
    input:
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        network_p=solved_previous_horizon,  #solved network at previous time step
        costs="data/costs_{planning_horizons}.csv",
        cop_soil_total=RESOURCES + "cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_air_total=RESOURCES + "cop_air_total_elec_s{simpl}_{clusters}.nc",
    output:
        RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    threads: 4
    resources:
        mem_mb=10000,
    log:
        LOGS
        + "add_brownfield_elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            BENCHMARKS
            + "add_brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_brownfield.py"


ruleorder: add_existing_baseyear > add_brownfield


rule solve_sector_network_myopic:
    params:
        solving=config["solving"],
        foresight=config["foresight"],
        planning_horizons=config["scenario"]["planning_horizons"],
        co2_sequestration_potential=config["sector"].get(
            "co2_sequestration_potential", 200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        costs="data/costs_{planning_horizons}.csv",
        config=RESULTS + "config.yaml",
    output:
        RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    shadow:
        "shallow"
    log:
        solver=LOGS
        + "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        python=LOGS
        + "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: solver_threads
    resources:
        mem_mb=config["solving"]["mem"],
        walltime=config["solving"].get("walltime", "12:00:00"),
    benchmark:
        (
            BENCHMARKS
            + "solve_sector_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"
