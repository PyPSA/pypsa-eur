# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule add_existing_baseyear:
    params:
        baseyear=config["scenario"]["planning_horizons"][0],
        sector=config["sector"],
        existing_capacities=config["existing_capacities"],
        costs=config["costs"],
    input:
        network=RESULTS
        + "prenetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        powerplants=RESOURCES + "powerplants.csv",
        busmap_s=RESOURCES + "busmap_elec{weather_year}_s{simpl}.csv",
        busmap=RESOURCES + "busmap_elec{weather_year}_s{simpl}_{clusters}.csv",
        clustered_pop_layout=RESOURCES
        + "pop_layout_elec{weather_year}_s{simpl}_{clusters}.csv",
        costs="data/costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
        cop_soil_total=RESOURCES
        + "cop_soil_total_elec{weather_year}_s{simpl}_{clusters}.nc",
        cop_air_total=RESOURCES
        + "cop_air_total_elec{weather_year}_s{simpl}_{clusters}.nc",
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
        existing_solar="data/existing_infrastructure/solar_capacity_IRENA.csv",
        existing_onwind="data/existing_infrastructure/onwind_capacity_IRENA.csv",
        existing_offwind="data/existing_infrastructure/offwind_capacity_IRENA.csv",
    output:
        RESULTS
        + "prenetworks-brownfield/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    wildcard_constraints:
        planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
    threads: 1
    resources:
        mem_mb=2000,
    log:
        LOGS
        + "add_existing_baseyear_elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            BENCHMARKS
            + "add_existing_baseyear/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
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
        + "prenetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        network_p=solved_previous_horizon,  #solved network at previous time step
        costs="data/costs_{planning_horizons}.csv",
        cop_soil_total=RESOURCES
        + "cop_soil_total_elec{weather_year}_s{simpl}_{clusters}.nc",
        cop_air_total=RESOURCES
        + "cop_air_total_elec{weather_year}_s{simpl}_{clusters}.nc",
    output:
        RESULTS
        + "prenetworks-brownfield/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    threads: 4
    resources:
        mem_mb=10000,
    log:
        LOGS
        + "add_brownfield_elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            BENCHMARKS
            + "add_brownfield/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
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
    input:
        network=RESULTS
        + "prenetworks-brownfield/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        costs="data/costs_{planning_horizons}.csv",
        config=RESULTS + "config.yaml",
    output:
        RESULTS
        + "postnetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    shadow:
        "shallow"
    log:
        solver=LOGS
        + "elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        python=LOGS
        + "elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: 4
    resources:
        mem_mb=config["solving"]["mem"],
        walltime=config["solving"].get("walltime", "12:00:00"),
    benchmark:
        (
            BENCHMARKS
            + "solve_sector_network/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"
