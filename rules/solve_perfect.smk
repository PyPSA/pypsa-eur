# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
rule add_existing_baseyear:
    params:
        baseyear=config_provider("scenario", "planning_horizons", 0),
        sector=config_provider("sector"),
        existing_capacities=config_provider("existing_capacities"),
        costs=config_provider("costs"),
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
        cop_soil_total=resources("cop_soil_total_elec_s{simpl}_{clusters}.nc"),
        cop_air_total=resources("cop_air_total_elec_s{simpl}_{clusters}.nc"),
        existing_heating_distribution=resources(
            "existing_heating_distribution_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
        ),
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
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
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="24h"),
    log:
        logs(
            "add_existing_baseyear_elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "add_existing_baseyear/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_existing_baseyear.py"


def input_network_year(w):
    return {
        f"network_{year}": RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}"
        + f"_{year}.nc"
        for year in config_provider("scenario", "planning_horizons")(w)[1:]
    }


rule prepare_perfect_foresight:
    params:
        costs=config_provider("costs"),
        time_resolution=config_provider("clustering", "temporal", "sector"),
    input:
        unpack(input_network_year),
        brownfield_network=lambda w: (
            RESULTS
            + "prenetworks-brownfield/"
            + "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_"
            + "{}.nc".format(
                str(config_provider("scenario", "planning_horizons", 0)(w))
            )
        ),
    output:
        RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years.nc",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        logs(
            "prepare_perfect_foresight{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}.log"
        ),
    benchmark:
        benchmarks(
            "prepare_perfect_foresight{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/prepare_perfect_foresight.py"


rule solve_sector_network_perfect:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        sector=config_provider("sector"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years.nc",
        costs=resources("costs_2030.csv"),
    output:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years.nc",
        config=RESULTS
        + "configs/config.elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years.yaml",
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem"),
    shadow:
        "shallow"
    log:
        solver=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years_solver.log",
        python=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years_python.log",
        memory=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years_memory.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"


def input_networks_make_summary_perfect(w):
    return {
        f"networks_{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}": RESULTS
        + f"postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years.nc"
        for simpl in config_provider("scenario", "simpl")(w)
        for clusters in config_provider("scenario", "clusters")(w)
        for opts in config_provider("scenario", "opts")(w)
        for sector_opts in config_provider("scenario", "sector_opts")(w)
        for ll in config_provider("scenario", "ll")(w)
    }


rule make_summary_perfect:
    input:
        unpack(input_networks_make_summary_perfect),
        costs=resources("costs_2020.csv"),
    output:
        nodal_costs=RESULTS + "csvs/nodal_costs.csv",
        nodal_capacities=RESULTS + "csvs/nodal_capacities.csv",
        nodal_cfs=RESULTS + "csvs/nodal_cfs.csv",
        cfs=RESULTS + "csvs/cfs.csv",
        costs=RESULTS + "csvs/costs.csv",
        capacities=RESULTS + "csvs/capacities.csv",
        curtailment=RESULTS + "csvs/curtailment.csv",
        energy=RESULTS + "csvs/energy.csv",
        supply=RESULTS + "csvs/supply.csv",
        supply_energy=RESULTS + "csvs/supply_energy.csv",
        prices=RESULTS + "csvs/prices.csv",
        weighted_prices=RESULTS + "csvs/weighted_prices.csv",
        market_values=RESULTS + "csvs/market_values.csv",
        price_statistics=RESULTS + "csvs/price_statistics.csv",
        metrics=RESULTS + "csvs/metrics.csv",
        co2_emissions=RESULTS + "csvs/co2_emissions.csv",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        logs("make_summary_perfect.log"),
    benchmark:
        benchmarks("make_summary_perfect")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/make_summary_perfect.py"
