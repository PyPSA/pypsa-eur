# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
rule add_existing_baseyear:
    message:
        "Adding existing infrastructure for base year for {wildcards.clusters} clusters, {wildcards.planning_horizons} planning horizons, {wildcards.opts} electric options and {wildcards.sector_opts} sector options"
    params:
        baseyear=config_provider("scenario", "planning_horizons", 0),
        sector=config_provider("sector"),
        existing_capacities=config_provider("existing_capacities"),
        carriers=config_provider("electricity", "renewable_carriers"),
        costs=config_provider("costs"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    input:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        powerplants=resources("powerplants_s_{clusters}.csv"),
        busmap_s=resources("busmap_base_s.csv"),
        busmap=resources("busmap_base_s_{clusters}.csv"),
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
        costs=lambda w: resources(
            f"costs_{config_provider('scenario', 'planning_horizons',0)(w)}_processed.csv"
        ),
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
        existing_heating_distribution=resources(
            "existing_heating_distribution_base_s_{clusters}_{planning_horizons}.csv"
        ),
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
        heating_efficiencies=resources("heating_efficiencies.csv"),
    output:
        resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_brownfield.nc"
        ),
    wildcard_constraints:
        planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
    threads: 1
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="24h"),
    log:
        logs(
            "add_existing_baseyear_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "add_existing_baseyear/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "../scripts/add_existing_baseyear.py"


def input_network_year(w):
    return {
        f"network_{year}": resources("networks/base_s_{clusters}_{opts}_{sector_opts}")
        + f"_{year}.nc"
        for year in config_provider("scenario", "planning_horizons")(w)[1:]
    }


rule prepare_perfect_foresight:
    message:
        "Preparing data for perfect foresight optimization for {wildcards.clusters} clusters, {wildcards.opts} electric options and {wildcards.sector_opts} sector options"
    params:
        costs=config_provider("costs"),
        time_resolution=config_provider("clustering", "temporal", "sector"),
    input:
        unpack(input_network_year),
        brownfield_network=lambda w: (
            resources("networks/base_s_{clusters}_{opts}_{sector_opts}_")
            + "{}_brownfield.nc".format(
                str(config_provider("scenario", "planning_horizons", 0)(w))
            )
        ),
    output:
        resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.nc"
        ),
    threads: 2
    resources:
        mem_mb=10000,
    log:
        logs("prepare_perfect_foresight_{clusters}_{opts}_{sector_opts}.log"),
    benchmark:
        benchmarks("prepare_perfect_foresight_{clusters}_{opts}_{sector_opts}")
    script:
        "../scripts/prepare_perfect_foresight.py"


rule solve_sector_network_perfect:
    message:
        "Solving sector-coupled network with perfect foresight for {wildcards.clusters} clusters, {wildcards.opts} electric options and {wildcards.sector_opts} sector options"
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
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.nc"
        ),
        costs=resources("costs_2030_processed.csv"),
    output:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.nc",
        config=RESULTS
        + "configs/config.base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.yaml",
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem"),
    shadow:
        shadow_config
    log:
        solver=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years_solver.log",
        python=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years_python.log",
        memory=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years_memory.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years}"
        )
    script:
        "../scripts/solve_network.py"


def input_networks_make_summary_perfect(w):
    return {
        f"networks_s_{clusters}_{opts}_{sector_opts}": RESULTS
        + f"networks/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.nc"
        for clusters in config_provider("scenario", "clusters")(w)
        for opts in config_provider("scenario", "opts")(w)
        for sector_opts in config_provider("scenario", "sector_opts")(w)
    }


rule make_summary_perfect:
    message:
        "Creating summary for perfect foresight optimization results"
    input:
        unpack(input_networks_make_summary_perfect),
        costs=resources("costs_2020_processed.csv"),
    output:
        nodal_capacities=RESULTS + "csvs/nodal_capacities.csv",
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
    script:
        "../scripts/make_summary_perfect.py"
