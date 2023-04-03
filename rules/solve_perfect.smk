# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
rule add_existing_baseyear:
    input:
        overrides="data/override_component_attrs",
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        powerplants=RESOURCES + "powerplants.csv",
        busmap_s=RESOURCES + "busmap_elec_s{simpl}.csv",
        busmap=RESOURCES + "busmap_elec_s{simpl}_{clusters}.csv",
        clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}_{clusters}.csv",
        costs="data/costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
        cop_soil_total=RESOURCES + "cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_air_total=RESOURCES + "cop_air_total_elec_s{simpl}_{clusters}.nc",
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
    input:
        overrides="data/override_component_attrs",
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


rule prepare_perfect_foresight:
    input:
        **{
            f"network_{year}": RDIR
            + "/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_"
            + f"{year}.nc"
            for year in config["scenario"]["planning_horizons"][1:]
        },
        brownfield_network=lambda w: (
            RDIR
            + "/prenetworks-brownfield/"
            + "elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_"
            + "{}.nc".format(str(config["scenario"]["planning_horizons"][0]))
        ),
        overrides="data/override_component_attrs",
    output:
        RDIR
        + "/prenetworks-brownfield/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_brownfield_all_years.nc",
    threads: 2
    resources:
        mem_mb=10000,
    script:
        "../scripts/prepare_perfect_foresight.py"
    log:
        LOGS
        + "prepare_perfect_foresight{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}.log",
    benchmark:
        (
            BENCHMARKS
            + "prepare_perfect_foresight{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/prepare_perfect_foresight.py"


rule solve_network_perfect:
    input:
        overrides="data/override_component_attrs",
        network=RDIR
        + "/prenetworks-brownfield/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_brownfield_all_years.nc",
        costs="data/costs_2030.csv",
        config=SDIR + "/configs/config.yaml",
    output:
        RDIR
        + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_brownfield_all_years.nc",
    shadow:
        "shallow"
    log:
        solver=RDIR
        + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_brownfield_all_years_solver.log",
        python=RDIR
        + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_brownfield_all_years_python.log",
        memory=RDIR
        + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_brownfield_all_years_memory.log",
    threads: 4
    resources:
        mem_mb=config["solving"]["mem"],
    benchmark:
        (
            RDIR
            + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_brownfield_all_years"
        )
    script:
        "../scripts/solve_network.py"
    conda:
        "../envs/environment.yaml"


rule make_summary_perfect:
    input:
        **{
            f"networks_{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}": RDIR
            + f"/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_brownfield_all_years.nc"
            for simpl in config["scenario"]["simpl"]
            for clusters in config["scenario"]["clusters"]
            for opts in config["scenario"]["opts"]
            for sector_opts in config["scenario"]["sector_opts"]
            for lv in config["scenario"]["lv"]
        },
        costs="data/costs_2020.csv",
        overrides="data/override_component_attrs",
    output:
        nodal_costs="results" + "/" + config["run"] + "/csvs/nodal_costs.csv",
        nodal_capacities="results" + "/" + config["run"] + "/csvs/nodal_capacities.csv",
        nodal_cfs="results" + "/" + config["run"] + "/csvs/nodal_cfs.csv",
        cfs="results" + "/" + config["run"] + "/csvs/cfs.csv",
        costs="results" + "/" + config["run"] + "/csvs/costs.csv",
        capacities="results" + "/" + config["run"] + "/csvs/capacities.csv",
        curtailment="results" + "/" + config["run"] + "/csvs/curtailment.csv",
        capital_cost="results" + "/" + config["run"] + "/csvs/capital_cost.csv",
        energy="results" + "/" + config["run"] + "/csvs/energy.csv",
        supply="results" + "/" + config["run"] + "/csvs/supply.csv",
        supply_energy="results" + "/" + config["run"] + "/csvs/supply_energy.csv",
        prices="results" + "/" + config["run"] + "/csvs/prices.csv",
        weighted_prices="results" + "/" + config["run"] + "/csvs/weighted_prices.csv",
        market_values="results" + "/" + config["run"] + "/csvs/market_values.csv",
        price_statistics="results" + "/" + config["run"] + "/csvs/price_statistics.csv",
        metrics="results" + "/" + config["run"] + "/csvs/metrics.csv",
        co2_emissions="results" + "/" + config["run"] + "/csvs/co2_emissions.csv",
    threads: 2
    resources:
        mem_mb=10000,
    script:
        "../scripts/make_summary_perfect.py"
    log:
        LOGS + "make_summary_perfect{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}.log",
    benchmark:
        (
            BENCHMARKS
            + "make_summary_perfect{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}"
        )
    conda:
        "../envs/environment.yaml"


rule plot_summary_perfect:
    input:
        **{
            f"costs_{year}": f"data/costs_{year}.csv"
            for year in config["scenario"]["planning_horizons"]
        },
        costs_csv="results" + "/" + config["run"] + "/csvs/costs.csv",
        energy="results" + "/" + config["run"] + "/csvs/energy.csv",
        balances="results" + "/" + config["run"] + "/csvs/supply_energy.csv",
        eea="data/eea/UNFCCC_v24.csv",
        countries="results" + "/" + config["run"] + "/csvs/nodal_capacities.csv",
        co2_emissions="results" + "/" + config["run"] + "/csvs/co2_emissions.csv",
        capacities="results" + "/" + config["run"] + "/csvs/capacities.csv",
        capital_cost="results" + "/" + config["run"] + "/csvs/capital_cost.csv",
    output:
        costs1="results" + "/" + config["run"] + "/graphs/costs.pdf",
        costs2="results" + "/" + config["run"] + "/graphs/costs2.pdf",
        costs3="results" + "/" + config["run"] + "/graphs/total_costs_per_year.pdf",
        # energy="results"  + '/' + config['run'] + '/graphs/energy.pdf',
        balances="results" + "/" + config["run"] + "/graphs/balances-energy.pdf",
        co2_emissions="results" + "/" + config["run"] + "/graphs/carbon_budget_plot.pdf",
        capacities="results" + "/" + config["run"] + "/graphs/capacities.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    script:
        "../scripts/plot_summary_perfect.py"
    log:
        LOGS + "plot_summary_perfect.log",
    benchmark:
        (BENCHMARKS + "plot_summary_perfect")
    conda:
        "../envs/environment.yaml"
