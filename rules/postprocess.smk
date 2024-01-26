# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


localrules:
    copy_config,


if config["foresight"] != "perfect":

    rule plot_power_network_clustered:
        params:
            plotting=config["plotting"],
        input:
            network=RESOURCES + "networks/elec_s{simpl}_{clusters}.nc",
            regions_onshore=RESOURCES
            + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        output:
            map=RESULTS + "maps/power-network-s{simpl}-{clusters}.pdf",
        threads: 1
        resources:
            mem_mb=4000,
        benchmark:
            BENCHMARKS + "plot_power_network_clustered/elec_s{simpl}_{clusters}"
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_power_network_clustered.py"

    rule plot_power_network:
        params:
            plotting=config["plotting"],
        input:
            network=RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            regions=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        output:
            map=RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
        threads: 2
        resources:
            mem_mb=10000,
        benchmark:
            (
                BENCHMARKS
                + "plot_power_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
            )
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_power_network.py"

    rule plot_hydrogen_network:
        params:
            plotting=config["plotting"],
        input:
            network=RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            regions=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        output:
            map=RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-h2_network_{planning_horizons}.pdf",
        threads: 2
        resources:
            mem_mb=10000,
        benchmark:
            (
                BENCHMARKS
                + "plot_hydrogen_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
            )
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_hydrogen_network.py"

    rule plot_gas_network:
        params:
            plotting=config["plotting"],
        input:
            network=RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            regions=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        output:
            map=RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-ch4_network_{planning_horizons}.pdf",
        threads: 2
        resources:
            mem_mb=10000,
        benchmark:
            (
                BENCHMARKS
                + "plot_gas_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
            )
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_gas_network.py"


if config["foresight"] == "perfect":

    rule plot_power_network_perfect:
        params:
            plotting=config["plotting"],
        input:
            network=RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years.nc",
            regions=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        output:
            **{
                f"map_{year}": RESULTS
                + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_"
                + f"{year}.pdf"
                for year in config["scenario"]["planning_horizons"]
            },
        threads: 2
        resources:
            mem_mb=10000,
        benchmark:
            BENCHMARKS
            +"postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years_benchmark"
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_power_network_perfect.py"


rule copy_config:
    params:
        RDIR=RDIR,
    output:
        RESULTS + "config.yaml",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        BENCHMARKS + "copy_config"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/copy_config.py"


rule make_summary:
    params:
        foresight=config["foresight"],
        costs=config["costs"],
        snapshots={k: config["snapshots"][k] for k in ["start", "end", "inclusive"]},
        scenario=config["scenario"],
        RDIR=RDIR,
    input:
        expand(
            RESULTS + "maps/power-network-s{simpl}-{clusters}.pdf",
            **config["scenario"]
        ),
        networks=expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),
        costs="data/costs_{}.csv".format(config["costs"]["year"])
        if config["foresight"] == "overnight"
        else "data/costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
        ac_plot=expand(
            RESULTS + "maps/power-network-s{simpl}-{clusters}.pdf",
            **config["scenario"]
        ),
        costs_plot=expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"]
        ),
        h2_plot=expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-h2_network_{planning_horizons}.pdf"
            if config["sector"]["H2_network"]
            else [],
            **config["scenario"]
        ),
        ch4_plot=expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-ch4_network_{planning_horizons}.pdf"
            if config["sector"]["gas_network"]
            else [],
            **config["scenario"]
        ),
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
    threads: 2
    resources:
        mem_mb=10000,
    log:
        LOGS + "make_summary.log",
    benchmark:
        BENCHMARKS + "make_summary"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/make_summary.py"


rule plot_summary:
    params:
        countries=config["countries"],
        planning_horizons=config["scenario"]["planning_horizons"],
        sector_opts=config["scenario"]["sector_opts"],
        emissions_scope=config["energy"]["emissions"],
        eurostat_report_year=config["energy"]["eurostat_report_year"],
        plotting=config["plotting"],
        RDIR=RDIR,
    input:
        costs=RESULTS + "csvs/costs.csv",
        energy=RESULTS + "csvs/energy.csv",
        balances=RESULTS + "csvs/supply_energy.csv",
        eurostat=input_eurostat,
        co2="data/bundle-sector/eea/UNFCCC_v23.csv",
    output:
        costs=RESULTS + "graphs/costs.pdf",
        energy=RESULTS + "graphs/energy.pdf",
        balances=RESULTS + "graphs/balances-energy.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        LOGS + "plot_summary.log",
    benchmark:
        BENCHMARKS + "plot_summary"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plot_summary.py"


STATISTICS_BARPLOTS = [
    "capacity_factor",
    "installed_capacity",
    "optimal_capacity",
    "capital_expenditure",
    "operational_expenditure",
    "curtailment",
    "supply",
    "withdrawal",
    "market_value",
]


rule plot_elec_statistics:
    params:
        plotting=config["plotting"],
        barplots=STATISTICS_BARPLOTS,
    input:
        network=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        **{
            f"{plot}_bar": RESULTS
            + f"figures/statistics_{plot}_bar_elec_s{{simpl}}_{{clusters}}_ec_l{{ll}}_{{opts}}.pdf"
            for plot in STATISTICS_BARPLOTS
        },
        barplots_touch=RESULTS
        + "figures/.statistics_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
    script:
        "../scripts/plot_statistics.py"
