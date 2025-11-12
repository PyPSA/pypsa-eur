# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from numpy import atleast_1d


if config["foresight"] != "perfect":

    rule plot_base_network:
        params:
            plotting=config_provider("plotting"),
        input:
            network=resources("networks/base.nc"),
            regions_onshore=resources("regions_onshore.geojson"),
        output:
            map=resources("maps/base_network.pdf"),
        threads: 1
        resources:
            mem_mb=4000,
        benchmark:
            benchmarks("plot_base_network")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_base_network.py"

    rule plot_clustered_network:
        params:
            plotting=config_provider("plotting"),
            transmission_limit=config_provider("electricity", "transmission_limit"),
        input:
            network=resources("networks/clustered.nc"),
            regions=resources("regions_onshore.geojson"),
        output:
            map=resources("maps/clustered_network.pdf"),
        threads: 1
        resources:
            mem_mb=4000,
        benchmark:
            benchmarks("plot_power_network")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_power_network.py"

    rule plot_power_network:
        params:
            plotting=config_provider("plotting"),
            transmission_limit=config_provider("electricity", "transmission_limit"),
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            regions=resources("regions_onshore.geojson"),
        output:
            map=RESULTS + "maps/power_network_{horizon}.pdf",
        threads: 2
        resources:
            mem_mb=10000,
        log:
            RESULTS + "logs/plot_power_network_{horizon}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_power_network_{horizon}")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_power_network.py"

    rule plot_hydrogen_network:
        params:
            plotting=config_provider("plotting"),
            foresight=config_provider("foresight"),
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            regions=resources("regions_onshore.geojson"),
        output:
            map=RESULTS + "maps/h2_network_{horizon}.pdf",
        threads: 2
        resources:
            mem_mb=10000,
        log:
            RESULTS + "logs/plot_hydrogen_network_{horizon}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_hydrogen_network_{horizon}")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_hydrogen_network.py"

    rule plot_gas_network:
        params:
            plotting=config_provider("plotting"),
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            regions=resources("regions_onshore.geojson"),
        output:
            map=RESULTS + "maps/ch4_network_{horizon}.pdf",
        threads: 2
        resources:
            mem_mb=10000,
        log:
            RESULTS + "logs/plot_gas_network_{horizon}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_gas_network_{horizon}")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_gas_network.py"

    rule plot_balance_map:
        params:
            plotting=config_provider("plotting"),
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            regions=resources("regions_onshore.geojson"),
        output:
            RESULTS + "maps/{carrier}_balance_map_{horizon}.pdf",
        threads: 1
        resources:
            mem_mb=8000,
        log:
            RESULTS + "logs/plot_balance_map_{horizon}-{carrier}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_balance_map_{horizon}-{carrier}")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_balance_map.py"

    rule plot_heat_source_map:
        params:
            plotting=config_provider("plotting"),
            heat_sources=config_provider("sector", "heat_pump_sources"),
        input:
            regions=resources("regions_onshore.geojson"),
            heat_source_temperature=lambda w: (
                resources("temp_" + w.carrier + "_temporal_aggregate.nc")
                if w.carrier in ["river_water", "sea_water", "ambient_air"]
                else []
            ),
            heat_source_energy=lambda w: (
                resources("heat_source_energy_" + w.carrier + "_temporal_aggregate.nc")
                if w.carrier in ["river_water"]
                else []
            ),
        output:
            temp_map=RESULTS
            + "maps/{horizon}-heat_source_temperature_map_{carrier}.html",
            energy_map=RESULTS + "maps/{horizon}-heat_source_energy_map_{carrier}.html",
        threads: 1
        resources:
            mem_mb=150000,
        log:
            RESULTS + "logs/plot_heat_source_map/{horizon}_{carrier}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_heat_source_map/{horizon}_{carrier}")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_heat_source_map.py"


if config["foresight"] == "perfect" and config.get("enable", {}).get(
    "perfect_foresight_postprocess", False
):

    def output_map_year(w):
        return {
            f"map_{year}": RESULTS + "maps/costs-all_" + f"{year}.pdf"
            for year in config_provider("planning_horizons")(w)
        }

    rule plot_power_network_perfect:
        params:
            plotting=config_provider("plotting"),
        input:
            network=RESULTS + "networks_brownfield_all_years.nc",
            regions=resources("regions_onshore.geojson"),
        output:
            unpack(output_map_year),
        threads: 2
        resources:
            mem_mb=10000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/plot_power_network_perfect.py"


rule make_summary:
    params:
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("planning_horizons"),
    input:
        networks=lambda w: (
            [
                RESULTS
                + f"networks/solved_{atleast_1d(config['planning_horizons'])[-1]}.nc"
            ]
            if config["foresight"] == "perfect"
            else [
                RESULTS + f"networks/solved_{h}.nc"
                for h in atleast_1d(config["planning_horizons"])
            ]
        ),
    output:
        nodal_costs=RESULTS + "csvs/nodal_costs.csv",
        nodal_capacities=RESULTS + "csvs/nodal_capacities.csv",
        nodal_capacity_factors=RESULTS + "csvs/nodal_capacity_factors.csv",
        capacity_factors=RESULTS + "csvs/capacity_factors.csv",
        costs=RESULTS + "csvs/costs.csv",
        capacities=RESULTS + "csvs/capacities.csv",
        curtailment=RESULTS + "csvs/curtailment.csv",
        energy=RESULTS + "csvs/energy.csv",
        energy_balance=RESULTS + "csvs/energy_balance.csv",
        nodal_energy_balance=RESULTS + "csvs/nodal_energy_balance.csv",
        prices=RESULTS + "csvs/prices.csv",
        weighted_prices=RESULTS + "csvs/weighted_prices.csv",
        market_values=RESULTS + "csvs/market_values.csv",
        metrics=RESULTS + "csvs/metrics.csv",
        cumulative_costs=RESULTS + "csvs/cumulative_costs.csv",
    threads: 1
    resources:
        mem_mb=16000,
    log:
        RESULTS + "logs/make_summary.log",
    benchmark:
        RESULTS + "benchmarks/make_summary"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/make_summary.py"


rule plot_summary:
    params:
        countries=config_provider("countries"),
        planning_horizons=config_provider("planning_horizons"),
        emissions_scope=config_provider("co2_budget", "emissions_scope"),
        plotting=config_provider("plotting"),
        foresight=config_provider("foresight"),
        co2_budget=config_provider("co2_budget"),
        sector=config_provider("sector"),
        RDIR=RDIR,
    input:
        costs=RESULTS + "csvs/costs.csv",
        energy=RESULTS + "csvs/energy.csv",
        balances=RESULTS + "csvs/energy_balance.csv",
        eurostat="data/eurostat/Balances-April2023",
        co2="data/bundle/eea/UNFCCC_v23.csv",
    output:
        costs=RESULTS + "graphs/costs.svg",
        energy=RESULTS + "graphs/energy.svg",
        balances=RESULTS + "graphs/balances-energy.svg",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        RESULTS + "logs/plot_summary.log",
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plot_summary.py"


rule plot_balance_timeseries:
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        network=RESULTS + "networks/solved_{horizon}.nc",
        rc="matplotlibrc",
    threads: 16
    resources:
        mem_mb=10000,
    log:
        RESULTS + "logs/plot_balance_timeseries_{horizon}.log",
    benchmark:
        RESULTS
        +"benchmarks/plot_balance_timeseries_{horizon}"
    conda:
        "../envs/environment.yaml"
    output:
        directory(RESULTS + "graphs/balance_timeseries_{horizon}"),
    script:
        "../scripts/plot_balance_timeseries.py"


rule plot_heatmap_timeseries:
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        network=RESULTS + "networks/solved_{horizon}.nc",
        rc="matplotlibrc",
    threads: 16
    resources:
        mem_mb=10000,
    log:
        RESULTS + "logs/plot_heatmap_timeseries_{horizon}.log",
    benchmark:
        RESULTS
        +"benchmarks/plot_heatmap_timeseries_{horizon}"
    conda:
        "../envs/environment.yaml"
    output:
        directory(RESULTS + "graphs/heatmap_timeseries_{horizon}"),
        {},
    script:
        "../scripts/plot_heatmap_timeseries.py"


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


rule plot_base_statistics:
    params:
        plotting=config_provider("plotting"),
        barplots=STATISTICS_BARPLOTS,
    input:
        network=RESULTS + "networks/solved.nc",
    output:
        **{
            f"{plot}_bar": RESULTS + f"figures/statistics_{plot}_bar.pdf"
            for plot in STATISTICS_BARPLOTS
        },
        barplots_touch=RESULTS + "figures/.statistics_plots",
    script:
        "../scripts/plot_statistics.py"


rule build_ambient_air_temperature_yearly_average:
    input:
        cutout=lambda w: input_cutout(w),
        regions_onshore=resources("regions_onshore.geojson"),
    output:
        average_ambient_air_temperature=resources(
            "temp_ambient_air_temporal_aggregate.nc"
        ),
    threads: 1
    resources:
        mem_mb=5000,
    log:
        RESULTS + "logs/build_ambient_air_temperature_yearly_average.log",
    benchmark:
        (RESULTS + "benchmarks/build_ambient_air_temperature_yearly_average")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_ambient_air_temperature_yearly_average.py"


rule plot_cop_profiles:
    input:
        cop_profiles=resources("cop_profiles_{horizon}.nc"),
    output:
        html=RESULTS + "graphs/cop_profiles_{horizon}.html",
    log:
        RESULTS + "logs/plot_cop_profiles_{horizon}.log",
    benchmark:
        RESULTS + "benchmarks/plot_cop_profiles/{horizon}"
    resources:
        mem_mb=10000,
    script:
        "../scripts/plot_cop_profiles/plot_cop_profiles.py"


rule plot_interactive_bus_balance:
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        bus_name_pattern=config_provider(
            "plotting", "interactive_bus_balance", "bus_name_pattern"
        ),
    input:
        network=RESULTS + "networks/solved_{horizon}.nc",
        rc="matplotlibrc",
    output:
        directory=directory(RESULTS + "graphics/interactive_bus_balance/{horizon}"),
    log:
        RESULTS + "logs/plot_interactive_bus_balance/{horizon}.log",
    benchmark:
        RESULTS + "benchmarks/plot_interactive_bus_balance/{horizon}"
    resources:
        mem_mb=20000,
    script:
        "../scripts/plot_interactive_bus_balance.py"
