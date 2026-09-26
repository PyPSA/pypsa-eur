# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

if config["foresight"] != "perfect":

    rule plot_base_network:
        """Plots the base transmission network topology with line and link capacities on a map."""
        input:
            network=resources("networks/base.nc"),
            onshore_regions=resources("onshore_regions_base.geojson"),
        output:
            map=resources("maps/base_network.pdf"),
        benchmark:
            benchmarks("plot_base_network")
        threads: 1
        resources:
            mem_mb=4000,
        params:
            plotting=config_provider("plotting"),
        script:
            scripts("plot_base_network.py")

    rule plot_clustered_network:
        """Plots the clustered transmission network with existing and planned HVDC links on a map."""
        input:
            network=resources("networks/clustered.nc"),
            regions=resources("onshore_regions.geojson"),
        output:
            map=resources("maps/clustered_network.pdf"),
        benchmark:
            benchmarks("plot_clustered_network")
        threads: 1
        resources:
            mem_mb=4000,
        params:
            plotting=config_provider("plotting"),
        script:
            scripts("plot_power_network_clustered.py")

    rule plot_power_network:
        """Plots the optimised power network with regional technology costs and grid expansion."""
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            regions=resources("onshore_regions.geojson"),
        output:
            map=RESULTS + "maps/static/power_network_{horizon}.pdf",
        log:
            RESULTS + "logs/plot_power_network_{horizon}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_power_network_{horizon}")
        threads: 2
        resources:
            mem_mb=10000,
        params:
            plotting=config_provider("plotting"),
            transmission_limit=config_provider("electricity", "transmission_limit"),
        script:
            scripts("plot_power_network.py")

    rule plot_hydrogen_network:
        """Plots the optimised hydrogen pipelines, storage, electrolysis and fuel cells on a map."""
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            regions=resources("onshore_regions.geojson"),
        output:
            map=RESULTS + "maps/static/h2_network_{horizon}.pdf",
        log:
            RESULTS + "logs/plot_hydrogen_network_{horizon}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_hydrogen_network_{horizon}")
        threads: 2
        resources:
            mem_mb=10000,
        params:
            plotting=config_provider("plotting"),
            foresight=config_provider("foresight"),
        script:
            scripts("plot_hydrogen_network.py")

    rule plot_gas_network:
        """Plots the optimised methane pipeline network and regional gas sources on a map."""
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            regions=resources("onshore_regions.geojson"),
        output:
            map=RESULTS + "maps/static/ch4_network_{horizon}.pdf",
        log:
            RESULTS + "logs/plot_gas_network_{horizon}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_gas_network_{horizon}")
        threads: 2
        resources:
            mem_mb=10000,
        params:
            plotting=config_provider("plotting"),
        script:
            scripts("plot_gas_network.py")

    rule plot_balance_map:
        """Plots a static map of the nodal energy balance, flows and prices for one bus carrier."""
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            regions=resources("onshore_regions.geojson"),
        output:
            RESULTS + "maps/static/balance_map_{carrier}_{horizon}.pdf",
        log:
            RESULTS + "logs/plot_balance_map_{horizon}_{carrier}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_balance_map_{horizon}_{carrier}")
        threads: 1
        resources:
            mem_mb=8000,
        params:
            plotting=config_provider("plotting"),
            settings=lambda w: config_provider("plotting", "balance_map", w.carrier),
        script:
            scripts("plot_balance_map.py")

    rule plot_balance_map_interactive:
        """Plots an interactive map of the nodal energy balance, flows and prices for one bus carrier."""
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            regions=resources("onshore_regions.geojson"),
        output:
            RESULTS + "maps/interactive/balance_map_{carrier}_{horizon}.html",
        log:
            RESULTS + "logs/plot_balance_map_interactive/{horizon}_{carrier}.log",
        benchmark:
            RESULTS + "benchmarks/plot_balance_map_interactive/{horizon}_{carrier}"
        threads: 1
        resources:
            mem_mb=8000,
        params:
            settings=lambda w: config_provider(
                "plotting", "balance_map_interactive", w.carrier
            ),
        script:
            scripts("plot_balance_map_interactive.py")

    rule plot_heat_source_map:
        """Plots interactive maps of heat source temperatures and energy potentials by region."""
        input:
            regions=resources("onshore_regions.geojson"),
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
            + "maps/static/heat_source_temperature_map_{carrier}_{horizon}.html",
            energy_map=RESULTS
            + "maps/static/heat_source_energy_map_{carrier}_{horizon}.html",
        log:
            RESULTS + "logs/plot_heat_source_map/{carrier}_{horizon}.log",
        benchmark:
            (RESULTS + "benchmarks/plot_heat_source_map/{carrier}_{horizon}")
        threads: 1
        resources:
            mem_mb=150000,
        script:
            scripts("plot_heat_source_map.py")


rule make_summary:
    """Computes summary tables of costs, capacities, energy balances, prices and metrics."""
    input:
        networks=lambda w: (
            [RESULTS + f"networks/solved_{config['planning_horizons'][-1]}.nc"]
            if config["foresight"] == "perfect"
            else [
                RESULTS + f"networks/solved_{h}.nc"
                for h in config["planning_horizons"]
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
    log:
        RESULTS + "logs/make_summary.log",
    benchmark:
        RESULTS + "benchmarks/make_summary"
    threads: 1
    resources:
        mem_mb=16000,
    params:
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("planning_horizons"),
    script:
        scripts("make_summary.py")


rule plot_summary:
    """Plots stacked bar charts of system costs, energy and balances across planning horizons."""
    input:
        costs=RESULTS + "csvs/costs.csv",
        energy=RESULTS + "csvs/energy.csv",
        balances=RESULTS + "csvs/energy_balance.csv",
        eurostat=resources("eurostat_energy_balances.csv"),
        co2=rules.retrieve_ghg_emissions.output["csv"],
    output:
        costs=RESULTS + "graphs/costs.pdf",
        energy=RESULTS + "graphs/energy.pdf",
        balances=RESULTS + "graphs/balances_energy.pdf",
    log:
        RESULTS + "logs/plot_summary.log",
    threads: 2
    resources:
        mem_mb=10000,
    params:
        countries=config_provider("countries"),
        planning_horizons=config_provider("planning_horizons"),
        emissions_scope=config_provider("co2_budget", "emissions_scope"),
        plotting=config_provider("plotting"),
        foresight=config_provider("foresight"),
        sector=config_provider("sector"),
        RDIR=RDIR,
    script:
        scripts("plot_summary.py")


rule plot_balance_timeseries:
    """Plots stacked energy balance time series per carrier at annual and monthly resolution."""
    input:
        network=RESULTS + "networks/solved_{horizon}.nc",
        rc="matplotlibrc",
    output:
        directory(RESULTS + "graphs/balance_timeseries_{horizon}"),
    log:
        RESULTS + "logs/plot_balance_timeseries_{horizon}.log",
    benchmark:
        RESULTS + "benchmarks/plot_balance_timeseries_{horizon}"
    threads: 16
    resources:
        mem_mb=10000,
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    script:
        scripts("plot_balance_timeseries.py")


rule plot_heatmap_timeseries:
    """Plots hour-by-day heatmaps of utilisation rates, marginal prices and storage levels."""
    input:
        network=RESULTS + "networks/solved_{horizon}.nc",
        rc="matplotlibrc",
    output:
        directory(RESULTS + "graphs/heatmap_timeseries_{horizon}"),
    log:
        RESULTS + "logs/plot_heatmap_timeseries_{horizon}.log",
    benchmark:
        RESULTS + "benchmarks/plot_heatmap_timeseries_{horizon}"
    threads: 16
    resources:
        mem_mb=10000,
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    script:
        scripts("plot_heatmap_timeseries.py")


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
    """Plots bar charts of capacities, costs, curtailment, supply and market values per carrier."""
    input:
        network=RESULTS + "networks/solved_{horizon}.nc",
    output:
        **{
            f"{plot}_bar": RESULTS + f"figures/statistics_{plot}_bar_{{horizon}}.pdf"
            for plot in STATISTICS_BARPLOTS
        },
        barplots_touch=RESULTS + "figures/.statistics_plots_{horizon}",
    params:
        plotting=config_provider("plotting"),
        barplots=STATISTICS_BARPLOTS,
    script:
        scripts("plot_statistics.py")


rule build_ambient_air_temperature_yearly_average:
    """Averages ambient air temperature over the year for each onshore grid cell."""
    input:
        cutout=lambda w: input_cutout(w),
        onshore_regions=resources("onshore_regions.geojson"),
    output:
        average_ambient_air_temperature=resources(
            "temp_ambient_air_temporal_aggregate.nc"
        ),
    log:
        logs("logs/build_ambient_air_temperature_yearly_average.log"),
    benchmark:
        benchmarks("build_ambient_air_temperature_yearly_average")
    threads: 1
    resources:
        mem_mb=5000,
    script:
        scripts("build_ambient_air_temperature_yearly_average.py")


rule plot_cop_profiles:
    """Plots interactive heat pump coefficient of performance profiles per region."""
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
        scripts("plot_cop_profiles/plot_cop_profiles.py")


rule plot_interactive_bus_balance:
    """Plots interactive energy balance time series for each bus matching a name pattern."""
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
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        bus_name_pattern=config_provider(
            "plotting", "interactive_bus_balance", "bus_name_pattern"
        ),
    script:
        scripts("plot_interactive_bus_balance.py")
