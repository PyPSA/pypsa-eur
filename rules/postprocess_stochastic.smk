# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule plot_power_network_stochastic_scenario:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        map=RESULTS
        + "maps/static/stochastic_scenarios/{stoch_scenario}/base_s_{clusters}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
    log:
        RESULTS
        + "logs/plot_power_network/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_power_network/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    wildcard_constraints:
        stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
    threads: 2
    resources:
        mem_mb=10000,
    params:
        plotting=config_provider("plotting"),
        transmission_limit=config_provider("electricity", "transmission_limit"),
    message:
        "Plotting power network for stochastic scenario {wildcards.stoch_scenario}"
    script:
        scripts("plot_power_network.py")


rule plot_hydrogen_network_stochastic_scenario:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        map=RESULTS
        + "maps/static/stochastic_scenarios/{stoch_scenario}/base_s_{clusters}_{opts}_{sector_opts}-h2_network_{planning_horizons}.pdf",
    log:
        RESULTS
        + "logs/plot_hydrogen_network/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_hydrogen_network/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    wildcard_constraints:
        stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
    threads: 2
    resources:
        mem_mb=10000,
    params:
        plotting=config_provider("plotting"),
        foresight=config_provider("foresight"),
    message:
        "Plotting hydrogen network for stochastic scenario {wildcards.stoch_scenario}"
    script:
        scripts("plot_hydrogen_network.py")


rule plot_gas_network_stochastic_scenario:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        map=RESULTS
        + "maps/static/stochastic_scenarios/{stoch_scenario}/base_s_{clusters}_{opts}_{sector_opts}-ch4_network_{planning_horizons}.pdf",
    log:
        RESULTS
        + "logs/plot_gas_network/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_gas_network/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    wildcard_constraints:
        stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
    threads: 2
    resources:
        mem_mb=10000,
    params:
        plotting=config_provider("plotting"),
    message:
        "Plotting methane network for stochastic scenario {wildcards.stoch_scenario}"
    script:
        scripts("plot_gas_network.py")


rule plot_balance_map_stochastic_scenario:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        RESULTS
        + "maps/static/stochastic_scenarios/{stoch_scenario}/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-balance_map_{carrier}.pdf",
    log:
        RESULTS
        + "logs/plot_balance_map/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{carrier}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_balance_map/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{carrier}"
        )
    wildcard_constraints:
        stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
    threads: 1
    resources:
        mem_mb=8000,
    params:
        plotting=config_provider("plotting"),
        settings=lambda w: config_provider("plotting", "balance_map", w.carrier)(w),
    message:
        "Plotting balance map for stochastic scenario {wildcards.stoch_scenario} and carrier {wildcards.carrier}"
    script:
        scripts("plot_balance_map.py")


rule plot_balance_map_interactive_stochastic_scenario:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        RESULTS
        + "maps/interactive/stochastic_scenarios/{stoch_scenario}/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-balance_map_{carrier}.html",
    log:
        RESULTS
        + "logs/plot_balance_map_interactive/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{carrier}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_interactive_map/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{carrier}"
        )
    wildcard_constraints:
        stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
    threads: 1
    resources:
        mem_mb=8000,
    params:
        settings=lambda w: config_provider(
            "plotting", "balance_map_interactive", w.carrier
        )(w),
    message:
        "Plotting interactive balance map for stochastic scenario {wildcards.stoch_scenario} and carrier {wildcards.carrier}"
    script:
        scripts("plot_balance_map_interactive.py")


rule make_summary_stochastic_scenario:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
    output:
        nodal_costs=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/nodal_costs_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        nodal_capacities=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/nodal_capacities_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        nodal_capacity_factors=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/nodal_capacity_factors_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        capacity_factors=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/capacity_factors_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        costs=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/costs_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        capacities=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/capacities_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        curtailment=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/curtailment_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        energy=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/energy_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        energy_balance=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/energy_balance_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        nodal_energy_balance=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/nodal_energy_balance_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        prices=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/prices_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        weighted_prices=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/weighted_prices_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        market_values=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/market_values_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        metrics=RESULTS
        + "csvs/stochastic_scenarios/{stoch_scenario}/individual/metrics_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
    log:
        RESULTS
        + "logs/make_summary/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/make_summary/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    wildcard_constraints:
        stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
    threads: 1
    resources:
        mem_mb=8000,
    message:
        "Creating summary statistics for stochastic scenario {wildcards.stoch_scenario}"
    script:
        scripts("make_summary.py")


rule plot_balance_timeseries_stochastic_scenario:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/stochastic_scenarios/{stoch_scenario}/balance_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    log:
        RESULTS
        + "logs/plot_balance_timeseries/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_balance_timeseries/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    wildcard_constraints:
        stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
    threads: 16
    resources:
        mem_mb=10000,
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    message:
        "Plotting energy balance time series for stochastic scenario {wildcards.stoch_scenario}"
    script:
        scripts("plot_balance_timeseries.py")


rule plot_heatmap_timeseries_stochastic_scenario:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/stochastic_scenarios/{stoch_scenario}/heatmap_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    log:
        RESULTS
        + "logs/plot_heatmap_timeseries/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_heatmap_timeseries/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    wildcard_constraints:
        stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
    threads: 16
    resources:
        mem_mb=10000,
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    message:
        "Plotting heatmap time series for stochastic scenario {wildcards.stoch_scenario}"
    script:
        scripts("plot_heatmap_timeseries.py")


rule plot_interactive_bus_balance_stochastic_scenario:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}__sc-{stoch_scenario}.nc",
        rc="matplotlibrc",
    output:
        directory=directory(
            RESULTS
            + "graphics/stochastic_scenarios/{stoch_scenario}/interactive_bus_balance/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    log:
        RESULTS
        + "logs/plot_interactive_bus_balance/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_interactive_bus_balance/stochastic_scenarios/{stoch_scenario}_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    wildcard_constraints:
        stoch_scenario=STOCHASTIC_SCENARIO_PATTERN,
    resources:
        mem_mb=20000,
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        bus_name_pattern=config_provider(
            "plotting", "interactive_bus_balance", "bus_name_pattern"
        ),
    message:
        "Plotting interactive bus balance for stochastic scenario {wildcards.stoch_scenario}"
    script:
        scripts("plot_interactive_bus_balance.py")

