# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

def get_network(w):
    """Return the solved network used by post-processing rules."""
    network = (
        RESULTS
        + f"networks/base_s_{w.clusters}_{w.opts}_{w.sector_opts}_{w.planning_horizons}_{w.solver}.nc"
    )

    if (
        config["stochastic_scenarios"]["enable"]
        and config["stochastic_scenarios"]["postprocess"]["use_average"]
    ):
        return network.replace(".nc", "__avg.nc")

    return network


if config["foresight"] != "perfect":

    rule plot_base_network:
        input:
            network=resources("networks/base.nc"),
            regions_onshore=resources("regions_onshore.geojson"),
        output:
            map=resources("maps/power-network.pdf"),
        benchmark:
            benchmarks("plot_base_network/base")
        threads: 1
        resources:
            mem_mb=4000,
        params:
            plotting=config_provider("plotting"),
        message:
            "Plotting base power network"
        script:
            scripts("plot_base_network.py")

    rule plot_power_network_clustered:
        input:
            network=resources("networks/base_s_{clusters}.nc"),
            regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        output:
            map=resources("maps/power-network-s-{clusters}.pdf"),
        benchmark:
            benchmarks("plot_power_network_clustered/base_s_{clusters}")
        threads: 1
        resources:
            mem_mb=4000,
        params:
            plotting=config_provider("plotting"),
        message:
            "Plotting clustered power network for {wildcards.clusters} clusters"
        script:
            scripts("plot_power_network_clustered.py")

    rule plot_power_network:
        input:
            network=get_network,
            regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        output:
            map=RESULTS
            + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-costs-all_{planning_horizons}_{solver}.pdf",
        log:
            RESULTS
            + "logs/plot_power_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.log",
        benchmark:
            (
                RESULTS
                + "benchmarks/plot_power_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
            )
        threads: 2
        resources:
            mem_mb=10000,
        params:
            plotting=config_provider("plotting"),
            transmission_limit=config_provider("electricity", "transmission_limit"),
        message:
            "Plotting power network for {wildcards.clusters} clusters, {wildcards.opts} electric options, {wildcards.sector_opts} sector options and {wildcards.planning_horizons} planning horizons"
        script:
            scripts("plot_power_network.py")

    rule plot_hydrogen_network:
        input:
            network=get_network,
            regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        output:
            map=RESULTS
            + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-h2_network_{planning_horizons}_{solver}.pdf",
        log:
            RESULTS
            + "logs/plot_hydrogen_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.log",
        benchmark:
            (
                RESULTS
                + "benchmarks/plot_hydrogen_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
            )
        threads: 2
        resources:
            mem_mb=10000,
        params:
            plotting=config_provider("plotting"),
            foresight=config_provider("foresight"),
        message:
            "Plotting hydrogen network for {wildcards.clusters} clusters, {wildcards.opts} electric options, {wildcards.sector_opts} sector options and {wildcards.planning_horizons} planning horizons"
        script:
            scripts("plot_hydrogen_network.py")

    rule plot_gas_network:
        input:
            network=get_network,
            regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        output:
            map=RESULTS
            + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-ch4_network_{planning_horizons}_{solver}.pdf",
        log:
            RESULTS
            + "logs/plot_gas_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.log",
        benchmark:
            (
                RESULTS
                + "benchmarks/plot_gas_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
            )
        threads: 2
        resources:
            mem_mb=10000,
        params:
            plotting=config_provider("plotting"),
        message:
            "Plotting methane network for {wildcards.clusters} clusters, {wildcards.opts} electric options, {wildcards.sector_opts} sector options and {wildcards.planning_horizons} planning horizon"
        script:
            scripts("plot_gas_network.py")

    rule plot_balance_map:
        input:
            network=get_network,
            regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        output:
            RESULTS
            + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}-balance_map_{carrier}.pdf",
        log:
            RESULTS
            + "logs/plot_balance_map/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}_{carrier}.log",
        benchmark:
            (
                RESULTS
                + "benchmarks/plot_balance_map/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}_{carrier}"
            )
        threads: 1
        resources:
            mem_mb=8000,
        params:
            plotting=config_provider("plotting"),
            settings=lambda w: config_provider("plotting", "balance_map", w.carrier),
        message:
            "Plotting balance map for {wildcards.clusters} clusters, {wildcards.opts} electric options, {wildcards.sector_opts} sector options, {wildcards.planning_horizons} planning horizons and {wildcards.carrier} carrier"
        script:
            scripts("plot_balance_map.py")

    rule plot_balance_map_interactive:
        input:
            network=get_network,
            regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        output:
            RESULTS
            + "maps/interactive/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}-balance_map_{carrier}.html",
        log:
            RESULTS
            + "logs/plot_balance_map_interactive/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}_{carrier}.log",
        benchmark:
            (
                RESULTS
                + "benchmarks/plot_interactive_map/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}_{carrier}"
            )
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
        input:
            regions=resources("regions_onshore_base_s_{clusters}.geojson"),
            heat_source_temperature=lambda w: (
                resources(
                    "temp_" + w.carrier + "_base_s_{clusters}_temporal_aggregate.nc"
                )
                if w.carrier in ["river_water", "sea_water", "ambient_air"]
                else []
            ),
            heat_source_energy=lambda w: (
                resources(
                    "heat_source_energy_"
                    + w.carrier
                    + "_base_s_{clusters}_temporal_aggregate.nc"
                )
                if w.carrier in ["river_water"]
                else []
            ),
        output:
            temp_map=RESULTS
            + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}-heat_source_temperature_map_{carrier}.html",
            energy_map=RESULTS
            + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}-heat_source_energy_map_{carrier}.html",
        log:
            RESULTS
            + "logs/plot_heat_source_map/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}_{carrier}.log",
        benchmark:
            (
                RESULTS
                + "benchmarks/plot_heat_source_map/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}_{carrier}"
            )
        threads: 1
        resources:
            mem_mb=150000,
        params:
            plotting=config_provider("plotting"),
            heat_sources=config_provider("sector", "heat_pump_sources"),
        script:
            scripts("plot_heat_source_map.py")


if config["foresight"] == "perfect":

    def output_map_year(w):
        return {
            f"map_{year}": RESULTS
            + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-costs-all_"
            + f"{year}_{w.solver}.pdf"
            for year in config_provider("scenario", "planning_horizons")(w)
        }

    rule plot_power_network_perfect:
        input:
            network=RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years_{solver}.nc",
            regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        output:
            unpack(output_map_year),
        threads: 2
        resources:
            mem_mb=10000,
        params:
            plotting=config_provider("plotting"),
        message:
            "Plotting power network with perfect foresight for {wildcards.clusters} clusters, {wildcards.opts} electric options and {wildcards.sector_opts} sector options"
        script:
            scripts("plot_power_network_perfect.py")


rule make_solver_comparison_elec:
    input:
        networks=lambda w: expand(
            RESULTS + "networks/base_s_{clusters}_elec_{opts}_{solver}.nc",
            clusters=w.clusters,
            opts=w.opts,
            solver=solver_names(w),
            run=config["run"]["name"],
        ),
        benchmarks=lambda w: expand(
            RESULTS + "benchmarks/solve_network/base_s_{clusters}_elec_{opts}_{solver}",
            clusters=w.clusters,
            opts=w.opts,
            solver=solver_names(w),
            run=config["run"]["name"],
        ),
    output:
        summary=RESULTS
        + "csvs/solver_comparison/summary_s_{clusters}_elec_{opts}.csv",
        optimal_capacity=RESULTS
        + "csvs/solver_comparison/optimal_capacity_s_{clusters}_elec_{opts}.csv",
        energy_balance=RESULTS
        + "csvs/solver_comparison/energy_balance_s_{clusters}_elec_{opts}.csv",
        benchmarks=RESULTS
        + "csvs/solver_comparison/benchmarks_s_{clusters}_elec_{opts}.csv",
    log:
        RESULTS + "logs/make_solver_comparison/base_s_{clusters}_elec_{opts}.log",
    benchmark:
        RESULTS + "benchmarks/make_solver_comparison/base_s_{clusters}_elec_{opts}"
    threads: 1
    resources:
        mem_mb=4000,
    params:
        solver_specs=solver_run_specs,
    message:
        "Comparing solved electricity networks across configured solvers for {wildcards.clusters} clusters and {wildcards.opts} electric options"
    script:
        scripts("compare_solver_results.py")


rule make_solver_comparison_sector:
    input:
        networks=lambda w: expand(
            RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.nc",
            clusters=w.clusters,
            opts=w.opts,
            sector_opts=w.sector_opts,
            planning_horizons=w.planning_horizons,
            solver=solver_names(w),
            run=config["run"]["name"],
        ),
        benchmarks=lambda w: expand(
            RESULTS
            + "benchmarks/solve_sector_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}",
            clusters=w.clusters,
            opts=w.opts,
            sector_opts=w.sector_opts,
            planning_horizons=w.planning_horizons,
            solver=solver_names(w),
            run=config["run"]["name"],
        ),
    output:
        summary=RESULTS
        + "csvs/solver_comparison/summary_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        optimal_capacity=RESULTS
        + "csvs/solver_comparison/optimal_capacity_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        energy_balance=RESULTS
        + "csvs/solver_comparison/energy_balance_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        benchmarks=RESULTS
        + "csvs/solver_comparison/benchmarks_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
    log:
        RESULTS
        + "logs/make_solver_comparison/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        RESULTS
        + "benchmarks/make_solver_comparison/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
    threads: 1
    resources:
        mem_mb=4000,
    params:
        solver_specs=solver_run_specs,
    message:
        "Comparing solved sector-coupled networks across configured solvers for {wildcards.clusters} clusters, {wildcards.planning_horizons} planning horizon, {wildcards.opts} electric options and {wildcards.sector_opts} sector options"
    script:
        scripts("compare_solver_results.py")


rule make_solver_comparison_sector_perfect:
    input:
        networks=lambda w: expand(
            RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years_{solver}.nc",
            clusters=w.clusters,
            opts=w.opts,
            sector_opts=w.sector_opts,
            solver=solver_names(w),
            run=config["run"]["name"],
        ),
        benchmarks=lambda w: expand(
            RESULTS
            + "benchmarks/solve_sector_network/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years_{solver}",
            clusters=w.clusters,
            opts=w.opts,
            sector_opts=w.sector_opts,
            solver=solver_names(w),
            run=config["run"]["name"],
        ),
    output:
        summary=RESULTS
        + "csvs/solver_comparison/summary_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.csv",
        optimal_capacity=RESULTS
        + "csvs/solver_comparison/optimal_capacity_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.csv",
        energy_balance=RESULTS
        + "csvs/solver_comparison/energy_balance_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.csv",
        benchmarks=RESULTS
        + "csvs/solver_comparison/benchmarks_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.csv",
    log:
        RESULTS
        + "logs/make_solver_comparison/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years.log",
    benchmark:
        RESULTS
        + "benchmarks/make_solver_comparison/base_s_{clusters}_{opts}_{sector_opts}_brownfield_all_years"
    threads: 1
    resources:
        mem_mb=4000,
    params:
        solver_specs=solver_run_specs,
    message:
        "Comparing solved perfect-foresight sector-coupled networks across configured solvers for {wildcards.clusters} clusters, {wildcards.opts} electric options and {wildcards.sector_opts} sector options"
    script:
        scripts("compare_solver_results.py")


rule make_summary:
    input:
        network=get_network,
    output:
        nodal_costs=RESULTS
        + "csvs/individual/nodal_costs_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        nodal_capacities=RESULTS
        + "csvs/individual/nodal_capacities_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        nodal_capacity_factors=RESULTS
        + "csvs/individual/nodal_capacity_factors_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        capacity_factors=RESULTS
        + "csvs/individual/capacity_factors_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        costs=RESULTS
        + "csvs/individual/costs_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        capacities=RESULTS
        + "csvs/individual/capacities_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        curtailment=RESULTS
        + "csvs/individual/curtailment_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        energy=RESULTS
        + "csvs/individual/energy_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        energy_balance=RESULTS
        + "csvs/individual/energy_balance_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        nodal_energy_balance=RESULTS
        + "csvs/individual/nodal_energy_balance_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        prices=RESULTS
        + "csvs/individual/prices_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        weighted_prices=RESULTS
        + "csvs/individual/weighted_prices_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        market_values=RESULTS
        + "csvs/individual/market_values_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
        metrics=RESULTS
        + "csvs/individual/metrics_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
    log:
        RESULTS
        + "logs/make_summary_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/make_summary_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
        )
    threads: 1
    resources:
        mem_mb=8000,
    message:
        "Creating optimization results summary statistics"
    script:
        scripts("make_summary.py")


rule make_global_summary:
    input:
        nodal_costs=expand(
            RESULTS
            + "csvs/individual/nodal_costs_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        nodal_capacities=expand(
            RESULTS
            + "csvs/individual/nodal_capacities_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        nodal_capacity_factors=expand(
            RESULTS
            + "csvs/individual/nodal_capacity_factors_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        capacity_factors=expand(
            RESULTS
            + "csvs/individual/capacity_factors_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        costs=expand(
            RESULTS
            + "csvs/individual/costs_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        capacities=expand(
            RESULTS
            + "csvs/individual/capacities_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        curtailment=expand(
            RESULTS
            + "csvs/individual/curtailment_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        energy=expand(
            RESULTS
            + "csvs/individual/energy_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        energy_balance=expand(
            RESULTS
            + "csvs/individual/energy_balance_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        nodal_energy_balance=expand(
            RESULTS
            + "csvs/individual/nodal_energy_balance_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        prices=expand(
            RESULTS
            + "csvs/individual/prices_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        weighted_prices=expand(
            RESULTS
            + "csvs/individual/weighted_prices_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        market_values=expand(
            RESULTS
            + "csvs/individual/market_values_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
        metrics=expand(
            RESULTS
            + "csvs/individual/metrics_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.csv",
            **config["scenario"],
            solver=solver_names(),
            allow_missing=True,
        ),
    output:
        costs=RESULTS + "csvs/costs.csv",
        capacities=RESULTS + "csvs/capacities.csv",
        energy=RESULTS + "csvs/energy.csv",
        energy_balance=RESULTS + "csvs/energy_balance.csv",
        capacity_factors=RESULTS + "csvs/capacity_factors.csv",
        metrics=RESULTS + "csvs/metrics.csv",
        curtailment=RESULTS + "csvs/curtailment.csv",
        prices=RESULTS + "csvs/prices.csv",
        weighted_prices=RESULTS + "csvs/weighted_prices.csv",
        market_values=RESULTS + "csvs/market_values.csv",
        nodal_costs=RESULTS + "csvs/nodal_costs.csv",
        nodal_capacities=RESULTS + "csvs/nodal_capacities.csv",
        nodal_energy_balance=RESULTS + "csvs/nodal_energy_balance.csv",
        nodal_capacity_factors=RESULTS + "csvs/nodal_capacity_factors.csv",
    log:
        RESULTS + "logs/make_global_summary.log",
    benchmark:
        RESULTS + "benchmarks/make_global_summary"
    threads: 1
    resources:
        mem_mb=8000,
    params:
        scenario=config_provider("scenario"),
        solvers=solver_names,
        RDIR=RDIR,
    message:
        "Creating global summary of optimization results for all scenarios"
    script:
        scripts("make_global_summary.py")


rule make_cumulative_costs:
    input:
        costs=RESULTS + "csvs/costs.csv",
    output:
        cumulative_costs=RESULTS + "csvs/cumulative_costs.csv",
    log:
        RESULTS + "logs/make_cumulative_costs.log",
    benchmark:
        RESULTS + "benchmarks/make_cumulative_costs"
    threads: 1
    resources:
        mem_mb=4000,
    params:
        scenario=config_provider("scenario"),
    message:
        "Calculating cumulative costs over time horizon"
    script:
        scripts("make_cumulative_costs.py")


rule plot_summary:
    input:
        costs=RESULTS + "csvs/costs.csv",
        energy=RESULTS + "csvs/energy.csv",
        balances=RESULTS + "csvs/energy_balance.csv",
        eurostat=resources("eurostat_energy_balances.csv"),
        co2=rules.retrieve_ghg_emissions.output["csv"],
    output:
        costs=RESULTS + "graphs/costs.pdf",
        energy=RESULTS + "graphs/energy.pdf",
        balances=RESULTS + "graphs/balances-energy.pdf",
    log:
        RESULTS + "logs/plot_summary.log",
    threads: 2
    resources:
        mem_mb=10000,
    params:
        countries=config_provider("countries"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        emissions_scope=config_provider("energy", "emissions"),
        plotting=config_provider("plotting"),
        foresight=config_provider("foresight"),
        co2_budget=config_provider("co2_budget"),
        sector=config_provider("sector"),
        RDIR=RDIR,
    message:
        "Plotting summary statistics and results"
    script:
        scripts("plot_summary.py")


rule plot_balance_timeseries:
    input:
        network=get_network,
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/balance_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
        ),
    log:
        RESULTS
        + "logs/plot_balance_timeseries/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.log",
    benchmark:
        RESULTS
        + "benchmarks/plot_balance_timeseries/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
    threads: 16
    resources:
        mem_mb=10000,
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    message:
        "Plotting energy balance time series for {wildcards.clusters} clusters, {wildcards.opts} electric options, {wildcards.sector_opts} sector options and {wildcards.planning_horizons} planning horizons"
    script:
        scripts("plot_balance_timeseries.py")


rule plot_heatmap_timeseries:
    input:
        network=get_network,
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/heatmap_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
        ),
    log:
        RESULTS
        + "logs/plot_heatmap_timeseries/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.log",
    benchmark:
        RESULTS
        + "benchmarks/plot_heatmap_timeseries/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
    threads: 16
    resources:
        mem_mb=10000,
    params:
        plotting=config_provider("plotting"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    message:
        "Plotting heatmap time series visualization for {wildcards.clusters} clusters, {wildcards.opts} electric options, {wildcards.sector_opts} sector options and {wildcards.planning_horizons} planning horizons"
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
    input:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}_{solver}.nc",
    output:
        **{
            f"{plot}_bar": RESULTS
            + f"figures/statistics_{plot}_bar_base_s_{{clusters}}_elec_{{opts}}_{{solver}}.pdf"
            for plot in STATISTICS_BARPLOTS
        },
        barplots_touch=RESULTS
        + "figures/.statistics_plots_base_s_{clusters}_elec_{opts}_{solver}",
    params:
        plotting=config_provider("plotting"),
        barplots=STATISTICS_BARPLOTS,
    message:
        "Plotting base scenario statistics for {wildcards.clusters} clusters and {wildcards.opts} electric options"
    script:
        scripts("plot_statistics.py")


rule build_ambient_air_temperature_yearly_average:
    input:
        cutout=lambda w: input_cutout(w),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        average_ambient_air_temperature=resources(
            "temp_ambient_air_base_s_{clusters}_temporal_aggregate.nc"
        ),
    log:
        RESULTS + "logs/build_ambient_air_temperature_yearly_average/base_s_{clusters}",
    benchmark:
        (
            RESULTS
            + "benchmarks/build_ambient_air_temperature_yearly_average/base_s_{clusters}"
        )
    threads: 1
    resources:
        mem_mb=5000,
    script:
        scripts("build_ambient_air_temperature_yearly_average.py")


rule plot_cop_profiles:
    input:
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
    output:
        html=RESULTS + "graphs/cop_profiles_s_{clusters}_{planning_horizons}.html",
    log:
        RESULTS + "logs/plot_cop_profiles_s_{clusters}_{planning_horizons}.log",
    benchmark:
        RESULTS + "benchmarks/plot_cop_profiles/s_{clusters}_{planning_horizons}"
    resources:
        mem_mb=10000,
    script:
        scripts("plot_cop_profiles/plot_cop_profiles.py")


rule plot_interactive_bus_balance:
    input:
        network=get_network,
        rc="matplotlibrc",
    output:
        directory=directory(
            RESULTS
            + "graphics/interactive_bus_balance/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
        ),
    log:
        RESULTS
        + "logs/plot_interactive_bus_balance/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}.log",
    benchmark:
        RESULTS
        + "benchmarks/plot_interactive_bus_balance/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{solver}"
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
