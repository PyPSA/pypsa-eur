# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

"""
Production implementation of streamlined PyPSA-EUR workflow compose rules.

This file implements the new 4-step workflow structure:
base → clustered → composed → solved

All configuration is now driven by config sections rather than wildcards.
"""


from pathlib import Path
from scripts._helpers import get_snapshots


def get_compose_inputs(wildcards):
    """Determine inputs for compose rule based on foresight and horizon."""
    temporal = config["temporal"]
    foresight = temporal["foresight"]
    horizon = int(wildcards.horizon)

    # Handle both single value and list for planning_horizons
    planning_horizons = temporal["planning_horizons"]
    if isinstance(planning_horizons, (int, str)):
        horizons = [int(planning_horizons)]
    else:
        horizons = [int(h) for h in planning_horizons]

    # Get cluster count from config using config_provider
    clusters = config_provider(
        "clustering", "cluster_network", "n_clusters", default=50
    )(wildcards)

    inputs = {
        "clustered": resources("networks/clustered.nc"),
        "busmap": resources(f"busmap_elec_s_{clusters}.csv"),
        "tech_costs": resources(f"costs_{horizon}.csv"),
        "load": resources("electricity_demand.csv"),
        "powerplants": resources(f"powerplants_s_{clusters}.csv"),
    }

    # Add renewable profiles
    if config.get("atlite", {}).get("default_cutout", False):
        cutout_year = config["atlite"]["default_cutout"].split("-")[1][:4]
        inputs.update(
            {
                "profile_onwind": resources(f"profile_onwind-{cutout_year}.nc"),
                "profile_offwind-ac": resources(f"profile_offwind-ac-{cutout_year}.nc"),
                "profile_offwind-dc": resources(f"profile_offwind-dc-{cutout_year}.nc"),
                "profile_offwind-float": resources(
                    f"profile_offwind-float-{cutout_year}.nc"
                ),
                "profile_solar": resources(f"profile_solar-{cutout_year}.nc"),
            }
        )
    else:
        # Fallback to generic profiles without year
        inputs.update(
            {
                "profile_onwind": resources("profile_onwind.nc"),
                "profile_offwind-ac": resources("profile_offwind-ac.nc"),
                "profile_offwind-dc": resources("profile_offwind-dc.nc"),
                "profile_offwind-float": resources("profile_offwind-float.nc"),
                "profile_solar": resources("profile_solar.nc"),
            }
        )

    # Add hydro if enabled
    renewable_carriers = config_provider(
        "electricity", "renewable_carriers", default=[]
    )(wildcards)
    if "hydro" in renewable_carriers:
        inputs.update(
            {
                "profile_hydro": resources("profile_hydro.nc"),
                "hydro_capacities": resources("hydro_capacities.csv"),
            }
        )

    # Add conventional inputs
    unit_commitment = config_provider("conventional", "unit_commitment", default=False)(
        wildcards
    )
    if unit_commitment:
        inputs["unit_commitment"] = resources("unit_commitment.csv")

    dynamic_fuel_price = config_provider(
        "conventional", "dynamic_fuel_price", default=False
    )(wildcards)
    if dynamic_fuel_price:
        inputs["fuel_price"] = resources("fuel_price.csv")

    # Add sector-specific inputs if sector coupling is enabled
    sector_enabled = config_provider("sector", "enabled", default=True)(wildcards)
    if sector_enabled:
        inputs.update(
            {
                "clustered_pop_layout": resources(f"pop_layout_elec_s_{clusters}.csv"),
                "pop_weighted_energy_totals": resources(
                    "pop_weighted_energy_totals.csv"
                ),
                "gas_input_nodes_simplified": resources(
                    "gas_input_nodes_simplified.csv"
                ),
                "industrial_demand": resources(f"industrial_demand_{horizon}.csv"),
                "shipping_demand": resources(f"shipping_demand_{horizon}.csv"),
                "aviation_demand": resources(f"aviation_demand_{horizon}.csv"),
                "transport_demand": resources(f"transport_demand_{horizon}.csv"),
                "transport_data": resources("transport_data.csv"),
                "temp_air_total": resources("temp_air_total.nc"),
                "retro_cost": resources("retrofitting_cost_energy.csv"),
                "floor_area": resources("floor_area.csv"),
                "biomass_potentials": resources("biomass_potentials.csv"),
                "biomass_transport_costs": resources("biomass_transport_costs.csv"),
                "co2_price": resources(f"co2_price_{horizon}.csv"),
            }
        )

        # Add heat-related inputs if heating is enabled
        heating_enabled = config_provider("sector", "heating", default=False)(wildcards)
        if heating_enabled:
            inputs.update(
                {
                    "pop_weighted_heat_totals": resources(
                        "pop_weighted_heat_totals.csv"
                    ),
                    "heating_efficiencies": resources("heating_efficiencies.csv"),
                    "cop_profiles": resources(f"cop_profiles_elec_s_{clusters}.nc"),
                    "hourly_heat_demand_total": resources(
                        f"hourly_heat_demand_total_elec_s_{clusters}.nc"
                    ),
                    "district_heat_share": resources(
                        f"district_heat_share_elec_s_{clusters}.csv"
                    ),
                    "solar_thermal_total": (
                        resources(f"solar_thermal_total_elec_s_{clusters}.csv")
                        if config_provider("sector", "solar_thermal", default=False)(
                            wildcards
                        )
                        else []
                    ),
                }
            )

            # Add heat source profiles
            limited_heat_sources = config_provider("limited_heat_sources", default=[])(
                wildcards
            )
            for source in limited_heat_sources:
                inputs[source] = resources(f"{source}_elec_s_{clusters}.csv")

        # Add transport profiles if transport is enabled
        transport_enabled = config_provider("sector", "transport", default=False)(
            wildcards
        )
        if transport_enabled:
            inputs.update(
                {
                    "avail_profile": resources(f"avail_profile_elec_s_{clusters}.csv"),
                    "dsm_profile": resources(f"dsm_profile_elec_s_{clusters}.csv"),
                }
            )

        # Add hydrogen infrastructure if enabled
        h2_network_enabled = config_provider("sector", "H2_network", default=True)(
            wildcards
        )
        if h2_network_enabled:
            inputs["h2_cavern"] = resources("salt_cavern_potentials.csv")
            inputs["clustered_gas_network"] = resources(
                f"gas_network_elec_s_{clusters}.csv"
            )

        # Add CO2 infrastructure if enabled
        co2_network_enabled = config_provider("sector", "co2_network", default=False)(
            wildcards
        )
        if co2_network_enabled:
            inputs["sequestration_potential"] = resources(
                "co2_sequestration_potential.csv"
            )

    # Add brownfield inputs for non-first horizons
    if foresight == "overnight" and len(horizons) > 1:
        raise ValueError(
            "Overnight optimization can only be run for a single planning horizon."
        )

    if horizon != horizons[0]:
        # Not first horizon - need previous network
        prev_horizon = horizons[horizons.index(horizon) - 1]

        if foresight == "myopic":
            # Myopic uses solved network from previous horizon
            inputs["network_previous"] = RESULTS + f"networks/solved_{prev_horizon}.nc"
        else:  # perfect foresight
            # Perfect foresight uses composed network from previous horizon
            inputs["network_previous"] = (
                RESULTS + f"networks/composed_{prev_horizon}.nc"
            )

    return inputs


def get_final_horizon():
    """Get the final planning horizon, handling both single values and lists."""
    planning_horizons = config["temporal"]["planning_horizons"]
    if isinstance(planning_horizons, (int, str)):
        return int(planning_horizons)
    else:
        return int(planning_horizons[-1])


def get_solve_inputs(wildcards):
    """Get inputs for solve rule - just the composed network."""
    return {"network": RESULTS + f"networks/composed_{wildcards.horizon}.nc"}


# Main composition rule - combines all network building steps
rule compose_network:
    input:
        unpack(get_compose_inputs),
    output:
        RESULTS + "networks/composed_{horizon}.nc",
    params:
        # Pass entire config sections using config_provider
        temporal=config_provider("temporal"),
        electricity=config_provider("electricity", default={}),
        sector=config_provider("sector", default={}),
        clustering=config_provider("clustering", default={}),
        clustering_temporal=config_provider("clustering", "temporal", default={}),
        existing_capacities=config_provider("existing_capacities", default={}),
        pypsa_eur=config_provider("pypsa_eur", default=["AC", "DC"]),
        renewable=config_provider("renewable", default={}),
        conventional=config_provider("conventional", default={}),
        costs=config_provider("costs", default={}),
        load=config_provider("load", default={}),
        lines=config_provider("lines", default={}),
        links=config_provider("links", default={}),
        industry=config_provider("industry", default={}),
        limited_heat_sources=config_provider("limited_heat_sources", default=[]),
        # Direct parameters
        countries=config_provider("countries", default=[]),
        snapshots=config_provider("snapshots", default={}),
        drop_leap_day=config_provider("enable", "drop_leap_day", default=False),
        # Run identification
        run_name=lambda w: w.run,
        # Derived parameters
        energy_totals_year=config_provider("energy", "energy_totals_year", default=2019),
        # Brownfield settings (for myopic)
        h2_retrofit=config_provider("sector", "H2_retrofit", default=False),
        h2_retrofit_capacity_per_ch4=config_provider(
            "sector", "H2_retrofit_capacity_per_CH4", default=0.6
        ),
        capacity_threshold=config_provider(
            "existing_capacities", "threshold_capacity", default=10
        ),
        # CO2 budget handling
        co2_budget=config_provider("co2_budget", "budget_national", default=None),
        emissions_scope=config_provider(
            "co2_budget", "emissions_scope", default="total"
        ),
    log:
        logs("compose_network_{horizon}.log"),
    threads: 1
    resources:
        mem_mb=20000,
    script:
        "../scripts/compose_network.py"


# Unified solving rule for all foresight approaches in streamlined workflow
rule solve_network_streamlined:
    input:
        unpack(get_solve_inputs),
    output:
        RESULTS + "networks/solved_{horizon}.nc",
        RESULTS + "logs/solver_{horizon}.log",
    params:
        temporal=config_provider("temporal"),
        solver=config_provider("solving", "solver"),
        solver_options=lambda w: config_provider(
            "solving",
            "solver_options",
            config_provider("solving", "solver", "name")(w),
        )(w),
        solving=config_provider("solving"),
    log:
        python=logs("solve_network_{horizon}.log"),
    threads: lambda w: config_provider("solving", "solver", "threads")(w)
    resources:
        mem_mb=lambda w: config_provider("solving", "mem")(w),
    shadow:
        "minimal"
    script:
        "../scripts/solve_network.py"


# Collection rules for different workflow targets
rule prepare_networks_streamlined:
    """Target rule to prepare all composed networks up to final horizon."""
    input:
        lambda w: expand(
            RESULTS + "networks/composed_{horizon}.nc",
            horizon=get_final_horizon(),
        ),


rule solve_networks_streamlined:
    """Target rule to solve all networks up to final horizon."""
    input:
        lambda w: expand(
            RESULTS + "networks/solved_{horizon}.nc",
            horizon=get_final_horizon(),
        ),


# Postprocessing rules
rule make_summary_streamlined:
    """Create summary of solved network."""
    input:
        network=RESULTS + "networks/solved_{horizon}.nc",
        tech_costs=lambda w: resources(f"costs_{w.horizon}.csv"),
    output:
        nodal_costs=RESULTS + "csvs/nodal_costs_{horizon}.csv",
        nodal_capacities=RESULTS + "csvs/nodal_capacities_{horizon}.csv",
        nodal_supply=RESULTS + "csvs/nodal_supply_{horizon}.csv",
        nodal_supply_energy=RESULTS + "csvs/nodal_supply_energy_{horizon}.csv",
        supply=RESULTS + "csvs/supply_{horizon}.csv",
        supply_energy=RESULTS + "csvs/supply_energy_{horizon}.csv",
        prices=RESULTS + "csvs/prices_{horizon}.csv",
        weighted_prices=RESULTS + "csvs/weighted_prices_{horizon}.csv",
        market_values=RESULTS + "csvs/market_values_{horizon}.csv",
        price_statistics=RESULTS + "csvs/price_statistics_{horizon}.csv",
        costs=RESULTS + "csvs/costs_{horizon}.csv",
        capacities=RESULTS + "csvs/capacities_{horizon}.csv",
        curtailment=RESULTS + "csvs/curtailment_{horizon}.csv",
        metrics=RESULTS + "csvs/metrics_{horizon}.csv",
    log:
        logs("make_summary_{horizon}.log"),
    script:
        "../scripts/make_summary.py"


rule plot_network_streamlined:
    """Plot the solved network."""
    input:
        network=RESULTS + "networks/solved_{horizon}.nc",
        regions="resources/regions_onshore_{run}.geojson",
    output:
        map=RESULTS + "plots/network_{horizon}.pdf",
    log:
        logs("plot_network_{horizon}.log"),
    script:
        "../scripts/plot_network.py"


# Validation rule to compare old vs new workflow results
rule validate_results_streamlined:
    """Compare results between old wildcard-based and new config-driven workflow."""
    input:
        old_network="networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",  # Old format
        new_network=RESULTS + "networks/solved_{horizon}.nc",  # New format
    output:
        validation_report=RESULTS + "validation_report_{horizon}.html",
    params:
        clusters=config_provider("clustering", "cluster_network", "n_clusters"),
        opts=config_provider("opts", default=""),
        sector_opts=config_provider("sector_opts", default=""),
        planning_horizons=lambda w: w.horizon,
    log:
        logs("validate_results_{horizon}.log"),
    script:
        "../scripts/validate_results.py"


# No custom resources function needed - use existing path provider
