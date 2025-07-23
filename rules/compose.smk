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

    inputs = {
        "clustered": f"networks/{wildcards.run}/clustered.nc",
        "busmap": resources(f"busmap_{wildcards.run}.csv"),
        "tech_costs": resources(f"costs_{horizon}.csv"),
        "load": resources("electricity_demand.csv"),
        "powerplants": resources("powerplants.csv"),
    }

    # Add renewable profiles
    if config.get("atlite", {}).get("default_cutout", False):
        cutout_year = config["atlite"]["default_cutout"].split("-")[1][:4]
        inputs.update(
            {
                "profile_onwind": resources(
                    f"profile_{wildcards.run}/onwind-{cutout_year}.nc"
                ),
                "profile_offwind-ac": resources(
                    f"profile_{wildcards.run}/offwind-ac-{cutout_year}.nc"
                ),
                "profile_offwind-dc": resources(
                    f"profile_{wildcards.run}/offwind-dc-{cutout_year}.nc"
                ),
                "profile_offwind-float": resources(
                    f"profile_{wildcards.run}/offwind-float-{cutout_year}.nc"
                ),
                "profile_solar": resources(
                    f"profile_{wildcards.run}/solar-{cutout_year}.nc"
                ),
            }
        )
    else:
        # Fallback to generic profiles without year
        inputs.update(
            {
                "profile_onwind": resources(f"profile_{wildcards.run}/onwind.nc"),
                "profile_offwind-ac": resources(
                    f"profile_{wildcards.run}/offwind-ac.nc"
                ),
                "profile_offwind-dc": resources(
                    f"profile_{wildcards.run}/offwind-dc.nc"
                ),
                "profile_offwind-float": resources(
                    f"profile_{wildcards.run}/offwind-float.nc"
                ),
                "profile_solar": resources(f"profile_{wildcards.run}/solar.nc"),
            }
        )

    # Add hydro if enabled
    if "hydro" in config.get("electricity", {}).get("renewable_carriers", []):
        inputs.update(
            {
                "profile_hydro": resources(f"profile_{wildcards.run}/hydro.nc"),
                "hydro_capacities": resources("hydro_capacities.csv"),
            }
        )

    # Add conventional inputs
    if config.get("conventional", {}).get("unit_commitment", False):
        inputs["unit_commitment"] = resources("unit_commitment.csv")

    if config.get("conventional", {}).get("dynamic_fuel_price", False):
        inputs["fuel_price"] = resources("fuel_price.csv")

    # Add sector-specific inputs if sector coupling is enabled
    if config.get("sector", {}).get("enabled", True):
        inputs.update(
            {
                "clustered_pop_layout": resources(f"pop_layout_{wildcards.run}.csv"),
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
                "temp_air_total": resources(f"temp_air_total_{wildcards.run}.nc"),
                "retro_cost": resources("retrofitting_cost_energy.csv"),
                "floor_area": resources("floor_area.csv"),
                "biomass_potentials": resources(
                    f"biomass_potentials_{wildcards.run}.csv"
                ),
                "biomass_transport_costs": resources("biomass_transport_costs.csv"),
                "co2_price": resources(f"co2_price_{horizon}.csv"),
            }
        )

        # Add heat-related inputs if heating is enabled
        if config.get("sector", {}).get("heating", False):
            inputs.update(
                {
                    "pop_weighted_heat_totals": resources(
                        "pop_weighted_heat_totals.csv"
                    ),
                    "heating_efficiencies": resources("heating_efficiencies.csv"),
                    "cop_profiles": resources(f"cop_profiles_{wildcards.run}.nc"),
                    "hourly_heat_demand_total": resources(
                        f"hourly_heat_demand_total_{wildcards.run}.nc"
                    ),
                    "district_heat_share": resources(
                        f"district_heat_share_{wildcards.run}.csv"
                    ),
                    "solar_thermal_total": (
                        resources(f"solar_thermal_total_{wildcards.run}.csv")
                        if config.get("sector", {}).get("solar_thermal", False)
                        else []
                    ),
                }
            )

            # Add heat source profiles
            for source in config.get("limited_heat_sources", []):
                inputs[source] = resources(f"{source}_{wildcards.run}.csv")

        # Add transport profiles if transport is enabled
        if config.get("sector", {}).get("transport", False):
            inputs.update(
                {
                    "avail_profile": resources(f"avail_profile_{wildcards.run}.csv"),
                    "dsm_profile": resources(f"dsm_profile_{wildcards.run}.csv"),
                }
            )

        # Add hydrogen infrastructure if enabled
        if config.get("sector", {}).get("H2_network", True):
            inputs["h2_cavern"] = resources("salt_cavern_potentials.csv")
            inputs["clustered_gas_network"] = resources(
                f"gas_network_{wildcards.run}.csv"
            )

        # Add CO2 infrastructure if enabled
        if config.get("sector", {}).get("co2_network", False):
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
            inputs["network_previous"] = (
                f"networks/{wildcards.run}/solved_{prev_horizon}.nc"
            )
        else:  # perfect foresight
            # Perfect foresight uses composed network from previous horizon
            inputs["network_previous"] = (
                f"networks/{wildcards.run}/composed_{prev_horizon}.nc"
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
    return {"network": f"networks/{wildcards.run}/composed_{wildcards.horizon}.nc"}


# Main composition rule - combines all network building steps
rule compose_network:
    input:
        unpack(get_compose_inputs),
    output:
        "networks/{run}/composed_{horizon}.nc",
    params:
        # Pass entire config sections
        temporal=lambda w: config["temporal"],
        electricity=lambda w: config.get("electricity", {}),
        sector=lambda w: config.get("sector", {}),
        clustering=lambda w: config.get("clustering", {}),
        clustering_temporal=lambda w: config.get("clustering", {}).get("temporal", {}),
        existing_capacities=lambda w: config.get("existing_capacities", {}),
        pypsa_eur=lambda w: config.get("pypsa_eur", ["AC", "DC"]),
        renewable=lambda w: config.get("renewable", {}),
        conventional=lambda w: config.get("conventional", {}),
        costs=lambda w: config.get("costs", {}),
        load=lambda w: config.get("load", {}),
        lines=lambda w: config.get("lines", {}),
        links=lambda w: config.get("links", {}),
        industry=lambda w: config.get("industry", {}),
        limited_heat_sources=lambda w: config.get("limited_heat_sources", []),
        # Direct parameters
        countries=lambda w: config.get("countries", []),
        snapshots=lambda w: config.get("snapshots", {}),
        drop_leap_day=lambda w: config.get("enable", {}).get("drop_leap_day", False),
        # Run identification
        run_name=lambda w: wildcards.run,
        # Derived parameters
        energy_totals_year=lambda w: config.get("energy", {}).get(
            "energy_totals_year", 2019
        ),
        # Brownfield settings (for myopic)
        h2_retrofit=lambda w: config.get("sector", {}).get("H2_retrofit", False),
        h2_retrofit_capacity_per_ch4=lambda w: config.get("sector", {}).get(
            "H2_retrofit_capacity_per_CH4", 0.6
        ),
        capacity_threshold=lambda w: config.get("existing_capacities", {}).get(
            "threshold_capacity", 10
        ),
        # CO2 budget handling
        co2_budget=lambda w: config.get("co2_budget", {}).get("budget_national", None),
        emissions_scope=lambda w: config.get("co2_budget", {}).get(
            "emissions_scope", "total"
        ),
    log:
        "logs/{run}/compose_network_{horizon}.log",
    threads: 1
    resources:
        mem_mb=20000,
    script:
        "../scripts/compose_network.py"


# Unified solving rule for all foresight approaches
rule solve_network:
    input:
        unpack(get_solve_inputs),
    output:
        "networks/{run}/solved_{horizon}.nc",
        "logs/{run}/solver_{horizon}.log",
    params:
        temporal=lambda w: config["temporal"],
        solver=lambda w: config["solving"]["solver"],
        solver_options=lambda w: config["solving"]["solver_options"][
            config["solving"]["solver"]["name"]
        ],
        solving=lambda w: config["solving"],
    log:
        python="logs/{run}/solve_network_{horizon}.log",
    threads: config["solving"]["solver"]["threads"]
    resources:
        mem_mb=config["solving"]["mem"],
    shadow:
        "minimal"
    script:
        "../scripts/solve_network.py"


# Collection rules for different workflow targets
rule prepare_networks:
    """Target rule to prepare all composed networks up to final horizon."""
    input:
        lambda w: expand(
            "networks/{run}/composed_{horizon}.nc",
            run=(
                config["run"]["name"]
                if isinstance(config["run"]["name"], str)
                else config["run"]["name"]
            ),
            horizon=get_final_horizon(),
        ),


rule solve_networks:
    """Target rule to solve all networks up to final horizon."""
    input:
        lambda w: expand(
            "networks/{run}/solved_{horizon}.nc",
            run=(
                config["run"]["name"]
                if isinstance(config["run"]["name"], str)
                else config["run"]["name"]
            ),
            horizon=get_final_horizon(),
        ),


# Postprocessing rules
rule make_summary:
    """Create summary of solved network."""
    input:
        network="networks/{run}/solved_{horizon}.nc",
        tech_costs=lambda w: resources(f"costs_{w.horizon}.csv"),
    output:
        nodal_costs="results/{run}/csvs/nodal_costs_{horizon}.csv",
        nodal_capacities="results/{run}/csvs/nodal_capacities_{horizon}.csv",
        nodal_supply="results/{run}/csvs/nodal_supply_{horizon}.csv",
        nodal_supply_energy="results/{run}/csvs/nodal_supply_energy_{horizon}.csv",
        supply="results/{run}/csvs/supply_{horizon}.csv",
        supply_energy="results/{run}/csvs/supply_energy_{horizon}.csv",
        prices="results/{run}/csvs/prices_{horizon}.csv",
        weighted_prices="results/{run}/csvs/weighted_prices_{horizon}.csv",
        market_values="results/{run}/csvs/market_values_{horizon}.csv",
        price_statistics="results/{run}/csvs/price_statistics_{horizon}.csv",
        costs="results/{run}/csvs/costs_{horizon}.csv",
        capacities="results/{run}/csvs/capacities_{horizon}.csv",
        curtailment="results/{run}/csvs/curtailment_{horizon}.csv",
        metrics="results/{run}/csvs/metrics_{horizon}.csv",
    log:
        "logs/{run}/make_summary_{horizon}.log",
    script:
        "../scripts/make_summary.py"


rule plot_network:
    """Plot the solved network."""
    input:
        network="networks/{run}/solved_{horizon}.nc",
        regions="resources/regions_onshore_{run}.geojson",
    output:
        map="results/{run}/plots/network_{horizon}.pdf",
    log:
        "logs/{run}/plot_network_{horizon}.log",
    script:
        "../scripts/plot_network.py"


# Validation rule to compare old vs new workflow results
rule validate_results:
    """Compare results between old wildcard-based and new config-driven workflow."""
    input:
        old_network="networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",  # Old format
        new_network="networks/{run}/solved_{horizon}.nc",  # New format
    output:
        validation_report="results/{run}/validation_report_{horizon}.html",
    params:
        clusters=lambda w: config["clustering"]["cluster_network"]["n_clusters"],
        opts=lambda w: config.get("opts", ""),
        sector_opts=lambda w: config.get("sector_opts", ""),
        planning_horizons=lambda w: w.horizon,
    log:
        "logs/{run}/validate_results_{horizon}.log",
    script:
        "../scripts/validate_results.py"


# Helper functions for resources
def resources(filename):
    """Return path to resource file."""
    run_name = config["run"]["name"]
    if isinstance(run_name, list):
        run_name = run_name[0]  # Use first run for shared resources
    return f"resources/{run_name}/{filename}"
