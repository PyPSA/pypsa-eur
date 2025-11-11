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


def get_compose_inputs(w):
    """Determine inputs for compose rule based on foresight and horizon."""
    cfg = get_config(w)
    foresight = cfg["foresight"]
    horizon = int(w.horizon)

    # Handle both single value and list for planning_horizons
    planning_horizons = cfg["planning_horizons"]
    if isinstance(planning_horizons, (int, str)):
        horizons = [int(planning_horizons)]
    else:
        horizons = [int(h) for h in planning_horizons]

    # Removed unused clusters configuration

    # Start with empty dict and build it up properly
    inputs = {
        **input_profile_tech(w),
        **input_class_regions(w),
        **input_conventional(w),
        **input_profile_offwind(w),
        **input_heat_source_power(w),
        **rules.cluster_gas_network.output,
        **rules.build_gas_input_locations.output,
        "base_network": resources("networks/simplified.nc"),
        "tech_costs": resources(f"costs_{horizon}.csv"),
        "regions": resources("regions_onshore.geojson"),
        "powerplants": resources("powerplants_s.csv"),
        "hydro_capacities": ancient("data/hydro_capacities.csv"),
        "unit_commitment": "data/unit_commitment.csv",
        "fuel_price": (
            resources("monthly_fuel_price.csv")
            if cfg["conventional"]["dynamic_fuel_price"]
            else []
        ),
        "co2_price": resources("co2_price.csv"),
        "load": resources("electricity_demand_simplified.nc"),
        "snapshot_weightings": resources("snapshot_weightings.csv"),
        "retro_cost": (
            resources("retro_cost.csv")
            if cfg["sector"]["retrofitting"]["retro_endogen"]
            else []
        ),
        "floor_area": (
            resources("floor_area.csv")
            if cfg["sector"]["retrofitting"]["retro_endogen"]
            else []
        ),
        "biomass_transport_costs": (
            resources("biomass_transport_costs.csv")
            if cfg["sector"]["biomass_transport"] or cfg["sector"]["biomass_spatial"]
            else []
        ),
        "sequestration_potential": (
            resources("co2_sequestration_potential.csv")
            if cfg["sector"]["regional_co2_sequestration_potential"]["enable"]
            else []
        ),
        "clustered": resources("networks/clustered.nc"),
        "network": resources("networks/clustered.nc"),
        "eurostat": "data/eurostat/Balances-April2023",
        "pop_weighted_energy_totals": resources("pop_weighted_energy_totals.csv"),
        "pop_weighted_heat_totals": resources("pop_weighted_heat_totals.csv"),
        "shipping_demand": resources("shipping_demand.csv"),
        "transport_demand": resources("transport_demand.csv"),
        "transport_data": resources("transport_data.csv"),
        "avail_profile": resources("avail_profile.csv"),
        "dsm_profile": resources("dsm_profile.csv"),
        "co2_totals_name": resources("co2_totals.csv"),
        "co2": "data/bundle/eea/UNFCCC_v23.csv",
        "co2_budget_distribution": resources("co2_budget_distribution.csv"),
        "biomass_potentials": resources("biomass_potentials_{horizon}.csv"),
        "costs": (
            resources("costs_{}.csv".format(cfg["costs"]["year"]))
            if foresight == "overnight"
            else resources("costs_{horizon}.csv")
        ),
        "h2_cavern": resources("salt_cavern_potentials.csv"),
        "busmap_s": resources("busmap_simplified.csv"),
        "busmap": resources("busmap.csv"),
        "clustered_pop_layout": resources("pop_layout.csv"),
        "industrial_demand": resources("industrial_energy_demand_{horizon}.csv"),
        "hourly_heat_demand_total": resources("hourly_heat_demand_total.nc"),
        "industrial_production": resources("industrial_production_{horizon}.csv"),
        "district_heat_share": resources("district_heat_share_{horizon}.csv"),
        "heating_efficiencies": resources("heating_efficiencies.csv"),
        "existing_heating_distribution": resources(
            "existing_heating_distribution_{horizon}.csv"
        ),
        "temp_soil_total": resources("temp_soil_total.nc"),
        "temp_air_total": resources("temp_air_total.nc"),
        "cop_profiles": resources("cop_profiles_{horizon}.nc"),
        "ptes_e_max_pu_profiles": (
            resources("ptes_e_max_pu_profiles_{horizon}.nc")
            if cfg["sector"]["district_heating"]["ptes"]["dynamic_capacity"]
            else []
        ),
        "ptes_direct_utilisation_profiles": (
            resources("ptes_direct_utilisation_profiles_{horizon}.nc")
            if cfg["sector"]["district_heating"]["ptes"]["supplemental_heating"][
                "enable"
            ]
            else []
        ),
        "solar_thermal_total": (
            resources("solar_thermal_total.nc")
            if cfg["sector"]["solar_thermal"]
            else []
        ),
        "solar_rooftop_potentials": (
            resources("solar_rooftop_potentials.csv")
            if "solar" in cfg["electricity"]["renewable_carriers"]
            else []
        ),
        "egs_potentials": (
            resources("egs_potentials.csv")
            if cfg["sector"]["enhanced_geothermal"]["enable"]
            else []
        ),
        "egs_overlap": (
            resources("egs_overlap.csv")
            if cfg["sector"]["enhanced_geothermal"]["enable"]
            else []
        ),
        "egs_capacity_factors": (
            resources("egs_capacity_factors.csv")
            if cfg["sector"]["enhanced_geothermal"]["enable"]
            else []
        ),
        "direct_heat_source_utilisation_profiles": resources(
            "direct_heat_source_utilisation_profiles_{horizon}.nc"
        ),
        "ates_potentials": (
            resources("ates_potentials_{horizon}.csv")
            if cfg["sector"]["district_heating"]["ates"]["enable"]
            else []
        ),
    }

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
        elif foresight == "perfect":
            # Perfect foresight uses composed network from previous horizon
            inputs["network_previous"] = resources(
                f"networks/composed_{prev_horizon}.nc"
            )
        else:
            raise ValueError(f"Invalid foresight type: {foresight}")

    return inputs


# Main composition rule - combines all network building steps
rule compose_network:
    input:
        unpack(get_compose_inputs),
    output:
        resources("networks/composed_{horizon}.nc"),
    params:
        # Pass entire config sections using config_provider
        foresight=config_provider("foresight"),
        electricity=config_provider("electricity"),
        sector=config_provider("sector"),
        clustering=config_provider("clustering"),
        clustering_temporal=config_provider("clustering", "temporal"),
        existing_capacities=config_provider("existing_capacities"),
        pypsa_eur=config_provider("pypsa_eur"),
        renewable=config_provider("renewable"),
        conventional=config_provider("conventional"),
        costs=config_provider("costs"),
        load=config_provider("load"),
        lines=config_provider("lines"),
        links=config_provider("links"),
        industry=config_provider("industry"),
        limited_heat_sources=config_provider(
            "sector", "district_heating", "limited_heat_sources"
        ),
        # Direct parameters
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        # Derived parameters
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        # Parameters from add_existing_baseyear
        baseyear=config_provider("planning_horizons"),
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        # Brownfield settings (for myopic)
        h2_retrofit=config_provider("sector", "H2_retrofit"),
        h2_retrofit_capacity_per_ch4=config_provider(
            "sector", "H2_retrofit_capacity_per_CH4"
        ),
        capacity_threshold=config_provider("existing_capacities", "threshold_capacity"),
        tes=config_provider("sector", "tes"),
        dynamic_ptes_capacity=config_provider(
            "sector", "district_heating", "ptes", "dynamic_capacity"
        ),
        # CO2 budget handling (pass entire co2_budget config section)
        co2_budget=config_provider("co2_budget"),
    log:
        logs("compose_network_{horizon}.log"),
    threads: 1
    resources:
        mem_mb=20000,
    script:
        "../scripts/compose_network.py"
