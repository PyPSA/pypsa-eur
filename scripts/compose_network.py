# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Compose network by combining all electricity and sector components.

This script combines functionality from add_electricity.py, prepare_sector_network.py,
add_existing_baseyear.py, add_brownfield.py, and prepare_network.py into a unified
composition step for the streamlined workflow.

All function calls are made directly in the main section without defining
additional functions, following the additive approach of importing existing
functions and calling them in sequence.
"""

import logging
import os

import pandas as pd
import pypsa
import xarray as xr

# Import all required functions from existing scripts
from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
    update_p_nom_max,
)
from scripts.add_brownfield import add_brownfield
from scripts.add_electricity import (
    attach_conventional_generators,
    attach_GEM_renewables,
    attach_hydro,
    attach_load,
    attach_storageunits,
    attach_stores,
    attach_wind_and_solar,
    estimate_renewable_capacities,
    load_and_aggregate_powerplants,
    load_costs,
    sanitize_carriers,
    sanitize_locations,
    set_transmission_costs,
)
from scripts.add_existing_baseyear import (
    add_build_year_to_new_assets,
    add_heating_capacities_installed_before_baseyear,
    add_power_capacities_installed_before_baseyear,
)
from scripts.prepare_network import (
    add_co2limit,
    add_dynamic_emission_prices,
    add_emission_prices,
    add_gaslimit,
    apply_time_segmentation,
    average_every_nhours,
    set_line_s_max_pu,
    set_transmission_limit,
)
from scripts.prepare_sector_network import (
    add_agriculture,
    add_allam_gas,
    add_ammonia,
    add_aviation,
    add_biomass,
    add_carrier_buses,
    add_co2_network,
    add_co2_tracking,
    add_dac,
    add_eu_bus,
    add_generation,
    add_heat,
    add_industry,
    add_land_transport,
    add_lifetime_wind_solar,
    add_methanol,
    add_shipping,
    add_storage_and_grids,
    add_waste_heat,
    cluster_heat_buses,
    decentral,
    define_spatial,
    patch_electricity_network,
    remove_h2_network,
    set_temporal_aggregation,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "compose_network",
            run="test-run",
            horizon="2050",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # Extract configuration and parameters
    config = snakemake.config
    params = snakemake.params
    foresight = config["foresight"]
    current_horizon = int(snakemake.wildcards.horizon)

    # Handle both single value and list for planning_horizons
    planning_horizons = config["planning_horizons"]
    if isinstance(planning_horizons, (int, str)):
        horizons = [int(planning_horizons)]
    else:
        horizons = [int(h) for h in planning_horizons]

    # Determine if this is the first horizon
    is_first_horizon = current_horizon == horizons[0]

    # Load appropriate base network based on foresight mode and horizon
    if is_first_horizon:
        # First horizon always starts from clustered network
        n = pypsa.Network(snakemake.input.clustered)
        logger.info(f"Loading clustered network for first horizon {current_horizon}")
    else:
        # Multi-horizon case: handle based on foresight mode
        if foresight == "myopic":
            # For myopic: load previous solved network and apply brownfield constraints
            n_previous = pypsa.Network(snakemake.input.network_previous)
            n = pypsa.Network(snakemake.input.clustered)
            logger.info(
                f"Applying brownfield constraints from horizon {horizons[horizons.index(current_horizon) - 1]}"
            )

            # Apply brownfield constraints
            add_brownfield(
                n,
                n_previous,
                current_horizon,
                h2_retrofit=params["h2_retrofit"],
                h2_retrofit_capacity_per_ch4=params["h2_retrofit_capacity_per_ch4"],
                capacity_threshold=params["capacity_threshold"],
            )
        elif foresight == "perfect":
            logger.info(
                f"Loading clustered network for perfect foresight horizon {current_horizon}"
            )
            n = pypsa.Network(snakemake.input.clustered)

            # Retain reference to the previously composed network for future multi-period coupling
            n.meta["previous_composed_network"] = snakemake.input.network_previous
        else:  # overnight
            # Should not reach here for overnight with multiple horizons
            n = pypsa.Network(snakemake.input.clustered)

    # Calculate year weighting
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0

    # Load costs with proper horizon
    costs = load_costs(
        snakemake.input.tech_costs,
        params.costs,
        params.electricity["max_hours"],
        Nyears,
    )

    # ========== ELECTRICITY COMPONENTS (from add_electricity.py) ==========

    # Load and aggregate powerplants
    ppl = load_and_aggregate_powerplants(
        snakemake.input.powerplants,
        costs,
        consider_efficiency_classes=params.clustering["consider_efficiency_classes"],
        aggregation_strategies=params.clustering["aggregation_strategies"],
        exclude_carriers=params.electricity["exclude_carriers"],
    )

    # Attach load
    attach_load(
        n,
        snakemake.input.load,
        snakemake.input.busmap,
        params.load["scaling_factor"],
    )

    # Set transmission costs
    set_transmission_costs(
        n,
        costs,
        params.lines["length_factor"],
        params.links["length_factor"],
    )

    # Define carrier sets
    renewable_carriers = set(params.electricity["renewable_carriers"])
    extendable_carriers = params.electricity["extendable_carriers"]
    conventional_carriers = params.electricity["conventional_carriers"]

    # Prepare conventional inputs
    conventional_inputs = {
        k: v for k, v in snakemake.input.items() if k.startswith("conventional_")
    }

    # Load unit commitment if enabled
    if params.conventional["unit_commitment"]:
        unit_commitment = pd.read_csv(snakemake.input.unit_commitment, index_col=0)
    else:
        unit_commitment = None

    # Load dynamic fuel prices if enabled
    if params.conventional["dynamic_fuel_price"]:
        fuel_price = pd.read_csv(
            snakemake.input.fuel_price, index_col=0, header=0, parse_dates=True
        )
        fuel_price = fuel_price.reindex(n.snapshots).ffill()
    else:
        fuel_price = None

    # Attach conventional generators
    attach_conventional_generators(
        n,
        costs,
        ppl,
        conventional_carriers,
        extendable_carriers,
        params.conventional,
        conventional_inputs,
        unit_commitment=unit_commitment,
        fuel_price=fuel_price,
    )

    # Prepare landfall lengths for renewables
    landfall_lengths = {
        tech: settings["landfall_length"]
        for tech, settings in params.renewable.items()
        if "landfall_length" in settings.keys()
    }

    # Attach wind and solar
    attach_wind_and_solar(
        n,
        costs,
        snakemake.input,
        renewable_carriers,
        extendable_carriers,
        params.lines["length_factor"],
        landfall_lengths,
    )

    # Attach hydro if included
    if "hydro" in renewable_carriers:
        hydro_params = params.renewable["hydro"]
        carriers = hydro_params.pop("carriers", [])
        attach_hydro(
            n,
            costs,
            ppl,
            snakemake.input.profile_hydro,
            snakemake.input.hydro_capacities,
            carriers,
            **hydro_params,
        )

    # Estimate renewable capacities if enabled and overnight mode
    estimate_renewable_caps = params.electricity["estimate_renewable_capacities"]
    if estimate_renewable_caps["enable"]:
        if foresight != "overnight":
            logger.info(
                "Skipping renewable capacity estimation because they are added later "
                "in add_existing_baseyear with foresight mode 'myopic'."
            )
        else:
            tech_map = estimate_renewable_caps["technology_mapping"]
            expansion_limit = estimate_renewable_caps["expansion_limit"]
            year = estimate_renewable_caps["year"]

            if estimate_renewable_caps["from_gem"]:
                attach_GEM_renewables(n, tech_map, snakemake.input)

            estimate_renewable_capacities(
                n, year, tech_map, expansion_limit, params.countries
            )

    # Update p_nom_max
    update_p_nom_max(n)

    # Attach storage units and stores
    max_hours = params.electricity["max_hours"]
    attach_storageunits(n, costs, extendable_carriers, max_hours)
    attach_stores(n, costs, extendable_carriers)

    # ========== SECTOR COMPONENTS (from prepare_sector_network.py) ==========

    if params.sector["enabled"]:
        logger.info("Adding sector components")

        # Load additional data for sectors
        pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
        pop_weighted_energy_totals = (
            pd.read_csv(snakemake.input.pop_weighted_energy_totals, index_col=0)
            * Nyears
        )

        # Load heat totals
        pop_weighted_heat_totals = (
            pd.read_csv(snakemake.input.pop_weighted_heat_totals, index_col=0) * Nyears
        )
        pop_weighted_energy_totals.update(pop_weighted_heat_totals)

        # Gas input nodes
        gas_input_nodes = pd.read_csv(
            snakemake.input.gas_input_nodes_simplified, index_col=0
        )

        # Profiles for sector components
        profiles = {
            key: snakemake.input[key]
            for key in snakemake.input.keys()
            if key.startswith("profile")
        }

        # Patch electricity network for sector coupling
        carriers_to_keep = params.pypsa_eur
        patch_electricity_network(
            n, costs, carriers_to_keep, profiles, landfall_lengths
        )

        # Load heating efficiencies
        year = int(params["energy_totals_year"])
        heating_efficiencies = pd.read_csv(
            snakemake.input.heating_efficiencies, index_col=[1, 0]
        ).loc[year]

        # Define spatial scope
        spatial = define_spatial(pop_layout.index, params.sector)

        # Add lifetime for wind and solar in myopic/perfect foresight
        if foresight in ["myopic", "perfect"]:
            add_lifetime_wind_solar(n, costs)

            # Add carrier buses for conventional carriers
            for carrier in conventional_carriers:
                add_carrier_buses(
                    n=n,
                    carrier=carrier,
                    costs=costs,
                    spatial=spatial,
                    options=params.sector,
                    cf_industry=params.industry,
                )

        # Add EU bus
        add_eu_bus(n)

        # Add CO2 tracking
        add_co2_tracking(
            n,
            costs,
            params.sector,
            sequestration_potential_file=snakemake.input["sequestration_potential"],
        )

        # Add generation
        add_generation(
            n=n,
            costs=costs,
            pop_layout=pop_layout,
            conventionals=params.sector["conventional_generation"],
            spatial=spatial,
            options=params.sector,
            cf_industry=params.industry,
        )

        # Add storage and grids
        add_storage_and_grids(
            n=n,
            costs=costs,
            pop_layout=pop_layout,
            h2_cavern_file=snakemake.input["h2_cavern"],
            cavern_types=params.sector["hydrogen_underground_storage_locations"],
            clustered_gas_network_file=snakemake.input["clustered_gas_network"],
            gas_input_nodes=gas_input_nodes,
            spatial=spatial,
            options=params.sector,
        )

        # Add transport if enabled
        if params.sector["transport"]:
            add_land_transport(
                n=n,
                costs=costs,
                transport_demand_file=snakemake.input.transport_demand,
                transport_data_file=snakemake.input.transport_data,
                avail_profile_file=snakemake.input.avail_profile,
                dsm_profile_file=snakemake.input.dsm_profile,
                temp_air_total_file=snakemake.input.temp_air_total,
                cf_industry=params.industry,
                options=params.sector,
                investment_year=current_horizon,
                nodes=spatial.nodes,
            )

        # Add heating if enabled
        if params.sector["heating"]:
            add_heat(
                n=n,
                costs=costs,
                cop_profiles_file=snakemake.input.cop_profiles,
                direct_heat_source_utilisation_profile_file=snakemake.input[
                    "direct_heat_source_utilisation_profiles"
                ],
                hourly_heat_demand_total_file=snakemake.input.hourly_heat_demand_total,
                ptes_e_max_pu_file=snakemake.input["ptes_e_max_pu_profiles"],
                ates_e_nom_max=snakemake.input["ates_potentials"],
                ates_capex_as_fraction_of_geothermal_heat_source=params.sector[
                    "district_heating"
                ]["ates"]["capex_as_fraction_of_geothermal_heat_source"],
                ates_marginal_cost_charger=params.sector["district_heating"]["ates"][
                    "marginal_cost_charger"
                ],
                ates_recovery_factor=params.sector["district_heating"]["ates"][
                    "recovery_factor"
                ],
                enable_ates=params.sector["district_heating"]["ates"]["enable"],
                ptes_direct_utilisation_profile=snakemake.input[
                    "ptes_direct_utilisation_profiles"
                ],
                district_heat_share_file=snakemake.input.district_heat_share,
                solar_thermal_total_file=snakemake.input["solar_thermal_total"],
                retro_cost_file=snakemake.input.retro_cost,
                floor_area_file=snakemake.input.floor_area,
                heat_source_profile_files={
                    source: snakemake.input[source]
                    for source in params["limited_heat_sources"]
                    if source in snakemake.input.keys()
                },
                params=params,
                pop_weighted_energy_totals=pop_weighted_energy_totals,
                heating_efficiencies=heating_efficiencies,
                pop_layout=pop_layout,
                spatial=spatial,
                options=params.sector,
                investment_year=current_horizon,
            )

        # Add biomass if enabled
        if params.sector["biomass"]:
            add_biomass(
                n=n,
                costs=costs,
                options=params.sector,
                spatial=spatial,
                cf_industry=params.industry,
                pop_layout=pop_layout,
                biomass_potentials_file=snakemake.input.biomass_potentials,
                biomass_transport_costs_file=snakemake.input.biomass_transport_costs,
                nyears=Nyears,
            )

        # Add ammonia if enabled
        if params.sector["ammonia"]:
            add_ammonia(n, costs, pop_layout, spatial, params.industry)

        # Add methanol if enabled
        if params.sector["methanol"]:
            add_methanol(
                n, costs, options=params.sector, spatial=spatial, pop_layout=pop_layout
            )

        # Add industry if enabled
        if params.sector["industry"]:
            add_industry(
                n=n,
                costs=costs,
                industrial_demand_file=snakemake.input.industrial_demand,
                pop_layout=pop_layout,
                pop_weighted_energy_totals=pop_weighted_energy_totals,
                options=params.sector,
                spatial=spatial,
                cf_industry=params.industry,
                investment_year=current_horizon,
            )

        # Add shipping if enabled
        if params.sector["shipping"]:
            add_shipping(
                n=n,
                costs=costs,
                shipping_demand_file=snakemake.input.shipping_demand,
                pop_layout=pop_layout,
                pop_weighted_energy_totals=pop_weighted_energy_totals,
                options=params.sector,
                spatial=spatial,
                investment_year=current_horizon,
            )

        # Add aviation if enabled
        if params.sector["aviation"]:
            add_aviation(
                n=n,
                costs=costs,
                pop_layout=pop_layout,
                pop_weighted_energy_totals=pop_weighted_energy_totals,
                options=params.sector,
                spatial=spatial,
            )

        # Add waste heat if heating is enabled
        if params.sector["heating"]:
            add_waste_heat(n, costs, params.sector, params.industry)

        # Add agriculture if enabled (requires H and I)
        if params.sector["agriculture"]:
            add_agriculture(
                n,
                costs,
                pop_layout,
                pop_weighted_energy_totals,
                current_horizon,
                params.sector,
                spatial,
            )

        # Add DAC if enabled
        if params.sector["dac"]:
            add_dac(n, costs)

        # Remove electricity transmission grid if disabled
        if not params.sector["electricity_transmission_grid"]:
            decentral(n)

        # Remove H2 network if disabled
        if not params.sector["H2_network"]:
            remove_h2_network(n)

        # Add CO2 network if enabled
        if params.sector["co2_network"]:
            add_co2_network(
                n,
                costs,
                co2_network_cost_factor=params.sector["co2_network_cost_factor"],
            )

        # Add Allam cycle gas if enabled
        if params.sector["allam_cycle_gas"]:
            add_allam_gas(n, costs, pop_layout=pop_layout, spatial=spatial)

        # Cluster heat buses if enabled
        if params.sector["cluster_heat_buses"]:
            cluster_heat_buses(n)

    # ========== EXISTING CAPACITIES (from add_existing_baseyear.py) ==========

    if params.existing_capacities["enabled"] and is_first_horizon:
        logger.info("Adding existing capacities")

        baseyear = params.existing_capacities["baseyear"]
        existing_cfg = params.existing_capacities

        # Add build year to new assets
        add_build_year_to_new_assets(n, baseyear)

        # Add power capacities installed before baseyear
        grouping_years_power = existing_cfg["grouping_years_power"]
        if grouping_years_power:
            add_power_capacities_installed_before_baseyear(
                n=n,
                costs=costs,
                grouping_years=grouping_years_power,
                baseyear=baseyear,
                powerplants_file=snakemake.input.powerplants,
                countries=params.countries,
                capacity_threshold=existing_cfg["threshold_capacity"],
                lifetime_values=params.costs["fill_values"],
                renewable_carriers=list(renewable_carriers),
            )

        # Add heating capacities if sector coupling is enabled
        if params.sector["enabled"] and params.sector["heating"]:
            grouping_years_heat = existing_cfg["grouping_years_heat"]
            if grouping_years_heat:
                existing_heat = pd.read_csv(
                    snakemake.input.existing_heating_distribution,
                    header=[0, 1],
                    index_col=0,
                )
                cop_profiles = xr.open_dataarray(snakemake.input.cop_profiles)

                add_heating_capacities_installed_before_baseyear(
                    n=n,
                    costs=costs,
                    baseyear=baseyear,
                    grouping_years=grouping_years_heat,
                    existing_capacities=existing_heat,
                    heat_pump_cop=cop_profiles,
                    heat_pump_source_types=params.heat_pump_sources,
                    efficiency_file=snakemake.input.heating_efficiencies,
                    use_time_dependent_cop=params.sector["time_dep_hp_cop"],
                    default_lifetime=existing_cfg["default_heating_lifetime"],
                    energy_totals_year=int(params.energy_totals_year),
                    capacity_threshold=existing_cfg["threshold_capacity"],
                    use_electricity_distribution_grid=params.sector[
                        "electricity_distribution_grid"
                    ],
                )

    # ========== NETWORK PREPARATION (from prepare_network.py) ==========

    # Apply temporal resolution averaging if specified
    clustering_temporal_cfg = params.clustering_temporal
    time_resolution_elec = clustering_temporal_cfg["resolution_elec"]
    if isinstance(time_resolution_elec, str) and time_resolution_elec.lower().endswith(
        "h"
    ):
        n = average_every_nhours(n, time_resolution_elec, params.drop_leap_day)

    # Set temporal aggregation for sector components
    if params.sector["enabled"]:
        time_resolution = clustering_temporal_cfg["resolution_sector"]

        if time_resolution:
            n = set_temporal_aggregation(
                n,
                time_resolution,
                snakemake.input.snapshot_weightings,
            )

    # Add CO2 limit if specified
    electricity_cfg = params.electricity
    if electricity_cfg["co2limit_enable"]:
        co2limit = electricity_cfg["co2limit"]
        add_co2limit(n, co2limit, Nyears)

    # Add gas limit if specified
    if electricity_cfg["gaslimit_enable"]:
        gaslimit = electricity_cfg["gaslimit"]
        add_gaslimit(n, gaslimit, Nyears)

    # Set transmission limit if specified
    transmission_limit = electricity_cfg["transmission_limit"]
    if isinstance(transmission_limit, str):
        kind = transmission_limit[0]
        factor = transmission_limit[1:] or "opt"
    elif isinstance(transmission_limit, (list, tuple)) and len(transmission_limit) == 2:
        kind, factor = transmission_limit
    else:
        raise ValueError(
            "transmission_limit must be a string like 'c1.25' or a (kind, factor) pair"
        )

    set_transmission_limit(n, kind, factor, costs, Nyears)

    # Add emission prices
    emission_prices = params.costs["emission_prices"]
    if emission_prices["co2_monthly_prices"]:
        add_dynamic_emission_prices(n, snakemake.input.co2_price)
    elif emission_prices["enable"]:
        add_emission_prices(
            n,
            {"co2": emission_prices["co2"]},
            exclude_co2=False,
        )

    # Add dynamic emission prices if available
    elif os.path.exists(snakemake.input.co2_price):
        add_dynamic_emission_prices(n, snakemake.input.co2_price)

    # Set line s_max_pu if specified
    s_max_pu = params.lines["s_max_pu"]
    if s_max_pu:
        set_line_s_max_pu(n, s_max_pu)

    # Apply time segmentation if specified
    time_segmentation = clustering_temporal_cfg["time_segmentation"]
    if time_segmentation["enable"] and time_segmentation["segments"]:
        n = apply_time_segmentation(
            n,
            time_segmentation["resolution"],
            time_segmentation["segments"],
        )

    # Handle carbon budget for sector coupling
    if params.sector["enabled"]:
        co2_budget = params["co2_budget"]
        if co2_budget and isinstance(co2_budget, str) and co2_budget.startswith("cb"):
            # Skip complex carbon budget building for now in simplified workflow
            logger.info(
                f"Skipping carbon budget calculation for {co2_budget} - not yet implemented in streamlined workflow"
            )
        elif co2_budget and isinstance(co2_budget, (int, float)):
            # Direct CO2 limit
            add_co2limit(n, co2_budget, Nyears)

    # ========== FINAL CLEANUP ==========

    # Sanitize carriers
    sanitize_carriers(n, config)

    # Sanitize locations if present
    if "location" in n.buses:
        sanitize_locations(n)

    if "marginal_cost" in n.generators_t:
        n.generators_t["marginal_cost"] = n.generators_t["marginal_cost"].fillna(
            n.generators.marginal_cost
        )

    # Add metadata
    n.meta = dict(config, **dict(wildcards=dict(snakemake.wildcards)))

    # Export composed network
    logger.info(f"Exporting composed network for horizon {current_horizon}")
    n.export_to_netcdf(snakemake.output[0])
