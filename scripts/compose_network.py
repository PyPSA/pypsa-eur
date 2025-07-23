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

# Import all required functions from existing scripts
from scripts._helpers import (
    configure_logging,
    get_snapshots,
    load_costs,
    sanitize_carriers,
    sanitize_locations,
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
    build_carbon_budget,
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
    temporal_config = config["temporal"]
    foresight = temporal_config["foresight"]
    current_horizon = int(snakemake.wildcards.horizon)

    # Handle both single value and list for planning_horizons
    planning_horizons = temporal_config["planning_horizons"]
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
                h2_retrofit=params.get("h2_retrofit", False),
                h2_retrofit_capacity_per_ch4=params.get("h2_retrofit_capacity_per_ch4"),
                capacity_threshold=params.get("capacity_threshold"),
            )
        elif foresight == "perfect":
            # For perfect foresight: concatenate networks across horizons
            n_previous = pypsa.Network(snakemake.input.network_previous)
            n_current = pypsa.Network(snakemake.input.clustered)

            logger.info(
                f"Concatenating networks for perfect foresight up to horizon {current_horizon}"
            )

            # Create a new network to hold the concatenated result
            n = pypsa.Network()

            # Get all previous horizons from the previous network's snapshots
            if hasattr(n_previous, "investment_periods"):
                previous_horizons = list(n_previous.investment_periods)
            else:
                # First concatenation - previous network only has one horizon
                prev_horizon_idx = horizons.index(current_horizon) - 1
                previous_horizons = [horizons[prev_horizon_idx]]

            all_horizons = previous_horizons + [current_horizon]

            # Process each horizon's network
            networks_to_concat = []

            # Add the previous network(s)
            if len(previous_horizons) == 1:
                # Previous network is a single horizon
                networks_to_concat.append((previous_horizons[0], n_previous))
            else:
                # Previous network already contains multiple horizons
                # We'll extract each period's data during the concatenation
                networks_to_concat.append(("multi", n_previous))

            # Add current horizon's network
            networks_to_concat.append((current_horizon, n_current))

            # Set up multi-indexed snapshots for all horizons
            all_snapshots = []
            for horizon in all_horizons:
                horizon_snapshots = pd.MultiIndex.from_product(
                    [[horizon], n_current.snapshots], names=["period", "timestep"]
                )
                all_snapshots.extend(horizon_snapshots)

            n.set_snapshots(pd.MultiIndex.from_tuples(all_snapshots))

            # Process static components first
            for horizon, network in networks_to_concat:
                if horizon == "multi":
                    # Skip multi-network for now, handle it separately
                    continue

                # Add build year to new assets
                add_build_year_to_new_assets(network, horizon)

                # Copy static components
                for component in network.iterate_components(
                    [
                        "Bus",
                        "Carrier",
                        "Generator",
                        "Link",
                        "Store",
                        "Load",
                        "Line",
                        "StorageUnit",
                    ]
                ):
                    df = component.df.copy()

                    # Get components that don't exist in the target network yet
                    existing_idx = []
                    if hasattr(n, component.list_name):
                        existing_df = getattr(n, component.list_name)
                        existing_idx = existing_df.index

                    new_idx = df.index.difference(existing_idx)
                    if len(new_idx) > 0:
                        n.add(component.name, new_idx, **df.loc[new_idx])

            # Handle multi-horizon previous network
            if networks_to_concat[0][0] == "multi":
                # Copy all components from the multi-horizon network
                multi_network = networks_to_concat[0][1]

                # Static components
                for component in multi_network.iterate_components(
                    [
                        "Bus",
                        "Carrier",
                        "Generator",
                        "Link",
                        "Store",
                        "Load",
                        "Line",
                        "StorageUnit",
                    ]
                ):
                    df = component.df.copy()
                    existing_df = getattr(n, component.list_name)
                    new_idx = df.index.difference(existing_df.index)
                    if len(new_idx) > 0:
                        n.add(component.name, new_idx, **df.loc[new_idx])

                # Time-varying data from multi-horizon network
                for component in multi_network.iterate_components():
                    target_pnl = getattr(n, component.list_name + "_t")
                    source_pnl = getattr(multi_network, component.list_name + "_t")

                    for attr in source_pnl.keys():
                        if not source_pnl[attr].empty:
                            # Initialize the attribute if it doesn't exist
                            if attr not in target_pnl:
                                target_pnl[attr] = pd.DataFrame(index=n.snapshots)

                            # Copy data for previous horizons
                            for prev_horizon in previous_horizons:
                                horizon_data = source_pnl[attr].loc[prev_horizon]
                                target_pnl[attr].loc[prev_horizon] = horizon_data

            # Process time-varying data for new horizon
            network = n_current
            horizon_snapshots = n.snapshots[
                n.snapshots.get_level_values(0) == current_horizon
            ]

            for component in network.iterate_components():
                pnl = getattr(n, component.list_name + "_t")
                source_pnl = component.pnl

                for attr in source_pnl.keys():
                    if not source_pnl[attr].empty:
                        # Initialize if doesn't exist
                        if attr not in pnl:
                            pnl[attr] = pd.DataFrame(index=n.snapshots)

                        # For components from previous periods, extend their time series
                        if current_horizon != all_horizons[0]:
                            # Get components that existed in previous period
                            prev_horizon = all_horizons[
                                all_horizons.index(current_horizon) - 1
                            ]
                            if prev_horizon in pnl[attr].index.get_level_values(0):
                                prev_data = pnl[attr].loc[prev_horizon]
                                # Use previous period's values as starting point
                                pnl[attr].loc[current_horizon] = prev_data.values

                        # Now update with actual data for new/modified components
                        new_data = source_pnl[attr].reindex(network.snapshots)
                        if not new_data.empty:
                            # Map the data to the multi-indexed snapshots
                            for col in new_data.columns:
                                if col in pnl[attr].columns:
                                    pnl[attr].loc[
                                        (current_horizon, slice(None)), col
                                    ] = new_data[col].values

            # Set up snapshot weightings for the new horizon
            if not hasattr(n, "snapshot_weightings") or n.snapshot_weightings.empty:
                n.snapshot_weightings = pd.DataFrame(index=n.snapshots)
                n.snapshot_weightings["objective"] = 0.0
                n.snapshot_weightings["generators"] = 0.0
                n.snapshot_weightings["stores"] = 0.0

            # Add weightings for current horizon
            current_weightings = n_current.snapshot_weightings
            n.snapshot_weightings.loc[current_horizon] = current_weightings.values

            # Set investment periods and weightings
            n.investment_periods = all_horizons

            # Calculate investment period weightings
            if not hasattr(n, "investment_period_weightings"):
                n.investment_period_weightings = pd.DataFrame(index=all_horizons)

            # Simple year-based weightings (can be refined)
            period_years = pd.Series(all_horizons, index=all_horizons)
            time_diff = (
                period_years.diff().shift(-1).fillna(10)
            )  # Assume 10 years for last period
            n.investment_period_weightings["years"] = time_diff

            # Calculate objective weightings with social discount rate
            social_discount_rate = params.costs.get("social_discountrate", 0.01)
            objective_weightings = []
            for i, horizon in enumerate(all_horizons):
                years_from_start = horizon - all_horizons[0]
                weight = 1 / (1 + social_discount_rate) ** years_from_start
                objective_weightings.append(weight * time_diff[horizon])

            n.investment_period_weightings["objective"] = objective_weightings

            logger.info(
                f"Successfully concatenated networks for horizons: {all_horizons}"
            )
        else:  # overnight
            # Should not reach here for overnight with multiple horizons
            n = pypsa.Network(snakemake.input.clustered)

    # Set snapshots
    time = get_snapshots(
        params.get("snapshots", params.get("time")), params.get("drop_leap_day", False)
    )
    n.set_snapshots(time)

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
        params.conventional.get("consider_efficiency_classes", False),
        params.conventional.get("aggregation_strategies", {}),
        params.electricity.get("exclude_carriers", []),
    )

    # Attach load
    attach_load(
        n,
        snakemake.input.load,
        snakemake.input.busmap,
        params.load.get("scaling_factor", 1.0),
    )

    # Set transmission costs
    set_transmission_costs(
        n,
        costs,
        params.lines.get("length_factor", 1.0),
        params.links.get("length_factor", 1.0),
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
    if params.conventional.get("unit_commitment", False):
        unit_commitment = pd.read_csv(snakemake.input.unit_commitment, index_col=0)
    else:
        unit_commitment = None

    # Load dynamic fuel prices if enabled
    if params.conventional.get("dynamic_fuel_price", False):
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
        params.lines.get("length_factor", 1.0),
        landfall_lengths,
    )

    # Attach hydro if included
    if "hydro" in renewable_carriers:
        hydro_params = params.renewable["hydro"].copy()
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
    estimate_renewable_caps = params.electricity.get(
        "estimate_renewable_capacities", {}
    )
    if estimate_renewable_caps.get("enable", False):
        if foresight != "overnight":
            logger.info(
                "Skipping renewable capacity estimation because they are added later "
                "in add_existing_baseyear with foresight mode 'myopic'."
            )
        else:
            tech_map = estimate_renewable_caps["technology_mapping"]
            expansion_limit = estimate_renewable_caps["expansion_limit"]
            year = estimate_renewable_caps["year"]

            if estimate_renewable_caps.get("from_gem", False):
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

    if params.sector.get("enabled", True):
        logger.info("Adding sector components")

        # Load additional data for sectors
        pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
        pop_weighted_energy_totals = (
            pd.read_csv(snakemake.input.pop_weighted_energy_totals, index_col=0)
            * Nyears
        )

        # Check if heat totals input exists
        if hasattr(snakemake.input, "pop_weighted_heat_totals"):
            pop_weighted_heat_totals = (
                pd.read_csv(snakemake.input.pop_weighted_heat_totals, index_col=0)
                * Nyears
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

        # Load heating efficiencies if available
        if hasattr(snakemake.input, "heating_efficiencies"):
            year = int(params.get("energy_totals_year", current_horizon))
            heating_efficiencies = pd.read_csv(
                snakemake.input.heating_efficiencies, index_col=[1, 0]
            ).loc[year]
        else:
            heating_efficiencies = None

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
            sequestration_potential_file=snakemake.input.get("sequestration_potential"),
        )

        # Add generation
        add_generation(
            n=n,
            costs=costs,
            pop_layout=pop_layout,
            conventionals=params.sector.get("conventional_generation", {}),
            spatial=spatial,
            options=params.sector,
            cf_industry=params.industry,
        )

        # Add storage and grids
        add_storage_and_grids(
            n=n,
            costs=costs,
            pop_layout=pop_layout,
            h2_cavern_file=snakemake.input.get("h2_cavern"),
            cavern_types=params.sector.get(
                "hydrogen_underground_storage_locations", []
            ),
            clustered_gas_network_file=snakemake.input.get("clustered_gas_network"),
            gas_input_nodes=gas_input_nodes,
            spatial=spatial,
            options=params.sector,
        )

        # Add transport if enabled
        if params.sector.get("transport", False):
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
        if params.sector.get("heating", False) and heating_efficiencies is not None:
            add_heat(
                n=n,
                costs=costs,
                cop_profiles_file=snakemake.input.cop_profiles,
                direct_heat_source_utilisation_profile_file=snakemake.input.get(
                    "direct_heat_source_utilisation_profiles"
                ),
                hourly_heat_demand_total_file=snakemake.input.hourly_heat_demand_total,
                ptes_e_max_pu_file=snakemake.input.get("ptes_e_max_pu_profiles"),
                ates_e_nom_max=snakemake.input.get("ates_potentials"),
                ates_capex_as_fraction_of_geothermal_heat_source=params.sector.get(
                    "district_heating", {}
                )
                .get("ates", {})
                .get("capex_as_fraction_of_geothermal_heat_source", 0),
                ates_marginal_cost_charger=params.sector.get("district_heating", {})
                .get("ates", {})
                .get("marginal_cost_charger", 0),
                ates_recovery_factor=params.sector.get("district_heating", {})
                .get("ates", {})
                .get("recovery_factor", 0),
                enable_ates=params.sector.get("district_heating", {})
                .get("ates", {})
                .get("enable", False),
                ptes_direct_utilisation_profile=snakemake.input.get(
                    "ptes_direct_utilisation_profiles"
                ),
                district_heat_share_file=snakemake.input.district_heat_share,
                solar_thermal_total_file=snakemake.input.get("solar_thermal_total"),
                retro_cost_file=snakemake.input.retro_cost,
                floor_area_file=snakemake.input.floor_area,
                heat_source_profile_files={
                    source: snakemake.input[source]
                    for source in params.get("limited_heat_sources", [])
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
        if params.sector.get("biomass", False):
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
        if params.sector.get("ammonia", False):
            add_ammonia(n, costs, pop_layout, spatial, params.industry)

        # Add methanol if enabled
        if params.sector.get("methanol", False):
            add_methanol(
                n, costs, options=params.sector, spatial=spatial, pop_layout=pop_layout
            )

        # Add industry if enabled
        if params.sector.get("industry", False):
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
        if params.sector.get("shipping", False):
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
        if params.sector.get("aviation", False):
            add_aviation(
                n=n,
                costs=costs,
                pop_layout=pop_layout,
                pop_weighted_energy_totals=pop_weighted_energy_totals,
                options=params.sector,
                spatial=spatial,
            )

        # Add waste heat if heating is enabled
        if params.sector.get("heating", False):
            add_waste_heat(n, costs, params.sector, params.industry)

        # Add agriculture if enabled (requires H and I)
        if params.sector.get("agriculture", False):
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
        if params.sector.get("dac", False):
            add_dac(n, costs)

        # Remove electricity transmission grid if disabled
        if not params.sector.get("electricity_transmission_grid", True):
            decentral(n)

        # Remove H2 network if disabled
        if not params.sector.get("H2_network", True):
            remove_h2_network(n)

        # Add CO2 network if enabled
        if params.sector.get("co2_network", False):
            add_co2_network(
                n,
                costs,
                co2_network_cost_factor=params.sector.get(
                    "co2_network_cost_factor", 1.0
                ),
            )

        # Add Allam cycle gas if enabled
        if params.sector.get("allam_cycle_gas", False):
            add_allam_gas(n, costs, pop_layout=pop_layout, spatial=spatial)

        # Cluster heat buses if enabled
        if params.sector.get("cluster_heat_buses", False):
            cluster_heat_buses(n)

    # ========== EXISTING CAPACITIES (from add_existing_baseyear.py) ==========

    if params.existing_capacities.get("enabled", False) and is_first_horizon:
        logger.info("Adding existing capacities")

        baseyear = params.existing_capacities.get("baseyear", 2019)

        # Add build year to new assets
        add_build_year_to_new_assets(n, baseyear)

        # Add power capacities installed before baseyear
        add_power_capacities_installed_before_baseyear(
            n,
            baseyear,
            snakemake.input.powerplants,
            params.existing_capacities,
        )

        # Add heating capacities if sector coupling is enabled
        if params.sector.get("enabled", True) and params.sector.get("heating", False):
            add_heating_capacities_installed_before_baseyear(
                n,
                baseyear,
                snakemake.input,
                params.existing_capacities,
            )

    # ========== NETWORK PREPARATION (from prepare_network.py) ==========

    # Apply temporal resolution averaging if specified
    if hasattr(params, "time_resolution") and params.time_resolution:
        offset = params.get("time_resolution_offset", "0h")
        n = average_every_nhours(n, offset, params.get("drop_leap_day", False))

    # Set temporal aggregation for sector components
    if params.sector.get("enabled", True):
        n = set_temporal_aggregation(
            n,
            params.get(
                "time_resolution",
                params.clustering_temporal.get("resolution_sector", "24h"),
            ),
            snakemake.input.get("snapshot_weightings"),
        )

    # Add CO2 limit if specified
    co2limit = params.get("co2limit")
    if co2limit:
        add_co2limit(n, co2limit, Nyears)

    # Add gas limit if specified
    gaslimit = params.get("gaslimit")
    if gaslimit:
        add_gaslimit(n, gaslimit, Nyears)

    # Set transmission limit if specified
    transmission_limit = params.electricity.get("transmission_limit")
    if transmission_limit:
        set_transmission_limit(n, transmission_limit, "volume", costs, Nyears)

    # Add emission prices
    emission_prices = params.costs.get("emission_prices", {})
    if emission_prices:
        add_emission_prices(n, emission_prices, exclude_co2=False)

    # Add dynamic emission prices if available
    if hasattr(snakemake.input, "co2_price") and os.path.exists(
        snakemake.input.co2_price
    ):
        add_dynamic_emission_prices(n, snakemake.input.co2_price)

    # Set line s_max_pu if specified
    s_max_pu = params.lines.get("s_max_pu", 0.7)
    if s_max_pu:
        set_line_s_max_pu(n, s_max_pu)

    # Apply time segmentation if specified
    if hasattr(params, "time_segmentation") and params.time_segmentation:
        n = apply_time_segmentation(
            n,
            params.time_segmentation.get("resolution", "24h"),
            params.time_segmentation.get("segments", []),
        )

    # Handle carbon budget for sector coupling
    if params.sector.get("enabled", True):
        co2_budget = params.get("co2_budget")
        if isinstance(co2_budget, str) and co2_budget.startswith("cb"):
            fn = f"results/{params.run_name}/csvs/carbon_budget_distribution.csv"
            if not os.path.exists(fn):
                emissions_scope = params.get("emissions_scope", "total")
                input_co2 = snakemake.input.get("co2")
                build_carbon_budget(
                    o=fn,
                    emissions_scope=emissions_scope,
                    input_co2=input_co2,
                    fn_industry=snakemake.input.industrial_demand,
                    limit=co2_budget,
                    countries=params.countries,
                )
            co2_budget = pd.read_csv(fn, index_col=0).squeeze()
            co2_budget = co2_budget.loc[str(current_horizon)]

        # Add CO2 limit for sector
        if co2_budget and isinstance(co2_budget, (int, float)):
            add_co2limit(n, co2_budget, Nyears)

    # ========== FINAL CLEANUP ==========

    # Sanitize carriers
    sanitize_carriers(n, config)

    # Sanitize locations if present
    if "location" in n.buses:
        sanitize_locations(n)

    # Add metadata
    n.meta = dict(config, **dict(wildcards=dict(snakemake.wildcards)))

    # Export composed network
    logger.info(f"Exporting composed network for horizon {current_horizon}")
    n.export_to_netcdf(snakemake.output[0])
