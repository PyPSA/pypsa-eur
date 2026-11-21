# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Compose the network for one planning horizon from the electricity and
sector-coupling building blocks.

Starting from the clustered base network, the script adds the electricity
components (`add_electricity`), the sector-coupling components if enabled
(`prepare_sector_network`), applies temporal aggregation and the final
preparation for solving (`prepare_network`), and then handles the planning
horizon according to the foresight mode. In the first horizon of a myopic or
perfect foresight run, existing capacities are added (`add_existing_baseyear`).
In later myopic horizons, the optimised capacities of the previous horizon's
solved network are carried over as brownfield (`add_brownfield`). Under
perfect foresight, the horizon is stacked onto the previous composed network
to form a multi-period network (`prepare_perfect_foresight`). Finally, CO2
budget constraints and configured adjustments to costs and potentials are
applied and the network is exported for solving.
"""

import logging
import os

import pandas as pd
import pypsa
import xarray as xr

from scripts._helpers import (
    configure_logging,
    get,
    load_costs,
    sanitize_carriers,
    sanitize_custom_columns,
    set_scenario_config,
    update_p_nom_max,
)
from scripts.add_brownfield import (
    add_brownfield,
    adjust_renewable_capacity_limits,
    adjust_renewable_profiles,
    disable_grid_expansion_if_limit_hit,
    update_dynamic_ptes_capacity,
    update_heat_pump_efficiency,
    update_transmission_limit,
)
from scripts.add_electricity import (
    attach_conventional_generators,
    attach_existing_batteries,
    attach_hydro,
    attach_load,
    attach_renewable_powerplants,
    attach_storageunits,
    attach_stores,
    attach_wind_and_solar,
    drop_sector_carriers,
    estimate_renewable_capacities,
    load_and_aggregate_powerplants,
    sanitize_locations,
    set_transmission_costs,
)
from scripts.add_existing_baseyear import (
    add_build_year_to_new_assets,
    add_heating_capacities_installed_before_baseyear,
    add_power_capacities_installed_before_baseyear,
)
from scripts.definitions.heat_system import HeatSystem
from scripts.prepare_network import (
    add_dynamic_emission_prices,
    add_emission_prices,
    add_gaslimit,
    apply_co2_budget_constraints,
    cap_transmission_capacity,
    maybe_adjust_costs_and_potentials,
    set_transmission_limit,
)
from scripts.prepare_perfect_foresight import (
    add_H2_boilers,
    adjust_stores_for_perfect_foresight,
    apply_investment_period_weightings,
    apply_phase_outs,
    concatenate_network_with_previous,
)
from scripts.prepare_sector_network import (
    add_afforestation,
    add_agriculture,
    add_allam_gas,
    add_ammonia,
    add_aquifer_thermal_energy_storage,
    add_aviation,
    add_biochar,
    add_biogas_upgrading,
    add_biogas_upgrading_cc,
    add_biomass_boilers,
    add_biomass_chp,
    add_biomass_supply,
    add_biomass_to_hydrogen,
    add_biomass_to_liquid,
    add_biomass_to_liquid_cc,
    add_biomass_transport,
    add_biomass_transport_costs,
    add_biosng,
    add_biosng_cc,
    add_carrier_buses,
    add_ccgt_cc,
    add_co2_network,
    add_co2_tracking,
    add_coal_cc,
    add_dac,
    add_electricity_grid_connection,
    add_electrobiofuels,
    add_electrolysis,
    add_enhanced_geothermal,
    add_fischer_tropsch,
    add_gas_boilers,
    add_gas_network,
    add_generation,
    add_h2_fuel_cell,
    add_h2_network,
    add_h2_pipeline_retrofit,
    add_h2_storage,
    add_h2_turbine,
    add_heat_demand,
    add_heat_dsm,
    add_heat_pumps,
    add_heat_vent,
    add_import_options,
    add_industry_ammonia_demand,
    add_industry_biomass_demand,
    add_industry_coal_demand,
    add_industry_electricity_demand,
    add_industry_gas_demand,
    add_industry_heat,
    add_industry_hydrogen_demand,
    add_industry_low_temperature_heat_demand,
    add_industry_methanol_demand,
    add_industry_naphtha_demand,
    add_land_transport,
    add_methanation,
    add_methanol,
    add_methanolisation,
    add_micro_chp,
    add_oil_boilers,
    add_perennials,
    add_process_emissions,
    add_resistive_heaters,
    add_retrofitting,
    add_rock_weathering,
    add_sector_buses,
    add_shipping,
    add_smr,
    add_smr_cc,
    add_solar_thermal,
    add_solid_biomass_import,
    add_unsustainable_biomass,
    add_urban_central_chp,
    add_waste,
    add_waste_heat,
    add_water_pits,
    add_water_tanks,
    build_heat_demand,
    cluster_heat_buses,
    decentral,
    define_spatial,
    heat_system_nodes,
    insert_electricity_distribution_grid,
    insert_gas_distribution_costs,
    load_pop_weighted_energy_totals,
    lossy_bidirectional_links,
    read_transport_corridors,
    reduce_space_heat_demand,
    remove_h2_network,
    remove_lines_without_circuits,
    set_temporal_aggregation,
)

logger = logging.getLogger(__name__)

SECTOR_RENEWABLE_CARRIERS = {"solar rooftop", "afforestation", "biochar", "perennials"}


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

    # Extract configuration and parameters
    params = snakemake.params
    inputs = snakemake.input

    current_horizon = int(snakemake.wildcards.horizon)
    horizons = params.horizons
    foresight = params.foresight
    options = params.sector
    sector_mode = options["enabled"]

    is_first_horizon = current_horizon == horizons[0]

    phase_outs = params.existing_capacities["phase_outs"]
    if phase_outs and foresight != "perfect":
        logger.warning(
            "existing_capacities.phase_outs is only applied under perfect foresight; "
            f"the {len(phase_outs)} configured rule(s) are ignored for foresight "
            f"'{foresight}', where phase-outs follow the powerplant lifetimes instead."
        )

    n = pypsa.Network(snakemake.input.network)
    n_previous = None if is_first_horizon else pypsa.Network(inputs.network_previous)
    nyears = n.snapshot_weightings.objective.sum() / 8760.0
    costs = load_costs(snakemake.input.tech_costs)

    logger.info("Adding electricity components")

    landfall_lengths = {
        tech: settings["landfall_length"]
        for tech, settings in params.renewable.items()
        if "landfall_length" in settings
    }

    if inputs.powerplants:
        ppl = load_and_aggregate_powerplants(
            inputs.powerplants,
            costs,
            consider_efficiency_classes=params.clustering[
                "consider_efficiency_classes"
            ],
            aggregation_strategies=params.clustering["aggregation_strategies"],
            exclude_carriers=params.clustering["exclude_carriers"],
            estimate_efficiencies=params.conventional["estimate_efficiencies"],
        )
    else:
        ppl = pd.DataFrame()

    attach_load(
        n,
        inputs["load"],
        params.load["scaling_factor"],
    )

    set_transmission_costs(
        n,
        costs,
        params.lines["length_factor"],
        params.links["length_factor"],
    )

    renewable_carriers = set(params.electricity["renewable_carriers"])
    extendable_carriers = params.electricity["extendable_carriers"]
    conventional_carriers = params.electricity["conventional_carriers"]
    if params.sector["enabled"]:
        conventional_carriers, extendable_carriers = drop_sector_carriers(
            conventional_carriers, extendable_carriers, foresight
        )
    conventional_inputs = {
        k: v for k, v in inputs.items() if k.startswith("conventional_")
    }

    conventional = params.conventional
    attach_conventional_generators(
        n,
        costs,
        ppl,
        conventional_carriers,
        extendable_carriers,
        renewable_carriers,
        conventional_params=conventional,
        conventional_inputs=conventional_inputs,
        unit_commitment_file=inputs.unit_commitment
        if conventional["unit_commitment"]
        else None,
        fuel_price_file=inputs.fuel_price
        if conventional["dynamic_fuel_price"]
        else None,
    )

    attach_wind_and_solar(
        n,
        costs,
        inputs,
        renewable_carriers,
        extendable_carriers,
        params.lines["length_factor"],
        landfall_lengths,
    )

    if "hydro" in renewable_carriers:
        attach_hydro(
            n,
            costs,
            ppl,
            inputs.profile_hydro,
            inputs.hydro_capacities,
            **params.renewable["hydro"],
        )

    estimate_renewable_caps = params.electricity["estimate_renewable_capacities"]
    if estimate_renewable_caps["enable"]:
        if foresight != "overnight":
            logger.info(
                "Skipping renewable capacity estimation because they are added later "
                "in add_existing_baseyear with foresight mode 'myopic'."
            )
        else:
            tech_map = estimate_renewable_caps["technology_mapping"]
            tech_map = {k: v for k, v in tech_map.items() if v in renewable_carriers}
            expansion_limit = estimate_renewable_caps["expansion_limit"]
            year = estimate_renewable_caps["year"]

            if estimate_renewable_caps["from_powerplantmatching"]:
                attach_renewable_powerplants(n, tech_map, inputs)

            if estimate_renewable_caps["from_irenastat"]:
                estimate_renewable_capacities(
                    n, year, tech_map, expansion_limit, params.countries
                )

    update_p_nom_max(n)

    max_hours = params.electricity["max_hours"]
    attach_storageunits(
        n, costs, n.buses.index, extendable_carriers["StorageUnit"], max_hours
    )
    attach_stores(n, costs, n.buses.index, extendable_carriers["Store"])

    if params.electricity["estimate_battery_capacities"]:
        attach_existing_batteries(n, costs, ppl)

    logger.info("Completed electricity components")

    if sector_mode:
        logger.info("Adding sector components")

        cf_industry = params.industry

        pop_layout = pd.read_csv(inputs.clustered_pop_layout, index_col=0)

        pop_weighted_energy_totals = load_pop_weighted_energy_totals(
            inputs.pop_weighted_energy_totals,
            inputs.pop_weighted_heat_totals if options["heating"] else None,
            nyears,
        )

        fn = inputs.gas_input_nodes_simplified
        gas_input_nodes = pd.read_csv(fn, index_col=0) if fn else None

        fn = inputs.heating_efficiencies
        year = int(params["energy_totals_year"])
        heating_efficiencies = pd.read_csv(fn, index_col=[1, 0]).loc[year]

        spatial = define_spatial(
            pop_layout.index,
            options,
            inputs.district_heat_share if options["heating"] else None,
        )

        corridors = read_transport_corridors(
            inputs.transmission_candidates,
            params.lines["length_factor"],
            options["transport_corridors"],
        )

        if foresight in ["myopic", "perfect"]:
            for carrier in params.existing_capacities["conventional_carriers"]:
                add_carrier_buses(
                    n=n,
                    carrier=carrier,
                    costs=costs,
                    spatial=spatial,
                    options=options,
                    cf_industry=cf_industry,
                )

        add_sector_buses(n, spatial, options)

        emission_prices = params.costs["emission_prices"]
        co2_price = (
            get(emission_prices["co2"], current_horizon)
            if emission_prices["enable"]
            else 0.0
        )
        add_co2_tracking(
            n,
            costs,
            options,
            spatial,
            sequestration_potential_file=inputs["sequestration_potential"],
            co2_price=co2_price,
            co2_liquefaction=options["co2_network_liquefaction"],
        )

        add_generation(
            n=n,
            costs=costs,
            conventionals=options["conventional_generation"],
            spatial=spatial,
            options=options,
            cf_industry=cf_industry,
        )

        logger.info("Add hydrogen storage")
        if options["electrolysis"]:
            add_electrolysis(n, costs, spatial, options)

        if options["hydrogen_fuel_cell"]:
            add_h2_fuel_cell(n, costs, spatial)

        if options["hydrogen_turbine"]:
            add_h2_turbine(n, costs, spatial)

        add_h2_storage(
            n,
            costs,
            spatial,
            options,
            inputs.h2_cavern,
            options["hydrogen_underground_storage_locations"],
        )

        if options["gas_network"]:
            add_gas_network(
                n,
                costs,
                spatial,
                options,
                inputs.clustered_gas_network,
                gas_input_nodes,
            )

        if options["H2_retrofit"]:
            add_h2_pipeline_retrofit(
                n, costs, spatial, options, inputs.clustered_gas_network
            )

        if options["H2_network"]:
            add_h2_network(n, costs, corridors["H2"])

        if options["methanation"]:
            add_methanation(n, costs, spatial, options)

        if options["coal_cc"]:
            add_coal_cc(n, costs, spatial, options)

        if options["ccgt_cc"]:
            add_ccgt_cc(n, costs, spatial, options)

        if options["SMR_cc"]:
            add_smr_cc(n, costs, spatial, options)

        if options["SMR"]:
            add_smr(n, costs, spatial)

        if options["transport"]:
            add_land_transport(
                n=n,
                costs=costs,
                transport_demand_file=inputs.transport_demand,
                transport_data_file=inputs.transport_data,
                avail_profile_file=inputs.avail_profile,
                dsm_profile_file=inputs.dsm_profile,
                temp_air_total_file=inputs.temp_air_total,
                cf_industry=cf_industry,
                options=options,
                spatial=spatial,
                current_horizon=current_horizon,
                nodes=spatial.nodes,
            )

        if options["heating"]:
            ates_capex_as_fraction_of_geothermal_heat_source = options[
                "district_heating"
            ]["ates"]["capex_as_fraction_of_geothermal_heat_source"]
            ates_marginal_cost_charger = options["district_heating"]["ates"][
                "marginal_cost_charger"
            ]
            ates_recovery_factor = options["district_heating"]["ates"][
                "recovery_factor"
            ]
            enable_ates = options["district_heating"]["ates"]["enable"]
            heat_source_profile_files = {
                source: inputs[source]
                for source in options["district_heating"]["limited_heat_sources"]
                if source in inputs.keys()  # noqa: SIM118
            }
            logger.info("Add heat sector")

            enable_ptes = options["district_heating"]["ptes"]["enable"]

            heat_demand = build_heat_demand(
                n,
                inputs.hourly_heat_demand_total,
                pop_weighted_energy_totals,
                heating_efficiencies,
            )

            district_heat_info = pd.read_csv(inputs.district_heat_share, index_col=0)
            dist_fraction = district_heat_info["district fraction of node"]
            urban_fraction = district_heat_info["urban fraction"]

            # NB: must add costs of central heating afterwards (EUR 400 / kWpeak, 50a, 1% FOM from Fraunhofer ISE)

            if options["reduce_space_heat_exogenously"]:
                reduce_space_heat_demand(heat_demand, options, current_horizon)

            for heat_system in (
                HeatSystem
            ):  # this loops through all heat systems defined in _entities.HeatSystem
                overdim_factor = options["overdimension_heat_generators"][
                    heat_system.central_or_decentral
                ]
                heat_nodes = heat_system_nodes(heat_system, spatial)

                # if heat_system == HeatSystem.URBAN_CENTRAL and options["central_heat_vent"]:
                if options["heat_vent"][heat_system.system_type.value]:
                    add_heat_vent(n, heat_nodes, heat_system, options)

                add_heat_demand(
                    n,
                    heat_nodes,
                    heat_system,
                    options,
                    heat_demand,
                    urban_fraction,
                    dist_fraction,
                )

                if options["residential_heat"]["dsm"]["enable"] and heat_system in [
                    HeatSystem.RESIDENTIAL_RURAL,
                    HeatSystem.RESIDENTIAL_URBAN_DECENTRAL,
                    HeatSystem.URBAN_CENTRAL,
                ]:
                    add_heat_dsm(
                        n,
                        costs,
                        heat_nodes,
                        heat_system,
                        options,
                        heat_demand,
                        urban_fraction,
                        dist_fraction,
                        inputs.heat_dsm_profile,
                        current_horizon,
                    )

                if options["ttes"]:
                    add_water_tanks(n, costs, heat_nodes, heat_system)

                if enable_ptes and heat_system == HeatSystem.URBAN_CENTRAL:
                    add_water_pits(
                        n,
                        costs,
                        heat_nodes,
                        heat_system,
                        options,
                        inputs.ptes_direct_utilisation_profiles,
                        inputs.ptes_e_max_pu_profiles,
                    )

                if enable_ates and heat_system == HeatSystem.URBAN_CENTRAL:
                    add_aquifer_thermal_energy_storage(
                        n,
                        costs,
                        heat_nodes,
                        heat_system,
                        inputs.ates_potentials,
                        ates_capex_as_fraction_of_geothermal_heat_source,
                        ates_recovery_factor,
                        ates_marginal_cost_charger,
                    )

                add_heat_pumps(
                    n,
                    costs,
                    heat_nodes,
                    heat_system,
                    options,
                    inputs.cop_profiles,
                    inputs.direct_heat_source_utilisation_profiles,
                    heat_source_profile_files,
                    overdim_factor,
                    enable_ptes,
                )

                if options["resistive_heaters"]:
                    add_resistive_heaters(
                        n, costs, heat_nodes, heat_system, overdim_factor
                    )

                if options["boilers"]:
                    add_gas_boilers(
                        n, costs, heat_nodes, heat_system, spatial, overdim_factor
                    )

                if options["solar_thermal"]:
                    add_solar_thermal(
                        n,
                        costs,
                        heat_nodes,
                        heat_system,
                        options,
                        inputs.solar_thermal_total,
                        overdim_factor,
                    )

                if options["chp"]["enable"] and heat_system == HeatSystem.URBAN_CENTRAL:
                    add_urban_central_chp(n, costs, heat_nodes, spatial, options)

                if (
                    options["chp"]["enable"]
                    and options["chp"]["micro_chp"]
                    and heat_system.value != "urban central"
                ):
                    add_micro_chp(n, costs, heat_nodes, heat_system, spatial)

            if options["retrofitting"]["retro_endogen"]:
                add_retrofitting(
                    n,
                    pop_layout,
                    options,
                    heat_demand,
                    dist_fraction,
                    urban_fraction,
                    inputs.retro_cost,
                    inputs.floor_area,
                )

        if options["biomass"]:
            logger.info("Add biomass")

            biomass_potentials = (
                pd.read_csv(inputs.biomass_potentials, index_col=0) * nyears
            )

            add_biomass_supply(n, costs, spatial, options, biomass_potentials)

            if options["solid_biomass_import"]["enable"]:
                add_solid_biomass_import(n, costs, spatial, options, nyears)

            add_unsustainable_biomass(
                n, costs, spatial, options, cf_industry, biomass_potentials
            )

            if options["biogas_upgrading"]:
                add_biogas_upgrading(n, costs, spatial)

            if options["biogas_upgrading_cc"]:
                add_biogas_upgrading_cc(n, costs, spatial)

            if options["biomass_transport"]:
                add_biomass_transport(
                    n, corridors["biomass"], inputs.biomass_transport_costs
                )

            elif options["biomass_spatial"]:
                add_biomass_transport_costs(
                    n,
                    costs,
                    spatial,
                    biomass_potentials,
                    inputs.biomass_transport_costs,
                )

            if options["chp"]["enable"] and "solid biomass" in options["chp"]["fuel"]:
                add_biomass_chp(n, costs, spatial, options)

            if options["biomass_boiler"]:
                add_biomass_boilers(n, costs, spatial, options)

            # Solid biomass to liquid fuel
            if options["biomass_to_liquid"]:
                add_biomass_to_liquid(n, costs, spatial, options, cf_industry)

            # Solid biomass to liquid fuel with carbon capture
            if options["biomass_to_liquid_cc"]:
                add_biomass_to_liquid_cc(n, costs, spatial)

            # Electrobiofuels (BtL with hydrogen addition to make more use of biogenic carbon).
            # Combination of efuels and biomass to liquid, both based on Fischer-Tropsch.
            # Experimental version - use with caution
            if options["electrobiofuels"]:
                add_electrobiofuels(n, costs, spatial, options, cf_industry)

            # BioSNG from solid biomass
            if options["biosng"]:
                add_biosng(n, costs, spatial)

            # BioSNG from solid biomass with carbon capture
            if options["biosng_cc"]:
                add_biosng_cc(n, costs, spatial)

            if options["bioH2"]:
                add_biomass_to_hydrogen(n, costs, spatial, options)

        if cf_industry["waste_to_energy"] or cf_industry["waste_to_energy_cc"]:
            add_waste(
                n=n,
                costs=costs,
                options=options,
                cf_industry=cf_industry,
                spatial=spatial,
                pop_layout=pop_layout,
                current_horizon=current_horizon,
                industrial_demand_file=inputs.industrial_demand,
                biomass_potentials_file=inputs.biomass_potentials,
                biomass_transport_costs_file=inputs.biomass_transport_costs,
                transport_corridors=corridors["biomass"],
                nyears=nyears,
            )

        if options["ammonia"]:
            add_ammonia(n, costs, spatial, cf_industry)

        if options["methanol"]:
            add_methanol(n, costs, options=options, spatial=spatial)

        if options["industry"]:
            logger.info("Add industrial demand")

            # Ensure the gas carrier bus exists before adding any gas-for-industry links.
            add_carrier_buses(
                n=n,
                carrier="gas",
                costs=costs,
                spatial=spatial,
                options=options,
                cf_industry=None,
            )

            # add oil buses for shipping, aviation and naptha for industry
            add_carrier_buses(
                n,
                carrier="oil",
                costs=costs,
                spatial=spatial,
                options=options,
                cf_industry=cf_industry,
            )
            add_carrier_buses(
                n,
                carrier="methanol",
                costs=costs,
                spatial=spatial,
                options=options,
                cf_industry=cf_industry,
            )

            nhours = n.snapshot_weightings.generators.sum()

            # 1e6 to convert TWh to MWh
            industrial_demand = (
                pd.read_csv(inputs.industrial_demand, index_col=0) * 1e6 * nyears
            )

            if not options["biomass"]:
                raise ValueError(
                    "Industry demand includes solid biomass, but `sector.biomass` is disabled. "
                    "Enable `sector: {biomass: true}` in config."
                )

            if options["industry_heat"]["enable"]:
                add_industry_heat(
                    n,
                    spatial.nodes,
                    industrial_demand / nhours,
                    costs,
                    options,
                    cf_industry,
                    spatial,
                )

            add_industry_biomass_demand(
                n, costs, spatial, options, cf_industry, industrial_demand, nhours
            )
            add_industry_gas_demand(
                n, costs, spatial, options, cf_industry, industrial_demand, nhours
            )
            add_industry_hydrogen_demand(n, spatial, industrial_demand, nhours)
            add_industry_methanol_demand(
                n, costs, spatial, options, industrial_demand, nhours
            )
            add_methanolisation(n, costs, spatial, options)

            if options["oil_boilers"]:
                add_oil_boilers(n, costs, spatial, options)

            add_fischer_tropsch(n, costs, spatial, options)
            add_industry_naphtha_demand(
                n,
                costs,
                spatial,
                options,
                cf_industry,
                industrial_demand,
                nhours,
                current_horizon,
            )
            add_industry_low_temperature_heat_demand(
                n, spatial, industrial_demand, nhours
            )
            add_industry_electricity_demand(
                n,
                spatial,
                cf_industry,
                industrial_demand,
                nhours,
                inputs.industrial_electricity_profiles,
            )
            add_process_emissions(
                n, costs, spatial, options, cf_industry, industrial_demand, nhours
            )

            if options["ammonia"]:
                add_industry_ammonia_demand(
                    n, spatial, options, industrial_demand, nhours
                )

            add_industry_coal_demand(
                n, costs, spatial, options, cf_industry, industrial_demand, nhours
            )

        if options["shipping"]:
            add_shipping(
                n=n,
                costs=costs,
                shipping_demand_file=inputs.shipping_demand,
                pop_weighted_energy_totals=pop_weighted_energy_totals,
                options=options,
                spatial=spatial,
                current_horizon=current_horizon,
            )

        if options["aviation"]:
            add_aviation(
                n=n,
                costs=costs,
                pop_weighted_energy_totals=pop_weighted_energy_totals,
                options=options,
                spatial=spatial,
            )

        if options["heating"]:
            add_waste_heat(n, costs, options, cf_industry)

        if options["agriculture"]:
            add_agriculture(
                n,
                costs,
                pop_weighted_energy_totals,
                current_horizon,
                options,
                spatial,
            )

        if options["dac"]:
            add_dac(n, costs, spatial)

        if options["rock_weathering"]["enable"]:
            add_rock_weathering(
                n,
                costs,
                spatial,
                options,
                inputs.carbon_dioxide_removal_land_rock_weathering,
            )

        if options["biochar"]["enable"]:
            add_biochar(
                n, costs, spatial, options, inputs.carbon_dioxide_removal_land_biochar
            )

        if options["afforestation"]["enable"]:
            add_afforestation(
                n, costs, spatial, options, inputs.afforestation_potentials, nyears
            )

        if options["perennials"]["enable"]:
            add_perennials(n, costs, spatial, options, inputs.biomass_potentials)

        if not options["electricity_transmission_grid"]:
            decentral(n)

        if not options["H2_network"]:
            remove_h2_network(n)

        if options["co2_network"]:
            add_co2_network(
                n,
                costs,
                corridors["CO2"],
                co2_network_cost_factor=options["co2_network_cost_factor"],
                co2_liquefaction=options["co2_network_liquefaction"],
            )

        if options["allam_cycle_gas"]:
            add_allam_gas(n, costs, spatial=spatial)

        if options["electricity_distribution_grid"]:
            insert_electricity_distribution_grid(
                n, costs, options, inputs.solar_rooftop_potentials
            )

        if options["enhanced_geothermal"]["enable"]:
            logger.info("Adding Enhanced Geothermal Systems (EGS).")
            add_enhanced_geothermal(
                n,
                costs=costs,
                costs_config=params.costs,
                egs_potentials=inputs.egs_potentials,
                egs_overlap=inputs.egs_overlap,
                egs_config=options["enhanced_geothermal"],
                spatial=spatial,
                egs_capacity_factors=inputs.egs_capacity_factors,
            )

        if options["imports"]["enable"]:
            add_import_options(n, costs, options, spatial, gas_input_nodes)

        if options["gas_distribution_grid"]:
            insert_gas_distribution_costs(n, costs, options=options)

        if options["electricity_grid_connection"]:
            add_electricity_grid_connection(n, costs)

        for k, v in options["transmission_efficiency"].items():
            if k in options["transmission_efficiency"]["enable"]:
                lossy_bidirectional_links(n, k, v)

        if params.transmission_losses:
            remove_lines_without_circuits(n)

        logger.info("Completed sector components")

    logger.info("Applying temporal aggregation")
    n = set_temporal_aggregation(
        n, params.clustering["temporal"], inputs.snapshot_weightings
    )
    logger.info("Completed temporal aggregation")

    capacity_limited_carriers = sorted(
        set(params.electricity["renewable_carriers"])
        | (
            SECTOR_RENEWABLE_CARRIERS
            & (set(n.generators.carrier) | set(n.links.carrier))
        )
    )

    if foresight != "overnight" and is_first_horizon:
        logger.info("Adding existing capacities")
        baseyear = params.horizons[0]

        # define spatial resolution of carriers
        spatial = define_spatial(n.buses[n.buses.carrier == "AC"].index, options)
        add_build_year_to_new_assets(n, baseyear)

        grouping_years_power = params.existing_capacities["grouping_years_power"]
        grouping_years_heat = params.existing_capacities["grouping_years_heat"]
        add_power_capacities_installed_before_baseyear(
            n=n,
            costs=costs,
            grouping_years=grouping_years_power,
            baseyear=baseyear,
            powerplants_file=inputs.powerplants,
            countries=params.countries,
            capacity_threshold=params.existing_capacities["threshold_capacity"],
            renewable_carriers=params.electricity["renewable_carriers"],
            spatial=spatial,
            solar_rooftop_ratio=params.existing_capacities["solar_rooftop_ratio"],
        )

        if options["heating"]:
            add_heating_capacities_installed_before_baseyear(
                n=n,
                costs=costs,
                baseyear=baseyear,
                grouping_years=grouping_years_heat,
                heat_pump_cop=xr.open_dataarray(inputs.cop_profiles),
                use_time_dependent_cop=options["time_dep_hp_cop"],
                default_lifetime=params.existing_capacities["default_heating_lifetime"],
                existing_capacities=pd.read_csv(
                    inputs.existing_heating_distribution,
                    header=[0, 1],
                    index_col=0,
                ),
                heat_pump_source_types=options["heat_pump_sources"],
                efficiency_file=inputs.heating_efficiencies,
                energy_totals_year=params["energy_totals_year"],
                capacity_threshold=params.existing_capacities["threshold_capacity"],
                use_electricity_distribution_grid=options[
                    "electricity_distribution_grid"
                ],
                spatial=spatial,
            )
        if foresight == "myopic":
            adjust_renewable_capacity_limits(
                n, str(current_horizon), capacity_limited_carriers
            )

    if sector_mode and params.sector["cluster_heat_buses"]:
        cluster_heat_buses(n)

    logger.info("Preparing network for solving")

    electricity_cfg = params.electricity

    if electricity_cfg["gaslimit_enable"]:
        gaslimit = electricity_cfg["gaslimit"]
        add_gaslimit(n, gaslimit, nyears)

    emission_prices = params.costs["emission_prices"]

    if emission_prices["dynamic"]:
        if not os.path.exists(inputs.co2_price):
            raise ValueError("CO2 price file for monthly prices not found")
        add_dynamic_emission_prices(n, inputs.co2_price)
    elif emission_prices["enable"]:
        if isinstance(emission_prices["co2"], dict):
            logger.warning(
                "Not setting emission prices specified per planning horizon. "
                "Use dynamic emission prices instead."
            )
        else:
            add_emission_prices(n, {"co2": emission_prices["co2"]}, exclude_co2=False)

    transmission_limit = electricity_cfg["transmission_limit"]
    if isinstance(transmission_limit, dict) and params.foresight == "perfect":
        raise ValueError(
            "Per-horizon `transmission_limit` is not supported for perfect foresight."
        )
    transmission_limit = get(transmission_limit, current_horizon)
    set_transmission_limit(n, transmission_limit)

    cap_transmission_capacity(
        n,
        line_max=params.lines["s_nom_max"],
        link_max=params.links["p_nom_max"],
        line_max_extension=params.lines["s_nom_max_extension"],
        link_max_extension=params.links["p_nom_max_extension"],
        line_max_pu=params.lines["s_max_pu"],
        link_max_pu=params.links["p_max_pu"],
    )

    if foresight == "myopic" and not is_first_horizon:
        adjust_renewable_profiles(n, inputs, params, current_horizon)

        add_build_year_to_new_assets(n, current_horizon)

        update_heat_pump_efficiency(n, n_previous, current_horizon)

        ptes = params.sector["district_heating"]["ptes"]
        if ptes["enable"] and ptes["dynamic_capacity"]:
            update_dynamic_ptes_capacity(n, n_previous, current_horizon)

        logger.info(
            f"Applying brownfield constraints from horizon {horizons[horizons.index(current_horizon) - 1]}"
        )

        add_brownfield(
            n,
            n_previous,
            current_horizon,
            h2_retrofit=params.sector["H2_retrofit"],
            h2_retrofit_capacity_per_ch4=params.sector["H2_retrofit_capacity_per_CH4"],
            capacity_threshold=params.existing_capacities["threshold_capacity"],
        )
        update_transmission_limit(
            n, get(params.electricity["transmission_limit"], current_horizon)
        )
        disable_grid_expansion_if_limit_hit(n)

        adjust_renewable_capacity_limits(
            n, str(current_horizon), capacity_limited_carriers
        )

    if foresight == "perfect":
        logger.warning(
            "Running perfect foresight is not properly tested and may not work as expected. "
            "Use at your own risk!"
        )

        add_build_year_to_new_assets(n, current_horizon)

        add_H2_boilers(n)

        adjust_stores_for_perfect_foresight(n)

        logger.info(f"Converting horizon {current_horizon} to multi-period structure")
        n.set_investment_periods([current_horizon])

        if not is_first_horizon:
            logger.info(
                "Concatenating with previous composed network for perfect foresight"
            )
            n = concatenate_network_with_previous(
                n_previous=n_previous,
                n_current=n,
                current_horizon=current_horizon,
            )

            social_discountrate = params.costs["social_discountrate"]
            apply_investment_period_weightings(n, social_discountrate)
            logger.info(
                f"Investment period weightings applied for {len(n.investment_periods)} periods"
            )
        apply_phase_outs(n, phase_outs, horizons)

    apply_co2_budget_constraints(
        n, inputs=inputs, params=params, nyears=nyears, current_horizon=current_horizon
    )

    adjustments = params.adjustments
    maybe_adjust_costs_and_potentials(n, adjustments["electricity"], current_horizon)
    maybe_adjust_costs_and_potentials(n, adjustments["sector"], current_horizon)

    sanitize_custom_columns(n)
    sanitize_carriers(n, params.nice_names, params.tech_colors)
    sanitize_locations(n)
    n.meta["snakemake_config"] = dict(
        snakemake.config, wildcards=dict(snakemake.wildcards)
    )
    if n.investment_periods.empty:
        n.consistency_check()
    else:
        logger.warning(
            "Skipping consistency check for multi-period (perfect foresight) network: "
            "known pandas/PyPSA MultiIndex limitation."
        )

    logger.info(f"Exporting composed network for horizon {current_horizon}")
    n.export_to_netcdf(snakemake.output[0])
