<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!---->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Building Sector-Coupled Networks

The preparation process of the sector-coupled version of the PyPSA-Eur energy system model consists of a group of `snakemake` rules which are briefly outlined and explained in detail in the sections below.

Not all data dependencies are shipped with the git repository.
Instead we provide separate data bundles which can be obtained
using the `retrieve*` rules ([Retrieving Data](retrieve.md)).
Having downloaded the necessary data,

- [add_brownfield][] builds and stores the base network with all buses, HVAC lines and HVDC links, while


## Rule `add_brownfield`

::: add_brownfield

## Rule `add_existing_baseyear`

::: add_existing_baseyear

## Rule `build_existing_heating_distribution`

::: build_existing_heating_distribution


## Rule `build_ammonia_production`

::: build_ammonia_production

## Rule `build_biomass_potentials`

::: build_biomass_potentials

## Rule `build_egs_potentials`

::: build_egs_potentials

## Rule `build_biomass_transport_costs`

::: build_biomass_transport_costs

## Rule `build_clustered_population_layouts`

::: build_clustered_population_layouts

## Rule `build_simplified_population_layouts`

<!-- ::: build_simplified_population_layouts (module not found) -->

## Rule `build_clustered_solar_rooftop_potentials`

::: build_clustered_solar_rooftop_potentials

## Rule `build_cop_profiles`

<!-- ::: build_cop_profiles (directory module, not importable) -->

## Rule `build_direct_heat_source_utilisation_profiles`

::: build_direct_heat_source_utilisation_profiles

## Rule `build_central_heating_temperature_profiles`

<!-- ::: build_central_heating_temperature_profiles (directory module, not importable) -->

## Rule `build_geothermal_heat_potential`

::: build_geothermal_heat_potential

## Rule `build_ates_potentials`

::: build_ates_potentials

## Rule `build_dh_areas`

::: build_dh_areas

## Rule `build_river_heat_potential`

<!-- ::: build_river_heat_potential (module not found) -->

## Rule `build_sea_heat_potential`

<!-- ::: build_sea_heat_potential (module not found) -->

## Rule `build_ptes_operations`

<!-- ::: build_ptes_operations (directory module, not importable) -->

## Rule `build_tes_capacity_profiles`

<!-- ::: build_tes_capacity_profiles (module not found) -->

## Rule `build_eurostat_balances`

::: build_eurostat_balances

## Rule `build_energy_totals`

::: build_energy_totals

## Rule `build_heat_totals`

::: build_heat_totals

## Rule `build_gas_input_locations`

::: build_gas_input_locations

## Rule `build_gas_network`

::: build_gas_network

## Rule `build_daily_heat_demand`

::: build_daily_heat_demand

## Rule `build_hourly_heat_demand`

::: build_hourly_heat_demand

## Rule `build_district_heat_share`

::: build_district_heat_share

## Rule `build_industrial_distribution_key`

::: build_industrial_distribution_key

## Rule `build_industrial_energy_demand_per_country_today`

::: build_industrial_energy_demand_per_country_today

## Rule `build_industrial_energy_demand_per_node_today`

::: build_industrial_energy_demand_per_node_today

## Rule `build_industrial_energy_demand_per_node`

::: build_industrial_energy_demand_per_node

## Rule `build_industrial_production_per_country_tomorrow`

::: build_industrial_production_per_country_tomorrow

## Rule `build_industrial_production_per_country`

::: build_industrial_production_per_country

## Rule `build_industrial_production_per_node`

::: build_industrial_production_per_node

## Rule `build_industry_sector_ratios`

::: build_industry_sector_ratios

## Rule `build_industry_sector_ratios_intermediate`

::: build_industry_sector_ratios_intermediate

## Rule `build_population_layouts`

::: build_population_layouts

## Rule `build_population_weighted_energy_totals`

::: build_population_weighted_energy_totals

## Rule `build_retro_cost`

::: build_retro_cost

## Rule `build_salt_cavern_potentials`

::: build_salt_cavern_potentials

## Rule `build_co2_sequestration_potentials`

::: build_co2_sequestration_potentials

## Rule `build_clustered_co2_sequestration_potentials`

::: build_clustered_co2_sequestration_potentials

## Rule `build_shipping_demand`

::: build_shipping_demand

## Rule `build_solar_thermal_profiles`

::: build_solar_thermal_profiles

## Rule `build_temperature_profiles`

::: build_temperature_profiles

## Rule `build_mobility_profiles`

::: build_mobility_profiles

## Rule `build_transport_demand`

::: build_transport_demand

## Rule `cluster_gas_network`

::: cluster_gas_network

## Rule `time_aggregation`

::: time_aggregation

## Rule `prepare_sector_network`

::: prepare_sector_network
