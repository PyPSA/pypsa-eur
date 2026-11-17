<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Heating

The rules on this page prepare the heat demand of buildings, the district
heating shares and areas, the temperature-dependent performance of heat pumps
and heat sources, thermal storage and the cost of building renovation. What is
modelled is described under [Heating](../design/heat.md) in the design section.

## Demand and district heating

{{ rules("build_temperature_profiles", "build_daily_heat_demand", "build_hourly_heat_demand", "build_district_heat_share", "build_dh_areas", "build_retro_cost", "build_existing_heating_distribution") }}

## Heat pumps and heat sources

{{ rules("build_central_heating_temperature_profiles", "build_cop_profiles", "build_direct_heat_source_utilisation_profiles", "build_geothermal_heat_potential", "build_river_heat_potential", "build_sea_heat_potential", "build_egs_potentials", "build_solar_thermal_profiles") }}

## Thermal storage

{{ rules("build_ptes_operations", "build_ates_potentials") }}
