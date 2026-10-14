<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Industry

The rules on this page build the material-based energy demand of industry
bottom-up: energy intensities per subsector today and after transformation,
material production per country, and its distribution to regions. What is
modelled is described under [Industry](../design/industry.md) in the design
section.

## Energy intensities

{{ rules("build_industry_sector_ratios", "build_industry_sector_ratios_intermediate") }}

## Production

{{ rules("build_ammonia_production", "build_industrial_production_per_country", "build_industrial_production_per_country_tomorrow") }}

## Regional demand

{{ rules("build_industrial_distribution_key", "build_industrial_production_per_node", "build_industrial_energy_demand_per_country_today", "build_industrial_energy_demand_per_node_today", "build_industrial_energy_demand_per_node") }}
