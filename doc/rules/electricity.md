<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Electricity Demand and Supply

The rules on this page prepare the electricity demand time series, the existing
power plant fleet and the potentials and availability profiles of wind, solar
and hydro power. Their outputs are attached to the clustered network by
[compose_network][].

## Demand

{{ rules("build_electricity_demand_base", "build_electricity_demand", "cluster_electricity_demand") }}

## Power plants

{{ rules("build_powerplants") }}

## Renewable potentials and profiles

{{ rules("build_ship_raster", "build_natura_raster", "determine_availability_matrix", "determine_availability_matrix_MD_UA", "build_renewable_profiles", "build_solar_rooftop_potentials", "build_country_runoff", "build_hydro_profile") }}
