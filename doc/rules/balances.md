<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Population and Energy Balances

The rules on this page distribute population within countries and derive the
annual energy and emission totals per country and sector from statistical
energy balances. They are the common basis of the sectoral demand rules for
[heating](heating.md), [transport](transport.md) and [industry](industry.md).

## Population layouts

{{ rules("build_population_layouts", "build_clustered_population_layouts") }}

## Energy balances and totals

{{ rules("build_eurostat_balances", "build_swiss_energy_balances", "build_transformation_output_coke", "build_energy_totals", "build_heat_totals", "build_country_hdd", "build_co2_totals", "build_population_weighted_energy_totals") }}
