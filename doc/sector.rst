..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Building Sector-Coupled Networks
##########################################

The preparation process of the sector-coupled version of the PyPSA-Eur energy system model consists of a group of ``snakemake`` rules which are briefly outlined and explained in detail in the sections below.

Not all data dependencies are shipped with the git repository.
Instead we provide separate data bundles which can be obtained
using the ``retrieve*`` rules (:ref:`data`).
Having downloaded the necessary data,

- :mod:`add_brownfield` builds and stores the base network with all buses, HVAC lines and HVDC links, while


Rule ``add_brownfield``
==============================================================================

.. automodule:: add_brownfield

Rule ``add_existing_baseyear``
==============================================================================

.. automodule:: add_existing_baseyear

Rule ``build_existing_heating_distribution``
==============================================================================

.. automodule:: build_existing_heating_distribution


Rule ``build_ammonia_production``
==============================================================================

.. automodule:: build_ammonia_production

Rule ``build_biomass_potentials``
==============================================================================

.. automodule:: build_biomass_potentials

Rule ``build_egs_potentials``
==============================================================================

.. automodule:: build_egs_potentials

Rule ``build_biomass_transport_costs``
==============================================================================

.. automodule:: build_biomass_transport_costs

Rule ``build_clustered_population_layouts``
==============================================================================

.. automodule:: build_clustered_population_layouts

Rule ``build_cop_profiles``
==============================================================================

.. automodule:: build_cop_profiles

Rule ``build_central_heating_temperature_profiles``
==============================================================================

.. automodule:: build_central_heating_temperature_profiles

Rule ``build_energy_totals``
==============================================================================

.. automodule:: build_energy_totals

Rule ``build_heat_totals``
==============================================================================

.. automodule:: build_heat_totals

Rule ``build_gas_input_locations``
==============================================================================

.. automodule:: build_gas_input_locations

Rule ``build_gas_network``
==============================================================================

.. automodule:: build_gas_network

Rule ``build_daily_heat_demand``
==============================================================================

.. automodule:: build_daily_heat_demand

Rule ``build_hourly_heat_demand``
==============================================================================

.. automodule:: build_hourly_heat_demand

Rule ``build_district_heat_share``
==============================================================================

.. automodule:: build_district_heat_share

Rule ``build_industrial_distribution_key``
==============================================================================

.. automodule:: build_industrial_distribution_key

Rule ``build_industrial_energy_demand_per_country_today``
==============================================================================

.. automodule:: build_industrial_energy_demand_per_country_today

Rule ``build_industrial_energy_demand_per_node_today``
==============================================================================

.. automodule:: build_industrial_energy_demand_per_node_today

Rule ``build_industrial_energy_demand_per_node``
==============================================================================

.. automodule:: build_industrial_energy_demand_per_node

Rule ``build_industrial_production_per_country_tomorrow``
==============================================================================

.. automodule:: build_industrial_production_per_country_tomorrow

Rule ``build_industrial_production_per_country``
==============================================================================

.. automodule:: build_industrial_production_per_country

Rule ``build_industrial_production_per_node``
==============================================================================

.. automodule:: build_industrial_production_per_node

Rule ``build_industry_sector_ratios``
==============================================================================

.. automodule:: build_industry_sector_ratios

Rule ``build_population_layouts``
==============================================================================

.. automodule:: build_population_layouts

Rule ``build_population_weighted_energy_totals``
==============================================================================

.. automodule:: build_population_weighted_energy_totals

Rule ``build_retro_cost``
==============================================================================

.. automodule:: build_retro_cost

Rule ``build_salt_cavern_potentials``
==============================================================================

.. automodule:: build_salt_cavern_potentials

Rule ``build_sequestration_potentials``
==============================================================================

.. automodule:: build_sequestration_potentials

Rule ``build_shipping_demand``
==============================================================================

.. automodule:: build_shipping_demand

Rule ``build_solar_thermal_profiles``
==============================================================================

.. automodule:: build_solar_thermal_profiles

Rule ``build_temperature_profiles``
==============================================================================

.. automodule:: build_temperature_profiles

Rule ``build_transport_demand``
==============================================================================

.. automodule:: build_transport_demand

Rule ``cluster_gas_network``
==============================================================================

.. automodule:: cluster_gas_network

Rule ``time_aggregation``
==============================================================================

.. automodule:: time_aggregation

Rule ``prepare_sector_network``
==============================================================================

.. automodule:: prepare_sector_network
