..
  SPDX-FileCopyrightText: 2023-2024 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _tutorial_sector:

###############################
Tutorial: Sector-Coupled
###############################

.. note::
    If you have not done it yet, follow the :ref:`installation` steps first.

    Also, checkout the tutorial for electricity-only systems first at :ref:`tutorial`.

In this tutorial, we will add further sectors to the electricity-only model from
:ref:`tutorial`, namely industry, transport, and buildings. This
requires processing of a few more raw data sources.

The sector-coupling code can be run as an overnight / greenfield scenario or
with multi-horizon investment with myopic foresight. Pathway analysis with
perfect foresight is under development. See also the documentation on
:ref:`foresight`.

Overnight Scenarios
===========================

Configuration
-------------

The default configuration file (``config/config.default.yaml``) is set up for running
overnight scenarios. Running a sector-coupled model unlocks many further
configuration options. In the example below, we say that the gas network should
be added and spatially resolved. We also say that the existing gas network may
be retrofitted to transport hydrogen instead.

.. literalinclude:: ../config/test/config.overnight.yaml
   :language: yaml
   :start-at: sector:
   :end-before: solving:

Documentation for all options will be added successively to :ref:`config`.

Scenarios can be defined like for electricity-only studies, but with additional
wildcard options.

.. literalinclude:: ../config/test/config.overnight.yaml
   :language: yaml
   :start-at: scenario:
   :end-before: countries:

For allowed wildcard values, refer to :ref:`wildcards`.

Execution
---------

To run an overnight / greenfiled scenario with the specifications above, run

.. code:: bash

    snakemake -call all --configfile config/test/config.overnight.yaml

which will result in the following jobs ``snakemake`` wants to run, some of
which were already included in the electricity-only tutorial:

.. code:: bash

    job                                                 count
    ------------------------------------------------  -------
    add_electricity                                         1
    add_extra_components                                    1
    all                                                     1
    build_ammonia_production                                1
    build_biomass_potentials                                1
    build_clustered_population_layouts                      1
    build_cop_profiles                                      1
    build_daily_heat_demand                                 1
    build_district_heat_share                               1
    build_energy_totals                                     1
    build_gas_input_locations                               1
    build_gas_network                                       1
    build_heat_totals                                       1
    build_hourly_heat_demand                                1
    build_industrial_distribution_key                       1
    build_industrial_energy_demand_per_country_today        1
    build_industrial_energy_demand_per_node                 1
    build_industrial_energy_demand_per_node_today           1
    build_industrial_production_per_country                 1
    build_industrial_production_per_country_tomorrow        1
    build_industrial_production_per_node                    1
    build_industry_sector_ratios                            1
    build_industry_sector_ratios_intermediate               1
    build_population_layouts                                1
    build_population_weighted_energy_totals                 2
    build_renewable_profiles                                4
    build_salt_cavern_potentials                            1
    build_ship_raster                                       1
    build_shipping_demand                                   1
    build_simplified_population_layouts                     1
    build_solar_thermal_profiles                            3
    build_temperature_profiles                              3
    build_transport_demand                                  1
    cluster_gas_network                                     1
    cluster_network                                         1
    make_summary                                            1
    plot_gas_network                                        1
    plot_hydrogen_network                                   1
    plot_power_network                                      1
    plot_power_network_clustered                            1
    plot_summary                                            1
    prepare_network                                         1
    prepare_sector_network                                  1
    retrieve_cost_data                                      1
    retrieve_cutout                                         1
    retrieve_databundle                                     1
    retrieve_electricity_demand                             1
    retrieve_eurostat_data                                  1
    retrieve_gas_infrastructure_data                        1
    simplify_network                                        1
    solve_sector_network                                    1
    total                                                  63

This covers the retrieval of additional raw data from online resources and
preprocessing data about the transport, industry, and heating sectors as well as
additional rules about geological storage and sequestration potentials, gas
infrastructure, and biomass potentials. The collection rule ``all`` will also
generate summary CSV files and plots after the network has been solved
successfully.



.. graphviz::
    :class: full-width
    :align: center

    digraph snakemake_dag {
        graph[bgcolor=white, margin=0];
        node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
        edge[penwidth=2, color=grey];
            0[label = "all", color = "0.44 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.63 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.59 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network", color = "0.47 0.6 0.85", style="rounded"];
            4[label = "prepare_sector_network\nsector_opts: CO2L0-24h-T-H-B-I-A-dist1", color = "0.47 0.6 0.85", style="rounded"];
            5[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.62 0.6 0.85", style="rounded"];
            6[label = "base_network", color = "0.27 0.6 0.85", style="rounded"];
            7[label = "build_shapes", color = "0.43 0.6 0.85", style="rounded"];
            8[label = "retrieve_databundle", color = "0.55 0.6 0.85", style="rounded"];
            9[label = "build_ship_raster", color = "0.29 0.6 0.85", style="rounded"];
            10[label = "retrieve_ship_raster", color = "0.33 0.6 0.85", style="rounded"];
            11[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.06 0.6 0.85", style="rounded"];
            12[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.62 0.6 0.85", style="rounded"];
            13[label = "cluster_gas_network", color = "0.60 0.6 0.85", style="rounded"];
            14[label = "build_gas_network", color = "0.34 0.6 0.85", style="rounded"];
            15[label = "retrieve_gas_infrastructure_data", color = "0.44 0.6 0.85", style="rounded"];
            16[label = "cluster_network\nclusters: 5", color = "0.03 0.6 0.85", style="rounded"];
            17[label = "simplify_network\nsimpl: ", color = "0.64 0.6 0.85", style="rounded"];
            18[label = "add_electricity", color = "0.51 0.6 0.85", style="rounded"];
            19[label = "build_renewable_profiles\ntechnology: solar", color = "0.62 0.6 0.85", style="rounded"];
            20[label = "build_renewable_profiles\ntechnology: onwind", color = "0.62 0.6 0.85", style="rounded"];
            21[label = "retrieve_cost_data\nyear: 2030", color = "0.57 0.6 0.85", style="rounded"];
            22[label = "build_powerplants", color = "0.39 0.6 0.85", style="rounded"];
            23[label = "build_electricity_demand", color = "0.22 0.6 0.85", style="rounded"];
            24[label = "retrieve_electricity_demand", color = "0.09 0.6 0.85", style="rounded"];
            25[label = "retrieve_synthetic_electricity_demand", color = "0.66 0.6 0.85", style="rounded"];
            26[label = "build_gas_input_locations", color = "0.50 0.6 0.85", style="rounded"];
            27[label = "prepare_network\nll: v1.5\nopts: ", color = "0.52 0.6 0.85", style="rounded"];
            28[label = "add_extra_components", color = "0.02 0.6 0.85", style="rounded"];
            29[label = "retrieve_eurostat_data", color = "0.04 0.6 0.85", style="rounded,dashed"];
            30[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.07 0.6 0.85", style="rounded"];
            31[label = "build_energy_totals", color = "0.24 0.6 0.85", style="rounded"];
            32[label = "build_clustered_population_layouts", color = "0.38 0.6 0.85", style="rounded"];
            33[label = "build_population_layouts", color = "0.61 0.6 0.85", style="rounded"];
            34[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.07 0.6 0.85", style="rounded"];
            35[label = "build_heat_totals", color = "0.22 0.6 0.85", style="rounded"];
            36[label = "build_shipping_demand", color = "0.19 0.6 0.85", style="rounded"];
            37[label = "build_transport_demand", color = "0.42 0.6 0.85", style="rounded"];
            38[label = "build_temperature_profiles\nscope: total", color = "0.58 0.6 0.85", style="rounded"];
            39[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.18 0.6 0.85", style="rounded"];
            40[label = "build_salt_cavern_potentials", color = "0.19 0.6 0.85", style="rounded"];
            41[label = "build_simplified_population_layouts", color = "0.37 0.6 0.85", style="rounded"];
            42[label = "build_industrial_energy_demand_per_node", color = "0.57 0.6 0.85", style="rounded"];
            43[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.06 0.6 0.85", style="rounded"];
            44[label = "build_industry_sector_ratios", color = "0.45 0.6 0.85", style="rounded"];
            45[label = "build_ammonia_production", color = "0.42 0.6 0.85", style="rounded"];
            46[label = "build_industrial_energy_demand_per_country_today", color = "0.12 0.6 0.85", style="rounded"];
            47[label = "build_industrial_production_per_country", color = "0.01 0.6 0.85", style="rounded"];
            48[label = "build_industrial_production_per_node", color = "0.30 0.6 0.85", style="rounded"];
            49[label = "build_industrial_distribution_key", color = "0.62 0.6 0.85", style="rounded"];
            50[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.13 0.6 0.85", style="rounded"];
            51[label = "build_industrial_energy_demand_per_node_today", color = "0.29 0.6 0.85", style="rounded"];
            52[label = "build_hourly_heat_demand", color = "0.24 0.6 0.85", style="rounded"];
            53[label = "build_daily_heat_demand\nscope: total", color = "0.27 0.6 0.85", style="rounded"];
            54[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.46 0.6 0.85", style="rounded"];
            55[label = "build_temperature_profiles\nscope: rural", color = "0.58 0.6 0.85", style="rounded"];
            56[label = "build_temperature_profiles\nscope: urban", color = "0.58 0.6 0.85", style="rounded"];
            57[label = "build_cop_profiles", color = "0.65 0.6 0.85", style="rounded"];
            58[label = "build_solar_thermal_profiles\nscope: total", color = "0.23 0.6 0.85", style="rounded"];
            59[label = "build_solar_thermal_profiles\nscope: urban", color = "0.23 0.6 0.85", style="rounded"];
            60[label = "build_solar_thermal_profiles\nscope: rural", color = "0.23 0.6 0.85", style="rounded"];
            61[label = "plot_power_network_clustered", color = "0.49 0.6 0.85", style="rounded"];
            62[label = "plot_power_network", color = "0.05 0.6 0.85", style="rounded"];
            63[label = "plot_hydrogen_network", color = "0.56 0.6 0.85", style="rounded"];
            64[label = "plot_gas_network", color = "0.59 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            29 -> 1
            8 -> 1
            3 -> 2
            21 -> 2
            61 -> 2
            62 -> 2
            63 -> 2
            64 -> 2
            4 -> 3
            5 -> 4
            12 -> 4
            13 -> 4
            26 -> 4
            27 -> 4
            29 -> 4
            30 -> 4
            34 -> 4
            36 -> 4
            37 -> 4
            31 -> 4
            8 -> 4
            39 -> 4
            21 -> 4
            40 -> 4
            17 -> 4
            16 -> 4
            32 -> 4
            41 -> 4
            42 -> 4
            52 -> 4
            54 -> 4
            38 -> 4
            55 -> 4
            56 -> 4
            57 -> 4
            58 -> 4
            59 -> 4
            60 -> 4
            6 -> 5
            8 -> 5
            9 -> 5
            7 -> 5
            11 -> 5
            7 -> 6
            8 -> 7
            10 -> 9
            11 -> 9
            6 -> 12
            8 -> 12
            9 -> 12
            7 -> 12
            11 -> 12
            14 -> 13
            16 -> 13
            15 -> 14
            17 -> 16
            21 -> 16
            18 -> 17
            21 -> 17
            6 -> 17
            19 -> 18
            20 -> 18
            5 -> 18
            12 -> 18
            6 -> 18
            21 -> 18
            22 -> 18
            23 -> 18
            7 -> 18
            6 -> 19
            8 -> 19
            7 -> 19
            11 -> 19
            6 -> 20
            8 -> 20
            7 -> 20
            11 -> 20
            6 -> 22
            24 -> 23
            25 -> 23
            15 -> 26
            16 -> 26
            28 -> 27
            21 -> 27
            16 -> 28
            21 -> 28
            31 -> 30
            32 -> 30
            7 -> 31
            8 -> 31
            29 -> 31
            33 -> 32
            16 -> 32
            11 -> 32
            7 -> 33
            11 -> 33
            35 -> 34
            32 -> 34
            31 -> 35
            7 -> 36
            16 -> 36
            31 -> 36
            32 -> 37
            30 -> 37
            31 -> 37
            8 -> 37
            38 -> 37
            33 -> 38
            16 -> 38
            11 -> 38
            8 -> 39
            16 -> 39
            7 -> 39
            8 -> 40
            16 -> 40
            33 -> 41
            17 -> 41
            11 -> 41
            43 -> 42
            48 -> 42
            51 -> 42
            44 -> 43
            46 -> 43
            47 -> 43
            45 -> 44
            8 -> 44
            8 -> 45
            8 -> 46
            47 -> 46
            45 -> 47
            8 -> 47
            29 -> 47
            49 -> 48
            50 -> 48
            16 -> 49
            32 -> 49
            47 -> 50
            49 -> 51
            46 -> 51
            53 -> 52
            33 -> 53
            16 -> 53
            11 -> 53
            31 -> 54
            32 -> 54
            33 -> 55
            16 -> 55
            11 -> 55
            33 -> 56
            16 -> 56
            11 -> 56
            38 -> 57
            55 -> 57
            56 -> 57
            33 -> 58
            16 -> 58
            11 -> 58
            33 -> 59
            16 -> 59
            11 -> 59
            33 -> 60
            16 -> 60
            11 -> 60
            16 -> 61
            3 -> 62
            16 -> 62
            3 -> 63
            16 -> 63
            3 -> 64
            16 -> 64
    }

|

Myopic Foresight Scenarios
===================================

Configuration
-------------

To activate the myopic foresight mode, set

.. code:: yaml

    foresight: myopic

Scenarios can be defined like for electricity-only studies, but with additional
wildcard options. For the myopic foresight mode, the ``{planning_horizons}`` wildcard
defines the sequence of investment horizons.

.. literalinclude:: ../config/test/config.myopic.yaml
   :language: yaml
   :start-at: scenario:
   :end-before: countries:

For allowed wildcard values, refer to :ref:`wildcards`.

In the myopic foresight mode, you can tweak for instance exogenously given transition paths, like the one for
the share of primary steel production we change below:

.. literalinclude:: ../config/test/config.myopic.yaml
   :language: yaml
   :start-at: industry:
   :end-before: solving:

Documentation for all options will be added successively to :ref:`config`.

Execution
---------

To run a myopic foresight scenario with the specifications above, run

.. code:: bash

    snakemake -call all --configfile config/test/config.myopic.yaml

which will result in additional jobs ``snakemake`` wants to run, which
translates to the following workflow diagram which nicely outlines how the
sequential pathway optimisation with myopic foresight is implemented in the
workflow:

.. graphviz::
    :class: full-width
    :align: center

    digraph snakemake_dag {
        graph[bgcolor=white, margin=0];
        node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
        edge[penwidth=2, color=grey];
            0[label = "all", color = "0.30 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.58 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.56 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network_myopic", color = "0.40 0.6 0.85", style="rounded"];
            4[label = "add_existing_baseyear", color = "0.22 0.6 0.85", style="rounded"];
            5[label = "prepare_sector_network\nsector_opts: 24h-T-H-B-I-A-dist1", color = "0.05 0.6 0.85", style="rounded"];
            6[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.58 0.6 0.85", style="rounded"];
            7[label = "base_network", color = "0.50 0.6 0.85", style="rounded"];
            8[label = "build_shapes", color = "0.11 0.6 0.85", style="rounded"];
            9[label = "retrieve_databundle", color = "0.15 0.6 0.85", style="rounded"];
            10[label = "build_ship_raster", color = "0.07 0.6 0.85", style="rounded"];
            11[label = "retrieve_ship_raster", color = "0.26 0.6 0.85", style="rounded"];
            12[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.31 0.6 0.85", style="rounded"];
            13[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.58 0.6 0.85", style="rounded"];
            14[label = "cluster_gas_network", color = "0.13 0.6 0.85", style="rounded"];
            15[label = "build_gas_network", color = "0.22 0.6 0.85", style="rounded"];
            16[label = "retrieve_gas_infrastructure_data", color = "0.08 0.6 0.85", style="rounded"];
            17[label = "cluster_network\nclusters: 5", color = "0.13 0.6 0.85", style="rounded"];
            18[label = "simplify_network\nsimpl: ", color = "0.60 0.6 0.85", style="rounded"];
            19[label = "add_electricity", color = "0.61 0.6 0.85", style="rounded"];
            20[label = "build_renewable_profiles\ntechnology: solar", color = "0.58 0.6 0.85", style="rounded"];
            21[label = "build_renewable_profiles\ntechnology: onwind", color = "0.58 0.6 0.85", style="rounded"];
            22[label = "retrieve_cost_data\nyear: 2030", color = "0.09 0.6 0.85", style="rounded"];
            23[label = "build_powerplants", color = "0.19 0.6 0.85", style="rounded"];
            24[label = "build_electricity_demand", color = "0.54 0.6 0.85", style="rounded"];
            25[label = "retrieve_electricity_demand", color = "0.18 0.6 0.85", style="rounded"];
            26[label = "retrieve_synthetic_electricity_demand", color = "0.39 0.6 0.85", style="rounded"];
            27[label = "build_gas_input_locations", color = "0.37 0.6 0.85", style="rounded"];
            28[label = "prepare_network\nll: v1.5\nopts: ", color = "0.53 0.6 0.85", style="rounded"];
            29[label = "add_extra_components", color = "0.65 0.6 0.85", style="rounded"];
            30[label = "retrieve_eurostat_data", color = "0.35 0.6 0.85", style="rounded,dashed"];
            31[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.60 0.6 0.85", style="rounded"];
            32[label = "build_energy_totals", color = "0.04 0.6 0.85", style="rounded"];
            33[label = "build_clustered_population_layouts", color = "0.43 0.6 0.85", style="rounded"];
            34[label = "build_population_layouts", color = "0.66 0.6 0.85", style="rounded"];
            35[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.60 0.6 0.85", style="rounded"];
            36[label = "build_heat_totals", color = "0.62 0.6 0.85", style="rounded"];
            37[label = "build_shipping_demand", color = "0.36 0.6 0.85", style="rounded"];
            38[label = "build_transport_demand", color = "0.21 0.6 0.85", style="rounded"];
            39[label = "build_temperature_profiles\nscope: total", color = "0.14 0.6 0.85", style="rounded"];
            40[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.46 0.6 0.85", style="rounded"];
            41[label = "build_salt_cavern_potentials", color = "0.46 0.6 0.85", style="rounded"];
            42[label = "build_simplified_population_layouts", color = "0.55 0.6 0.85", style="rounded"];
            43[label = "build_industrial_energy_demand_per_node", color = "0.32 0.6 0.85", style="rounded"];
            44[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.27 0.6 0.85", style="rounded"];
            45[label = "build_industry_sector_ratios", color = "0.08 0.6 0.85", style="rounded"];
            46[label = "build_ammonia_production", color = "0.49 0.6 0.85", style="rounded"];
            47[label = "build_industrial_energy_demand_per_country_today", color = "0.29 0.6 0.85", style="rounded"];
            48[label = "build_industrial_production_per_country", color = "0.33 0.6 0.85", style="rounded"];
            49[label = "build_industrial_production_per_node", color = "0.34 0.6 0.85", style="rounded"];
            50[label = "build_industrial_distribution_key", color = "0.34 0.6 0.85", style="rounded"];
            51[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.62 0.6 0.85", style="rounded"];
            52[label = "build_industrial_energy_demand_per_node_today", color = "0.11 0.6 0.85", style="rounded"];
            53[label = "build_hourly_heat_demand", color = "0.16 0.6 0.85", style="rounded"];
            54[label = "build_daily_heat_demand\nscope: total", color = "0.28 0.6 0.85", style="rounded"];
            55[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.51 0.6 0.85", style="rounded"];
            56[label = "build_temperature_profiles\nscope: rural", color = "0.14 0.6 0.85", style="rounded"];
            57[label = "build_temperature_profiles\nscope: urban", color = "0.14 0.6 0.85", style="rounded"];
            58[label = "build_cop_profiles", color = "0.51 0.6 0.85", style="rounded"];
            59[label = "build_solar_thermal_profiles\nscope: total", color = "0.52 0.6 0.85", style="rounded"];
            60[label = "build_solar_thermal_profiles\nscope: urban", color = "0.52 0.6 0.85", style="rounded"];
            61[label = "build_solar_thermal_profiles\nscope: rural", color = "0.52 0.6 0.85", style="rounded"];
            62[label = "build_existing_heating_distribution", color = "0.17 0.6 0.85", style="rounded"];
            63[label = "solve_sector_network_myopic", color = "0.40 0.6 0.85", style="rounded"];
            64[label = "add_brownfield", color = "0.00 0.6 0.85", style="rounded"];
            65[label = "prepare_sector_network\nsector_opts: 24h-T-H-B-I-A-dist1", color = "0.05 0.6 0.85", style="rounded"];
            66[label = "build_biomass_potentials\nplanning_horizons: 2040", color = "0.46 0.6 0.85", style="rounded"];
            67[label = "retrieve_cost_data\nyear: 2040", color = "0.09 0.6 0.85", style="rounded"];
            68[label = "build_industrial_energy_demand_per_node", color = "0.32 0.6 0.85", style="rounded"];
            69[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2040", color = "0.27 0.6 0.85", style="rounded"];
            70[label = "build_industrial_production_per_node", color = "0.34 0.6 0.85", style="rounded"];
            71[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2040", color = "0.62 0.6 0.85", style="rounded"];
            72[label = "build_district_heat_share\nplanning_horizons: 2040", color = "0.51 0.6 0.85", style="rounded"];
            73[label = "solve_sector_network_myopic", color = "0.40 0.6 0.85", style="rounded"];
            74[label = "add_brownfield", color = "0.00 0.6 0.85", style="rounded"];
            75[label = "prepare_sector_network\nsector_opts: 24h-T-H-B-I-A-dist1", color = "0.05 0.6 0.85", style="rounded"];
            76[label = "build_biomass_potentials\nplanning_horizons: 2050", color = "0.46 0.6 0.85", style="rounded"];
            77[label = "retrieve_cost_data\nyear: 2050", color = "0.09 0.6 0.85", style="rounded"];
            78[label = "build_industrial_energy_demand_per_node", color = "0.32 0.6 0.85", style="rounded"];
            79[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2050", color = "0.27 0.6 0.85", style="rounded"];
            80[label = "build_industrial_production_per_node", color = "0.34 0.6 0.85", style="rounded"];
            81[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2050", color = "0.62 0.6 0.85", style="rounded"];
            82[label = "build_district_heat_share\nplanning_horizons: 2050", color = "0.51 0.6 0.85", style="rounded"];
            83[label = "plot_power_network_clustered", color = "0.02 0.6 0.85", style="rounded"];
            84[label = "plot_power_network", color = "0.45 0.6 0.85", style="rounded"];
            85[label = "plot_power_network", color = "0.45 0.6 0.85", style="rounded"];
            86[label = "plot_power_network", color = "0.45 0.6 0.85", style="rounded"];
            87[label = "plot_hydrogen_network", color = "0.38 0.6 0.85", style="rounded"];
            88[label = "plot_hydrogen_network", color = "0.38 0.6 0.85", style="rounded"];
            89[label = "plot_hydrogen_network", color = "0.38 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            30 -> 1
            9 -> 1
            3 -> 2
            63 -> 2
            73 -> 2
            22 -> 2
            83 -> 2
            84 -> 2
            85 -> 2
            86 -> 2
            87 -> 2
            88 -> 2
            89 -> 2
            4 -> 3
            22 -> 3
            5 -> 4
            23 -> 4
            18 -> 4
            17 -> 4
            33 -> 4
            22 -> 4
            58 -> 4
            62 -> 4
            6 -> 5
            13 -> 5
            14 -> 5
            27 -> 5
            28 -> 5
            30 -> 5
            31 -> 5
            35 -> 5
            37 -> 5
            38 -> 5
            32 -> 5
            9 -> 5
            40 -> 5
            22 -> 5
            41 -> 5
            18 -> 5
            17 -> 5
            33 -> 5
            42 -> 5
            43 -> 5
            53 -> 5
            55 -> 5
            39 -> 5
            56 -> 5
            57 -> 5
            58 -> 5
            59 -> 5
            60 -> 5
            61 -> 5
            7 -> 6
            9 -> 6
            10 -> 6
            8 -> 6
            12 -> 6
            8 -> 7
            9 -> 8
            11 -> 10
            12 -> 10
            7 -> 13
            9 -> 13
            10 -> 13
            8 -> 13
            12 -> 13
            15 -> 14
            17 -> 14
            16 -> 15
            18 -> 17
            22 -> 17
            19 -> 18
            22 -> 18
            7 -> 18
            20 -> 19
            21 -> 19
            6 -> 19
            13 -> 19
            7 -> 19
            22 -> 19
            23 -> 19
            24 -> 19
            8 -> 19
            7 -> 20
            9 -> 20
            8 -> 20
            12 -> 20
            7 -> 21
            9 -> 21
            8 -> 21
            12 -> 21
            7 -> 23
            25 -> 24
            26 -> 24
            16 -> 27
            17 -> 27
            29 -> 28
            22 -> 28
            17 -> 29
            22 -> 29
            32 -> 31
            33 -> 31
            8 -> 32
            9 -> 32
            30 -> 32
            34 -> 33
            17 -> 33
            12 -> 33
            8 -> 34
            12 -> 34
            36 -> 35
            33 -> 35
            32 -> 36
            8 -> 37
            17 -> 37
            32 -> 37
            33 -> 38
            31 -> 38
            32 -> 38
            9 -> 38
            39 -> 38
            34 -> 39
            17 -> 39
            12 -> 39
            9 -> 40
            17 -> 40
            8 -> 40
            9 -> 41
            17 -> 41
            34 -> 42
            18 -> 42
            12 -> 42
            44 -> 43
            49 -> 43
            52 -> 43
            45 -> 44
            47 -> 44
            48 -> 44
            46 -> 45
            9 -> 45
            9 -> 46
            9 -> 47
            48 -> 47
            46 -> 48
            9 -> 48
            30 -> 48
            50 -> 49
            51 -> 49
            17 -> 50
            33 -> 50
            48 -> 51
            50 -> 52
            47 -> 52
            54 -> 53
            34 -> 54
            17 -> 54
            12 -> 54
            32 -> 55
            33 -> 55
            34 -> 56
            17 -> 56
            12 -> 56
            34 -> 57
            17 -> 57
            12 -> 57
            39 -> 58
            56 -> 58
            57 -> 58
            34 -> 59
            17 -> 59
            12 -> 59
            34 -> 60
            17 -> 60
            12 -> 60
            34 -> 61
            17 -> 61
            12 -> 61
            33 -> 62
            31 -> 62
            55 -> 62
            64 -> 63
            67 -> 63
            20 -> 64
            21 -> 64
            6 -> 64
            13 -> 64
            18 -> 64
            17 -> 64
            65 -> 64
            3 -> 64
            67 -> 64
            58 -> 64
            6 -> 65
            13 -> 65
            14 -> 65
            27 -> 65
            28 -> 65
            30 -> 65
            31 -> 65
            35 -> 65
            37 -> 65
            38 -> 65
            32 -> 65
            9 -> 65
            66 -> 65
            67 -> 65
            41 -> 65
            18 -> 65
            17 -> 65
            33 -> 65
            42 -> 65
            68 -> 65
            53 -> 65
            72 -> 65
            39 -> 65
            56 -> 65
            57 -> 65
            58 -> 65
            59 -> 65
            60 -> 65
            61 -> 65
            9 -> 66
            17 -> 66
            8 -> 66
            69 -> 68
            70 -> 68
            52 -> 68
            45 -> 69
            47 -> 69
            48 -> 69
            50 -> 70
            71 -> 70
            48 -> 71
            32 -> 72
            33 -> 72
            74 -> 73
            77 -> 73
            20 -> 74
            21 -> 74
            6 -> 74
            13 -> 74
            18 -> 74
            17 -> 74
            75 -> 74
            63 -> 74
            77 -> 74
            58 -> 74
            6 -> 75
            13 -> 75
            14 -> 75
            27 -> 75
            28 -> 75
            30 -> 75
            31 -> 75
            35 -> 75
            37 -> 75
            38 -> 75
            32 -> 75
            9 -> 75
            76 -> 75
            77 -> 75
            41 -> 75
            18 -> 75
            17 -> 75
            33 -> 75
            42 -> 75
            78 -> 75
            53 -> 75
            82 -> 75
            39 -> 75
            56 -> 75
            57 -> 75
            58 -> 75
            59 -> 75
            60 -> 75
            61 -> 75
            9 -> 76
            17 -> 76
            8 -> 76
            79 -> 78
            80 -> 78
            52 -> 78
            45 -> 79
            47 -> 79
            48 -> 79
            50 -> 80
            81 -> 80
            48 -> 81
            32 -> 82
            33 -> 82
            17 -> 83
            3 -> 84
            17 -> 84
            63 -> 85
            17 -> 85
            73 -> 86
            17 -> 86
            3 -> 87
            17 -> 87
            63 -> 88
            17 -> 88
            73 -> 89
            17 -> 89
    }

|


Scaling-Up
==========

If you now feel confident and want to tackle runs with larger temporal, technological and
spatial scope, clean-up the repository and after modifying the ``config/config.yaml`` file
target the collection rule ``all`` again without providing the test
configuration file.

.. code:: bash

    snakemake -call purge
    snakemake -call all

.. note::

    It is good practice to perform a dry-run using the option `-n`, before you
    commit to a run:

    .. code:: bash

        snakemake -call all -n
