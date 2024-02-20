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
    base_network                                            1
    build_ammonia_production                                1
    build_biomass_potentials                                1
    build_bus_regions                                       1
    build_clustered_population_layouts                      1
    build_cop_profiles                                      1
    build_daily_heat_demand                                 1
    build_district_heat_share                               1
    build_electricity_demand                                1
    build_energy_totals                                     1
    build_gas_input_locations                               1
    build_gas_network                                       1
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
    build_population_weighted_energy_totals                 1
    build_powerplants                                       1
    build_renewable_profiles                                4
    build_salt_cavern_potentials                            1
    build_shapes                                            1
    build_ship_raster                                       1
    build_shipping_demand                                   1
    build_simplified_population_layouts                     1
    build_temperature_profiles                              3
    build_transport_demand                                  1
    cluster_gas_network                                     1
    cluster_network                                         1
    copy_config                                             1
    make_summary                                            1
    plot_gas_network                                        1
    plot_hydrogen_network                                   1
    plot_power_network                                      1
    plot_power_network_clustered                            1
    plot_summary                                            1
    prepare_network                                         1
    prepare_sector_network                                  1
    retrieve_cost_data                                      1
    retrieve_databundle                                     1
    retrieve_electricity_demand                             1
    retrieve_gas_infrastructure_data                        1
    retrieve_natura_raster                                  1
    retrieve_sector_databundle                              1
    simplify_network                                        1
    solve_sector_network                                    1
    total                                                  60

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
            0[label = "all", color = "0.55 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.31 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.37 0.6 0.85", style="rounded"];
            3[label = "plot_power_network_clustered", color = "0.50 0.6 0.85", style="rounded"];
            4[label = "cluster_network\nclusters: 5", color = "0.62 0.6 0.85", style="rounded"];
            5[label = "simplify_network\nsimpl: ", color = "0.18 0.6 0.85", style="rounded"];
            6[label = "add_electricity", color = "0.33 0.6 0.85", style="rounded"];
            7[label = "build_renewable_profiles\ntechnology: solar", color = "0.20 0.6 0.85", style="rounded"];
            8[label = "base_network", color = "0.31 0.6 0.85", style="rounded"];
            9[label = "build_shapes", color = "0.36 0.6 0.85", style="rounded"];
            10[label = "retrieve_databundle", color = "0.29 0.6 0.85", style="rounded"];
            11[label = "retrieve_natura_raster", color = "0.01 0.6 0.85", style="rounded"];
            12[label = "build_bus_regions", color = "0.10 0.6 0.85", style="rounded"];
            13[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.37 0.6 0.85", style="rounded,dashed"];
            14[label = "build_renewable_profiles\ntechnology: onwind", color = "0.20 0.6 0.85", style="rounded"];
            15[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.20 0.6 0.85", style="rounded"];
            16[label = "build_ship_raster", color = "0.64 0.6 0.85", style="rounded"];
            17[label = "retrieve_ship_raster", color = "0.64 0.6 0.85", style="rounded,dashed"];
            18[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.20 0.6 0.85", style="rounded"];
            19[label = "retrieve_cost_data\nyear: 2030", color = "0.12 0.6 0.85", style="rounded"];
            20[label = "build_powerplants", color = "0.23 0.6 0.85", style="rounded"];
            21[label = "build_electricity_demand", color = "0.54 0.6 0.85", style="rounded"];
            22[label = "retrieve_electricity_demand", color = "0.07 0.6 0.85", style="rounded"];
            23[label = "solve_sector_network", color = "0.41 0.6 0.85", style="rounded"];
            24[label = "prepare_sector_network\nsector_opts: CO2L0-24h-T-H-B-I-A-dist1", color = "0.22 0.6 0.85", style="rounded"];
            25[label = "cluster_gas_network", color = "0.24 0.6 0.85", style="rounded"];
            26[label = "build_gas_network", color = "0.10 0.6 0.85", style="rounded"];
            27[label = "retrieve_gas_infrastructure_data", color = "0.17 0.6 0.85", style="rounded"];
            28[label = "build_gas_input_locations", color = "0.16 0.6 0.85", style="rounded"];
            29[label = "prepare_network\nll: v1.5\nopts: ", color = "0.49 0.6 0.85", style="rounded"];
            30[label = "add_extra_components", color = "0.14 0.6 0.85", style="rounded"];
            31[label = "build_energy_totals", color = "0.39 0.6 0.85", style="rounded"];
            32[label = "retrieve_sector_databundle", color = "0.58 0.6 0.85", style="rounded"];
            33[label = "build_population_weighted_energy_totals", color = "0.56 0.6 0.85", style="rounded"];
            34[label = "build_clustered_population_layouts", color = "0.49 0.6 0.85", style="rounded"];
            35[label = "build_population_layouts", color = "0.06 0.6 0.85", style="rounded"];
            36[label = "build_shipping_demand", color = "0.47 0.6 0.85", style="rounded"];
            37[label = "build_transport_demand", color = "0.45 0.6 0.85", style="rounded"];
            38[label = "build_temperature_profiles\nscope: total", color = "0.04 0.6 0.85", style="rounded"];
            39[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.11 0.6 0.85", style="rounded"];
            40[label = "build_salt_cavern_potentials", color = "0.15 0.6 0.85", style="rounded"];
            41[label = "build_simplified_population_layouts", color = "0.46 0.6 0.85", style="rounded"];
            42[label = "build_industrial_energy_demand_per_node", color = "0.63 0.6 0.85", style="rounded"];
            43[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.07 0.6 0.85", style="rounded"];
            44[label = "build_industry_sector_ratios", color = "0.59 0.6 0.85", style="rounded"];
            45[label = "build_ammonia_production", color = "0.04 0.6 0.85", style="rounded"];
            46[label = "build_industrial_energy_demand_per_country_today", color = "0.44 0.6 0.85", style="rounded"];
            47[label = "build_industrial_production_per_country", color = "0.34 0.6 0.85", style="rounded"];
            48[label = "build_industrial_production_per_node", color = "0.26 0.6 0.85", style="rounded"];
            49[label = "build_industrial_distribution_key", color = "0.13 0.6 0.85", style="rounded"];
            50[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.32 0.6 0.85", style="rounded"];
            51[label = "build_industrial_energy_demand_per_node_today", color = "0.48 0.6 0.85", style="rounded"];
            52[label = "build_hourly_heat_demand", color = "0.28 0.6 0.85", style="rounded"];
            53[label = "build_daily_heat_demand\nscope: total", color = "0.28 0.6 0.85", style="rounded"];
            54[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.52 0.6 0.85", style="rounded"];
            55[label = "build_temperature_profiles\nscope: rural", color = "0.04 0.6 0.85", style="rounded"];
            56[label = "build_temperature_profiles\nscope: urban", color = "0.04 0.6 0.85", style="rounded"];
            57[label = "build_cop_profiles", color = "0.38 0.6 0.85", style="rounded"];
            58[label = "copy_config", color = "0.19 0.6 0.85", style="rounded"];
            59[label = "plot_power_network", color = "0.60 0.6 0.85", style="rounded"];
            60[label = "plot_hydrogen_network", color = "0.27 0.6 0.85", style="rounded"];
            61[label = "plot_gas_network", color = "0.08 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            32 -> 1
            3 -> 2
            23 -> 2
            19 -> 2
            59 -> 2
            60 -> 2
            61 -> 2
            4 -> 3
            5 -> 4
            19 -> 4
            6 -> 5
            19 -> 5
            12 -> 5
            7 -> 6
            14 -> 6
            15 -> 6
            18 -> 6
            8 -> 6
            19 -> 6
            12 -> 6
            20 -> 6
            10 -> 6
            21 -> 6
            9 -> 6
            8 -> 7
            10 -> 7
            11 -> 7
            9 -> 7
            12 -> 7
            13 -> 7
            9 -> 8
            10 -> 9
            9 -> 12
            8 -> 12
            8 -> 14
            10 -> 14
            11 -> 14
            9 -> 14
            12 -> 14
            13 -> 14
            8 -> 15
            10 -> 15
            11 -> 15
            16 -> 15
            9 -> 15
            12 -> 15
            13 -> 15
            17 -> 16
            13 -> 16
            8 -> 18
            10 -> 18
            11 -> 18
            16 -> 18
            9 -> 18
            12 -> 18
            13 -> 18
            8 -> 20
            22 -> 21
            24 -> 23
            58 -> 23
            25 -> 24
            28 -> 24
            29 -> 24
            31 -> 24
            32 -> 24
            33 -> 24
            36 -> 24
            37 -> 24
            39 -> 24
            19 -> 24
            15 -> 24
            18 -> 24
            40 -> 24
            5 -> 24
            4 -> 24
            34 -> 24
            41 -> 24
            42 -> 24
            52 -> 24
            54 -> 24
            38 -> 24
            55 -> 24
            56 -> 24
            57 -> 24
            26 -> 25
            4 -> 25
            27 -> 26
            27 -> 28
            4 -> 28
            30 -> 29
            19 -> 29
            4 -> 30
            19 -> 30
            9 -> 31
            32 -> 31
            31 -> 33
            34 -> 33
            35 -> 34
            4 -> 34
            13 -> 34
            9 -> 35
            13 -> 35
            9 -> 36
            4 -> 36
            31 -> 36
            34 -> 37
            33 -> 37
            31 -> 37
            32 -> 37
            38 -> 37
            35 -> 38
            4 -> 38
            13 -> 38
            32 -> 39
            4 -> 39
            10 -> 39
            9 -> 39
            32 -> 40
            4 -> 40
            35 -> 41
            5 -> 41
            13 -> 41
            43 -> 42
            48 -> 42
            51 -> 42
            44 -> 43
            46 -> 43
            47 -> 43
            45 -> 44
            32 -> 44
            32 -> 45
            32 -> 46
            47 -> 46
            45 -> 47
            32 -> 47
            49 -> 48
            50 -> 48
            4 -> 49
            34 -> 49
            32 -> 49
            47 -> 50
            49 -> 51
            46 -> 51
            53 -> 52
            35 -> 53
            4 -> 53
            13 -> 53
            31 -> 54
            34 -> 54
            35 -> 55
            4 -> 55
            13 -> 55
            35 -> 56
            4 -> 56
            13 -> 56
            38 -> 57
            55 -> 57
            56 -> 57
            23 -> 59
            4 -> 59
            23 -> 60
            4 -> 60
            23 -> 61
            4 -> 61
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
            0[label = "all", color = "0.46 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.40 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.59 0.6 0.85", style="rounded"];
            3[label = "plot_power_network_clustered", color = "0.17 0.6 0.85", style="rounded"];
            4[label = "cluster_network\nclusters: 5", color = "0.49 0.6 0.85", style="rounded"];
            5[label = "simplify_network\nsimpl: ", color = "0.16 0.6 0.85", style="rounded"];
            6[label = "add_electricity", color = "0.32 0.6 0.85", style="rounded"];
            7[label = "build_renewable_profiles\ntechnology: solar", color = "0.63 0.6 0.85", style="rounded"];
            8[label = "base_network", color = "0.12 0.6 0.85", style="rounded"];
            9[label = "build_shapes", color = "0.23 0.6 0.85", style="rounded"];
            10[label = "retrieve_databundle", color = "0.61 0.6 0.85", style="rounded"];
            11[label = "retrieve_natura_raster", color = "0.50 0.6 0.85", style="rounded"];
            12[label = "build_bus_regions", color = "0.51 0.6 0.85", style="rounded"];
            13[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.37 0.6 0.85", style="rounded,dashed"];
            14[label = "build_renewable_profiles\ntechnology: onwind", color = "0.63 0.6 0.85", style="rounded"];
            15[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.63 0.6 0.85", style="rounded"];
            16[label = "build_ship_raster", color = "0.24 0.6 0.85", style="rounded"];
            17[label = "retrieve_ship_raster", color = "0.14 0.6 0.85", style="rounded,dashed"];
            18[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.63 0.6 0.85", style="rounded"];
            19[label = "retrieve_cost_data\nyear: 2030", color = "0.04 0.6 0.85", style="rounded"];
            20[label = "build_powerplants", color = "0.58 0.6 0.85", style="rounded"];
            21[label = "build_electricity_demand", color = "0.04 0.6 0.85", style="rounded"];
            22[label = "retrieve_electricity_demand", color = "0.62 0.6 0.85", style="rounded"];
            23[label = "solve_sector_network_myopic", color = "0.30 0.6 0.85", style="rounded"];
            24[label = "add_existing_baseyear", color = "0.34 0.6 0.85", style="rounded"];
            25[label = "prepare_sector_network\nsector_opts: 24h-T-H-B-I-A-dist1", color = "0.42 0.6 0.85", style="rounded"];
            26[label = "cluster_gas_network", color = "0.39 0.6 0.85", style="rounded"];
            27[label = "build_gas_network", color = "0.59 0.6 0.85", style="rounded"];
            28[label = "retrieve_gas_infrastructure_data", color = "0.15 0.6 0.85", style="rounded"];
            29[label = "build_gas_input_locations", color = "0.07 0.6 0.85", style="rounded"];
            30[label = "prepare_network\nll: v1.5\nopts: ", color = "0.56 0.6 0.85", style="rounded"];
            31[label = "add_extra_components", color = "0.11 0.6 0.85", style="rounded"];
            32[label = "build_energy_totals", color = "0.18 0.6 0.85", style="rounded"];
            33[label = "retrieve_sector_databundle", color = "0.06 0.6 0.85", style="rounded"];
            34[label = "build_population_weighted_energy_totals", color = "0.03 0.6 0.85", style="rounded"];
            35[label = "build_clustered_population_layouts", color = "0.25 0.6 0.85", style="rounded"];
            36[label = "build_population_layouts", color = "0.57 0.6 0.85", style="rounded"];
            37[label = "build_shipping_demand", color = "0.45 0.6 0.85", style="rounded"];
            38[label = "build_transport_demand", color = "0.18 0.6 0.85", style="rounded"];
            39[label = "build_temperature_profiles\nscope: total", color = "0.54 0.6 0.85", style="rounded"];
            40[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.41 0.6 0.85", style="rounded"];
            41[label = "build_salt_cavern_potentials", color = "0.02 0.6 0.85", style="rounded"];
            42[label = "build_simplified_population_layouts", color = "0.15 0.6 0.85", style="rounded"];
            43[label = "build_industrial_energy_demand_per_node", color = "0.47 0.6 0.85", style="rounded"];
            44[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.31 0.6 0.85", style="rounded"];
            45[label = "build_industry_sector_ratios", color = "0.48 0.6 0.85", style="rounded"];
            46[label = "build_ammonia_production", color = "0.00 0.6 0.85", style="rounded"];
            47[label = "build_industrial_energy_demand_per_country_today", color = "0.32 0.6 0.85", style="rounded"];
            48[label = "build_industrial_production_per_country", color = "0.60 0.6 0.85", style="rounded"];
            49[label = "build_industrial_production_per_node", color = "0.05 0.6 0.85", style="rounded"];
            50[label = "build_industrial_distribution_key", color = "0.21 0.6 0.85", style="rounded"];
            51[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.33 0.6 0.85", style="rounded"];
            52[label = "build_industrial_energy_demand_per_node_today", color = "0.62 0.6 0.85", style="rounded"];
            53[label = "build_hourly_heat_demand", color = "0.28 0.6 0.85", style="rounded"];
            54[label = "build_daily_heat_demand\nscope: total", color = "0.22 0.6 0.85", style="rounded"];
            55[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.21 0.6 0.85", style="rounded"];
            56[label = "build_temperature_profiles\nscope: rural", color = "0.54 0.6 0.85", style="rounded"];
            57[label = "build_temperature_profiles\nscope: urban", color = "0.54 0.6 0.85", style="rounded"];
            58[label = "build_cop_profiles", color = "0.52 0.6 0.85", style="rounded"];
            59[label = "build_existing_heating_distribution", color = "0.09 0.6 0.85", style="rounded"];
            60[label = "copy_config", color = "0.42 0.6 0.85", style="rounded"];
            61[label = "solve_sector_network_myopic", color = "0.30 0.6 0.85", style="rounded"];
            62[label = "add_brownfield", color = "0.10 0.6 0.85", style="rounded"];
            63[label = "prepare_sector_network\nsector_opts: 24h-T-H-B-I-A-dist1", color = "0.42 0.6 0.85", style="rounded"];
            64[label = "build_biomass_potentials\nplanning_horizons: 2040", color = "0.41 0.6 0.85", style="rounded"];
            65[label = "retrieve_cost_data\nyear: 2040", color = "0.04 0.6 0.85", style="rounded"];
            66[label = "build_industrial_energy_demand_per_node", color = "0.47 0.6 0.85", style="rounded"];
            67[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2040", color = "0.31 0.6 0.85", style="rounded"];
            68[label = "build_industrial_production_per_node", color = "0.05 0.6 0.85", style="rounded"];
            69[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2040", color = "0.33 0.6 0.85", style="rounded"];
            70[label = "build_district_heat_share\nplanning_horizons: 2040", color = "0.21 0.6 0.85", style="rounded"];
            71[label = "solve_sector_network_myopic", color = "0.30 0.6 0.85", style="rounded"];
            72[label = "add_brownfield", color = "0.10 0.6 0.85", style="rounded"];
            73[label = "prepare_sector_network\nsector_opts: 24h-T-H-B-I-A-dist1", color = "0.42 0.6 0.85", style="rounded"];
            74[label = "build_biomass_potentials\nplanning_horizons: 2050", color = "0.41 0.6 0.85", style="rounded"];
            75[label = "retrieve_cost_data\nyear: 2050", color = "0.04 0.6 0.85", style="rounded"];
            76[label = "build_industrial_energy_demand_per_node", color = "0.47 0.6 0.85", style="rounded"];
            77[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2050", color = "0.31 0.6 0.85", style="rounded"];
            78[label = "build_industrial_production_per_node", color = "0.05 0.6 0.85", style="rounded"];
            79[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2050", color = "0.33 0.6 0.85", style="rounded"];
            80[label = "build_district_heat_share\nplanning_horizons: 2050", color = "0.21 0.6 0.85", style="rounded"];
            81[label = "plot_power_network", color = "0.48 0.6 0.85", style="rounded"];
            82[label = "plot_power_network", color = "0.48 0.6 0.85", style="rounded"];
            83[label = "plot_power_network", color = "0.48 0.6 0.85", style="rounded"];
            84[label = "plot_hydrogen_network", color = "0.37 0.6 0.85", style="rounded"];
            85[label = "plot_hydrogen_network", color = "0.37 0.6 0.85", style="rounded"];
            86[label = "plot_hydrogen_network", color = "0.37 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            33 -> 1
            3 -> 2
            23 -> 2
            61 -> 2
            71 -> 2
            19 -> 2
            81 -> 2
            82 -> 2
            83 -> 2
            84 -> 2
            85 -> 2
            86 -> 2
            4 -> 3
            5 -> 4
            19 -> 4
            6 -> 5
            19 -> 5
            12 -> 5
            7 -> 6
            14 -> 6
            15 -> 6
            18 -> 6
            8 -> 6
            19 -> 6
            12 -> 6
            20 -> 6
            10 -> 6
            21 -> 6
            9 -> 6
            8 -> 7
            10 -> 7
            11 -> 7
            9 -> 7
            12 -> 7
            13 -> 7
            9 -> 8
            10 -> 9
            9 -> 12
            8 -> 12
            8 -> 14
            10 -> 14
            11 -> 14
            9 -> 14
            12 -> 14
            13 -> 14
            8 -> 15
            10 -> 15
            11 -> 15
            16 -> 15
            9 -> 15
            12 -> 15
            13 -> 15
            17 -> 16
            13 -> 16
            8 -> 18
            10 -> 18
            11 -> 18
            16 -> 18
            9 -> 18
            12 -> 18
            13 -> 18
            8 -> 20
            22 -> 21
            24 -> 23
            19 -> 23
            60 -> 23
            25 -> 24
            20 -> 24
            5 -> 24
            4 -> 24
            35 -> 24
            19 -> 24
            58 -> 24
            59 -> 24
            26 -> 25
            29 -> 25
            30 -> 25
            32 -> 25
            33 -> 25
            34 -> 25
            37 -> 25
            38 -> 25
            40 -> 25
            19 -> 25
            15 -> 25
            18 -> 25
            41 -> 25
            5 -> 25
            4 -> 25
            35 -> 25
            42 -> 25
            43 -> 25
            53 -> 25
            55 -> 25
            39 -> 25
            56 -> 25
            57 -> 25
            58 -> 25
            27 -> 26
            4 -> 26
            28 -> 27
            28 -> 29
            4 -> 29
            31 -> 30
            19 -> 30
            4 -> 31
            19 -> 31
            9 -> 32
            33 -> 32
            32 -> 34
            35 -> 34
            36 -> 35
            4 -> 35
            13 -> 35
            9 -> 36
            13 -> 36
            9 -> 37
            4 -> 37
            32 -> 37
            35 -> 38
            34 -> 38
            32 -> 38
            33 -> 38
            39 -> 38
            36 -> 39
            4 -> 39
            13 -> 39
            33 -> 40
            4 -> 40
            10 -> 40
            9 -> 40
            33 -> 41
            4 -> 41
            36 -> 42
            5 -> 42
            13 -> 42
            44 -> 43
            49 -> 43
            52 -> 43
            45 -> 44
            47 -> 44
            48 -> 44
            46 -> 45
            33 -> 45
            33 -> 46
            33 -> 47
            48 -> 47
            46 -> 48
            33 -> 48
            50 -> 49
            51 -> 49
            4 -> 50
            35 -> 50
            33 -> 50
            48 -> 51
            50 -> 52
            47 -> 52
            54 -> 53
            36 -> 54
            4 -> 54
            13 -> 54
            32 -> 55
            35 -> 55
            36 -> 56
            4 -> 56
            13 -> 56
            36 -> 57
            4 -> 57
            13 -> 57
            39 -> 58
            56 -> 58
            57 -> 58
            35 -> 59
            34 -> 59
            55 -> 59
            62 -> 61
            65 -> 61
            60 -> 61
            7 -> 62
            14 -> 62
            15 -> 62
            18 -> 62
            5 -> 62
            4 -> 62
            63 -> 62
            23 -> 62
            65 -> 62
            58 -> 62
            26 -> 63
            29 -> 63
            30 -> 63
            32 -> 63
            33 -> 63
            34 -> 63
            37 -> 63
            38 -> 63
            64 -> 63
            65 -> 63
            15 -> 63
            18 -> 63
            41 -> 63
            5 -> 63
            4 -> 63
            35 -> 63
            42 -> 63
            66 -> 63
            53 -> 63
            70 -> 63
            39 -> 63
            56 -> 63
            57 -> 63
            58 -> 63
            33 -> 64
            4 -> 64
            10 -> 64
            9 -> 64
            67 -> 66
            68 -> 66
            52 -> 66
            45 -> 67
            47 -> 67
            48 -> 67
            50 -> 68
            69 -> 68
            48 -> 69
            32 -> 70
            35 -> 70
            72 -> 71
            75 -> 71
            60 -> 71
            7 -> 72
            14 -> 72
            15 -> 72
            18 -> 72
            5 -> 72
            4 -> 72
            73 -> 72
            61 -> 72
            75 -> 72
            58 -> 72
            26 -> 73
            29 -> 73
            30 -> 73
            32 -> 73
            33 -> 73
            34 -> 73
            37 -> 73
            38 -> 73
            74 -> 73
            75 -> 73
            15 -> 73
            18 -> 73
            41 -> 73
            5 -> 73
            4 -> 73
            35 -> 73
            42 -> 73
            76 -> 73
            53 -> 73
            80 -> 73
            39 -> 73
            56 -> 73
            57 -> 73
            58 -> 73
            33 -> 74
            4 -> 74
            10 -> 74
            9 -> 74
            77 -> 76
            78 -> 76
            52 -> 76
            45 -> 77
            47 -> 77
            48 -> 77
            50 -> 78
            79 -> 78
            48 -> 79
            32 -> 80
            35 -> 80
            23 -> 81
            4 -> 81
            61 -> 82
            4 -> 82
            71 -> 83
            4 -> 83
            23 -> 84
            4 -> 84
            61 -> 85
            4 -> 85
            71 -> 86
            4 -> 86
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
