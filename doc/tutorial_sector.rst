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
    retrieve_natura_raster                                  1
    retrieve_sector_databundle                              1
    retrieve_ship_raster                                    1
    retrieve_synthetic_electricity_demand                   1
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
            0[label = "all", color = "0.66 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.20 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.02 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network", color = "0.11 0.6 0.85", style="rounded"];
            4[label = "prepare_sector_network\nsector_opts: CO2L0-24h-T-H-B-I-A-dist1", color = "0.22 0.6 0.85", style="rounded"];
            5[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.26 0.6 0.85", style="rounded"];
            6[label = "base_network", color = "0.53 0.6 0.85", style="rounded"];
            7[label = "build_shapes", color = "0.04 0.6 0.85", style="rounded"];
            8[label = "retrieve_databundle", color = "0.49 0.6 0.85", style="rounded"];
            9[label = "retrieve_natura_raster", color = "0.46 0.6 0.85", style="rounded"];
            10[label = "build_ship_raster", color = "0.29 0.6 0.85", style="rounded"];
            11[label = "retrieve_ship_raster", color = "0.42 0.6 0.85", style="rounded"];
            12[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.27 0.6 0.85", style="rounded"];
            13[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.26 0.6 0.85", style="rounded"];
            14[label = "cluster_gas_network", color = "0.48 0.6 0.85", style="rounded"];
            15[label = "build_gas_network", color = "0.27 0.6 0.85", style="rounded"];
            16[label = "retrieve_gas_infrastructure_data", color = "0.38 0.6 0.85", style="rounded"];
            17[label = "cluster_network\nclusters: 5", color = "0.58 0.6 0.85", style="rounded"];
            18[label = "simplify_network\nsimpl: ", color = "0.55 0.6 0.85", style="rounded"];
            19[label = "add_electricity", color = "0.37 0.6 0.85", style="rounded"];
            20[label = "build_renewable_profiles\ntechnology: solar", color = "0.26 0.6 0.85", style="rounded"];
            21[label = "build_renewable_profiles\ntechnology: onwind", color = "0.26 0.6 0.85", style="rounded"];
            22[label = "retrieve_cost_data\nyear: 2030", color = "0.14 0.6 0.85", style="rounded"];
            23[label = "build_powerplants", color = "0.64 0.6 0.85", style="rounded"];
            24[label = "build_electricity_demand", color = "0.61 0.6 0.85", style="rounded"];
            25[label = "retrieve_electricity_demand", color = "0.08 0.6 0.85", style="rounded"];
            26[label = "retrieve_synthetic_electricity_demand", color = "0.36 0.6 0.85", style="rounded"];
            27[label = "build_gas_input_locations", color = "0.44 0.6 0.85", style="rounded"];
            28[label = "prepare_network\nll: v1.5\nopts: ", color = "0.25 0.6 0.85", style="rounded"];
            29[label = "add_extra_components", color = "0.39 0.6 0.85", style="rounded"];
            30[label = "retrieve_eurostat_data", color = "0.20 0.6 0.85", style="rounded"];
            31[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.58 0.6 0.85", style="rounded"];
            32[label = "build_energy_totals", color = "0.44 0.6 0.85", style="rounded"];
            33[label = "retrieve_sector_databundle", color = "0.60 0.6 0.85", style="rounded"];
            34[label = "build_clustered_population_layouts", color = "0.46 0.6 0.85", style="rounded"];
            35[label = "build_population_layouts", color = "0.43 0.6 0.85", style="rounded"];
            36[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.58 0.6 0.85", style="rounded"];
            37[label = "build_heat_totals", color = "0.11 0.6 0.85", style="rounded"];
            38[label = "build_shipping_demand", color = "0.16 0.6 0.85", style="rounded"];
            39[label = "build_transport_demand", color = "0.04 0.6 0.85", style="rounded"];
            40[label = "build_temperature_profiles\nscope: total", color = "0.28 0.6 0.85", style="rounded"];
            41[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.07 0.6 0.85", style="rounded"];
            42[label = "build_salt_cavern_potentials", color = "0.47 0.6 0.85", style="rounded"];
            43[label = "build_simplified_population_layouts", color = "0.29 0.6 0.85", style="rounded"];
            44[label = "build_industrial_energy_demand_per_node", color = "0.39 0.6 0.85", style="rounded"];
            45[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.57 0.6 0.85", style="rounded"];
            46[label = "build_industry_sector_ratios", color = "0.55 0.6 0.85", style="rounded"];
            47[label = "build_ammonia_production", color = "0.00 0.6 0.85", style="rounded"];
            48[label = "build_industrial_energy_demand_per_country_today", color = "0.52 0.6 0.85", style="rounded"];
            49[label = "build_industrial_production_per_country", color = "0.19 0.6 0.85", style="rounded"];
            50[label = "build_industrial_production_per_node", color = "0.21 0.6 0.85", style="rounded"];
            51[label = "build_industrial_distribution_key", color = "0.10 0.6 0.85", style="rounded"];
            52[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.63 0.6 0.85", style="rounded"];
            53[label = "build_industrial_energy_demand_per_node_today", color = "0.18 0.6 0.85", style="rounded"];
            54[label = "build_hourly_heat_demand", color = "0.13 0.6 0.85", style="rounded"];
            55[label = "build_daily_heat_demand\nscope: total", color = "0.22 0.6 0.85", style="rounded"];
            56[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.34 0.6 0.85", style="rounded"];
            57[label = "build_temperature_profiles\nscope: rural", color = "0.28 0.6 0.85", style="rounded"];
            58[label = "build_temperature_profiles\nscope: urban", color = "0.28 0.6 0.85", style="rounded"];
            59[label = "build_cop_profiles", color = "0.65 0.6 0.85", style="rounded"];
            60[label = "build_solar_thermal_profiles\nscope: total", color = "0.54 0.6 0.85", style="rounded"];
            61[label = "build_solar_thermal_profiles\nscope: urban", color = "0.54 0.6 0.85", style="rounded"];
            62[label = "build_solar_thermal_profiles\nscope: rural", color = "0.54 0.6 0.85", style="rounded"];
            63[label = "plot_power_network_clustered", color = "0.15 0.6 0.85", style="rounded"];
            64[label = "plot_power_network", color = "0.56 0.6 0.85", style="rounded"];
            65[label = "plot_hydrogen_network", color = "0.60 0.6 0.85", style="rounded"];
            66[label = "plot_gas_network", color = "0.53 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            30 -> 1
            33 -> 1
            3 -> 2
            22 -> 2
            63 -> 2
            64 -> 2
            65 -> 2
            66 -> 2
            4 -> 3
            5 -> 4
            13 -> 4
            14 -> 4
            27 -> 4
            28 -> 4
            30 -> 4
            31 -> 4
            36 -> 4
            38 -> 4
            39 -> 4
            32 -> 4
            33 -> 4
            41 -> 4
            22 -> 4
            42 -> 4
            18 -> 4
            17 -> 4
            34 -> 4
            43 -> 4
            44 -> 4
            54 -> 4
            56 -> 4
            40 -> 4
            57 -> 4
            58 -> 4
            59 -> 4
            60 -> 4
            61 -> 4
            62 -> 4
            6 -> 5
            8 -> 5
            9 -> 5
            10 -> 5
            7 -> 5
            12 -> 5
            7 -> 6
            8 -> 7
            11 -> 10
            12 -> 10
            6 -> 13
            8 -> 13
            9 -> 13
            10 -> 13
            7 -> 13
            12 -> 13
            15 -> 14
            17 -> 14
            16 -> 15
            18 -> 17
            22 -> 17
            19 -> 18
            22 -> 18
            6 -> 18
            20 -> 19
            21 -> 19
            5 -> 19
            13 -> 19
            6 -> 19
            22 -> 19
            23 -> 19
            8 -> 19
            24 -> 19
            7 -> 19
            6 -> 20
            8 -> 20
            9 -> 20
            7 -> 20
            12 -> 20
            6 -> 21
            8 -> 21
            9 -> 21
            7 -> 21
            12 -> 21
            6 -> 23
            25 -> 24
            26 -> 24
            16 -> 27
            17 -> 27
            29 -> 28
            22 -> 28
            17 -> 29
            22 -> 29
            32 -> 31
            34 -> 31
            7 -> 32
            33 -> 32
            30 -> 32
            35 -> 34
            17 -> 34
            12 -> 34
            7 -> 35
            12 -> 35
            37 -> 36
            34 -> 36
            32 -> 37
            7 -> 38
            17 -> 38
            32 -> 38
            34 -> 39
            31 -> 39
            32 -> 39
            33 -> 39
            40 -> 39
            35 -> 40
            17 -> 40
            12 -> 40
            33 -> 41
            17 -> 41
            8 -> 41
            7 -> 41
            33 -> 42
            17 -> 42
            35 -> 43
            18 -> 43
            12 -> 43
            45 -> 44
            50 -> 44
            53 -> 44
            46 -> 45
            48 -> 45
            49 -> 45
            47 -> 46
            33 -> 46
            33 -> 47
            33 -> 48
            49 -> 48
            47 -> 49
            33 -> 49
            30 -> 49
            51 -> 50
            52 -> 50
            17 -> 51
            34 -> 51
            33 -> 51
            49 -> 52
            51 -> 53
            48 -> 53
            55 -> 54
            35 -> 55
            17 -> 55
            12 -> 55
            32 -> 56
            34 -> 56
            35 -> 57
            17 -> 57
            12 -> 57
            35 -> 58
            17 -> 58
            12 -> 58
            40 -> 59
            57 -> 59
            58 -> 59
            35 -> 60
            17 -> 60
            12 -> 60
            35 -> 61
            17 -> 61
            12 -> 61
            35 -> 62
            17 -> 62
            12 -> 62
            17 -> 63
            3 -> 64
            17 -> 64
            3 -> 65
            17 -> 65
            3 -> 66
            17 -> 66
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
            0[label = "all", color = "0.20 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.55 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.21 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network_myopic", color = "0.50 0.6 0.85", style="rounded"];
            4[label = "add_existing_baseyear", color = "0.38 0.6 0.85", style="rounded"];
            5[label = "prepare_sector_network\nsector_opts: 24h-T-H-B-I-A-dist1", color = "0.53 0.6 0.85", style="rounded"];
            6[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.54 0.6 0.85", style="rounded"];
            7[label = "base_network", color = "0.12 0.6 0.85", style="rounded"];
            8[label = "build_shapes", color = "0.21 0.6 0.85", style="rounded"];
            9[label = "retrieve_databundle", color = "0.41 0.6 0.85", style="rounded"];
            10[label = "retrieve_natura_raster", color = "0.38 0.6 0.85", style="rounded"];
            11[label = "build_ship_raster", color = "0.05 0.6 0.85", style="rounded"];
            12[label = "retrieve_ship_raster", color = "0.25 0.6 0.85", style="rounded"];
            13[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.23 0.6 0.85", style="rounded"];
            14[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.54 0.6 0.85", style="rounded"];
            15[label = "cluster_gas_network", color = "0.16 0.6 0.85", style="rounded"];
            16[label = "build_gas_network", color = "0.26 0.6 0.85", style="rounded"];
            17[label = "retrieve_gas_infrastructure_data", color = "0.04 0.6 0.85", style="rounded"];
            18[label = "cluster_network\nclusters: 5", color = "0.10 0.6 0.85", style="rounded"];
            19[label = "simplify_network\nsimpl: ", color = "0.02 0.6 0.85", style="rounded"];
            20[label = "add_electricity", color = "0.56 0.6 0.85", style="rounded"];
            21[label = "build_renewable_profiles\ntechnology: solar", color = "0.54 0.6 0.85", style="rounded"];
            22[label = "build_renewable_profiles\ntechnology: onwind", color = "0.54 0.6 0.85", style="rounded"];
            23[label = "retrieve_cost_data\nyear: 2030", color = "0.16 0.6 0.85", style="rounded"];
            24[label = "build_powerplants", color = "0.63 0.6 0.85", style="rounded"];
            25[label = "build_electricity_demand", color = "0.57 0.6 0.85", style="rounded"];
            26[label = "retrieve_electricity_demand", color = "0.27 0.6 0.85", style="rounded"];
            27[label = "retrieve_synthetic_electricity_demand", color = "0.58 0.6 0.85", style="rounded"];
            28[label = "build_gas_input_locations", color = "0.28 0.6 0.85", style="rounded"];
            29[label = "prepare_network\nll: v1.5\nopts: ", color = "0.14 0.6 0.85", style="rounded"];
            30[label = "add_extra_components", color = "0.14 0.6 0.85", style="rounded"];
            31[label = "retrieve_eurostat_data", color = "0.58 0.6 0.85", style="rounded"];
            32[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.36 0.6 0.85", style="rounded"];
            33[label = "build_energy_totals", color = "0.65 0.6 0.85", style="rounded"];
            34[label = "retrieve_sector_databundle", color = "0.46 0.6 0.85", style="rounded"];
            35[label = "build_clustered_population_layouts", color = "0.52 0.6 0.85", style="rounded"];
            36[label = "build_population_layouts", color = "0.13 0.6 0.85", style="rounded"];
            37[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.36 0.6 0.85", style="rounded"];
            38[label = "build_heat_totals", color = "0.31 0.6 0.85", style="rounded"];
            39[label = "build_shipping_demand", color = "0.01 0.6 0.85", style="rounded"];
            40[label = "build_transport_demand", color = "0.51 0.6 0.85", style="rounded"];
            41[label = "build_temperature_profiles\nscope: total", color = "0.00 0.6 0.85", style="rounded"];
            42[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.18 0.6 0.85", style="rounded"];
            43[label = "build_salt_cavern_potentials", color = "0.25 0.6 0.85", style="rounded"];
            44[label = "build_simplified_population_layouts", color = "0.27 0.6 0.85", style="rounded"];
            45[label = "build_industrial_energy_demand_per_node", color = "0.30 0.6 0.85", style="rounded"];
            46[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.41 0.6 0.85", style="rounded"];
            47[label = "build_industry_sector_ratios", color = "0.03 0.6 0.85", style="rounded"];
            48[label = "build_ammonia_production", color = "0.37 0.6 0.85", style="rounded"];
            49[label = "build_industrial_energy_demand_per_country_today", color = "0.10 0.6 0.85", style="rounded"];
            50[label = "build_industrial_production_per_country", color = "0.03 0.6 0.85", style="rounded"];
            51[label = "build_industrial_production_per_node", color = "0.63 0.6 0.85", style="rounded"];
            52[label = "build_industrial_distribution_key", color = "0.17 0.6 0.85", style="rounded"];
            53[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.06 0.6 0.85", style="rounded"];
            54[label = "build_industrial_energy_demand_per_node_today", color = "0.08 0.6 0.85", style="rounded"];
            55[label = "build_hourly_heat_demand", color = "0.08 0.6 0.85", style="rounded"];
            56[label = "build_daily_heat_demand\nscope: total", color = "0.60 0.6 0.85", style="rounded"];
            57[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.32 0.6 0.85", style="rounded"];
            58[label = "build_temperature_profiles\nscope: rural", color = "0.00 0.6 0.85", style="rounded"];
            59[label = "build_temperature_profiles\nscope: urban", color = "0.00 0.6 0.85", style="rounded"];
            60[label = "build_cop_profiles", color = "0.11 0.6 0.85", style="rounded"];
            61[label = "build_solar_thermal_profiles\nscope: total", color = "0.01 0.6 0.85", style="rounded"];
            62[label = "build_solar_thermal_profiles\nscope: urban", color = "0.01 0.6 0.85", style="rounded"];
            63[label = "build_solar_thermal_profiles\nscope: rural", color = "0.01 0.6 0.85", style="rounded"];
            64[label = "build_existing_heating_distribution", color = "0.40 0.6 0.85", style="rounded"];
            65[label = "solve_sector_network_myopic", color = "0.50 0.6 0.85", style="rounded"];
            66[label = "add_brownfield", color = "0.45 0.6 0.85", style="rounded"];
            67[label = "prepare_sector_network\nsector_opts: 24h-T-H-B-I-A-dist1", color = "0.53 0.6 0.85", style="rounded"];
            68[label = "build_biomass_potentials\nplanning_horizons: 2040", color = "0.18 0.6 0.85", style="rounded"];
            69[label = "retrieve_cost_data\nyear: 2040", color = "0.16 0.6 0.85", style="rounded"];
            70[label = "build_industrial_energy_demand_per_node", color = "0.30 0.6 0.85", style="rounded"];
            71[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2040", color = "0.41 0.6 0.85", style="rounded"];
            72[label = "build_industrial_production_per_node", color = "0.63 0.6 0.85", style="rounded"];
            73[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2040", color = "0.06 0.6 0.85", style="rounded"];
            74[label = "build_district_heat_share\nplanning_horizons: 2040", color = "0.32 0.6 0.85", style="rounded"];
            75[label = "solve_sector_network_myopic", color = "0.50 0.6 0.85", style="rounded"];
            76[label = "add_brownfield", color = "0.45 0.6 0.85", style="rounded"];
            77[label = "prepare_sector_network\nsector_opts: 24h-T-H-B-I-A-dist1", color = "0.53 0.6 0.85", style="rounded"];
            78[label = "build_biomass_potentials\nplanning_horizons: 2050", color = "0.18 0.6 0.85", style="rounded"];
            79[label = "retrieve_cost_data\nyear: 2050", color = "0.16 0.6 0.85", style="rounded"];
            80[label = "build_industrial_energy_demand_per_node", color = "0.30 0.6 0.85", style="rounded"];
            81[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2050", color = "0.41 0.6 0.85", style="rounded"];
            82[label = "build_industrial_production_per_node", color = "0.63 0.6 0.85", style="rounded"];
            83[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2050", color = "0.06 0.6 0.85", style="rounded"];
            84[label = "build_district_heat_share\nplanning_horizons: 2050", color = "0.32 0.6 0.85", style="rounded"];
            85[label = "plot_power_network_clustered", color = "0.09 0.6 0.85", style="rounded"];
            86[label = "plot_power_network", color = "0.43 0.6 0.85", style="rounded"];
            87[label = "plot_power_network", color = "0.43 0.6 0.85", style="rounded"];
            88[label = "plot_power_network", color = "0.43 0.6 0.85", style="rounded"];
            89[label = "plot_hydrogen_network", color = "0.33 0.6 0.85", style="rounded"];
            90[label = "plot_hydrogen_network", color = "0.33 0.6 0.85", style="rounded"];
            91[label = "plot_hydrogen_network", color = "0.33 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            31 -> 1
            34 -> 1
            3 -> 2
            65 -> 2
            75 -> 2
            23 -> 2
            85 -> 2
            86 -> 2
            87 -> 2
            88 -> 2
            89 -> 2
            90 -> 2
            91 -> 2
            4 -> 3
            23 -> 3
            5 -> 4
            24 -> 4
            19 -> 4
            18 -> 4
            35 -> 4
            23 -> 4
            60 -> 4
            64 -> 4
            6 -> 5
            14 -> 5
            15 -> 5
            28 -> 5
            29 -> 5
            31 -> 5
            32 -> 5
            37 -> 5
            39 -> 5
            40 -> 5
            33 -> 5
            34 -> 5
            42 -> 5
            23 -> 5
            43 -> 5
            19 -> 5
            18 -> 5
            35 -> 5
            44 -> 5
            45 -> 5
            55 -> 5
            57 -> 5
            41 -> 5
            58 -> 5
            59 -> 5
            60 -> 5
            61 -> 5
            62 -> 5
            63 -> 5
            7 -> 6
            9 -> 6
            10 -> 6
            11 -> 6
            8 -> 6
            13 -> 6
            8 -> 7
            9 -> 8
            12 -> 11
            13 -> 11
            7 -> 14
            9 -> 14
            10 -> 14
            11 -> 14
            8 -> 14
            13 -> 14
            16 -> 15
            18 -> 15
            17 -> 16
            19 -> 18
            23 -> 18
            20 -> 19
            23 -> 19
            7 -> 19
            21 -> 20
            22 -> 20
            6 -> 20
            14 -> 20
            7 -> 20
            23 -> 20
            24 -> 20
            9 -> 20
            25 -> 20
            8 -> 20
            7 -> 21
            9 -> 21
            10 -> 21
            8 -> 21
            13 -> 21
            7 -> 22
            9 -> 22
            10 -> 22
            8 -> 22
            13 -> 22
            7 -> 24
            26 -> 25
            27 -> 25
            17 -> 28
            18 -> 28
            30 -> 29
            23 -> 29
            18 -> 30
            23 -> 30
            33 -> 32
            35 -> 32
            8 -> 33
            34 -> 33
            31 -> 33
            36 -> 35
            18 -> 35
            13 -> 35
            8 -> 36
            13 -> 36
            38 -> 37
            35 -> 37
            33 -> 38
            8 -> 39
            18 -> 39
            33 -> 39
            35 -> 40
            32 -> 40
            33 -> 40
            34 -> 40
            41 -> 40
            36 -> 41
            18 -> 41
            13 -> 41
            34 -> 42
            18 -> 42
            9 -> 42
            8 -> 42
            34 -> 43
            18 -> 43
            36 -> 44
            19 -> 44
            13 -> 44
            46 -> 45
            51 -> 45
            54 -> 45
            47 -> 46
            49 -> 46
            50 -> 46
            48 -> 47
            34 -> 47
            34 -> 48
            34 -> 49
            50 -> 49
            48 -> 50
            34 -> 50
            31 -> 50
            52 -> 51
            53 -> 51
            18 -> 52
            35 -> 52
            34 -> 52
            50 -> 53
            52 -> 54
            49 -> 54
            56 -> 55
            36 -> 56
            18 -> 56
            13 -> 56
            33 -> 57
            35 -> 57
            36 -> 58
            18 -> 58
            13 -> 58
            36 -> 59
            18 -> 59
            13 -> 59
            41 -> 60
            58 -> 60
            59 -> 60
            36 -> 61
            18 -> 61
            13 -> 61
            36 -> 62
            18 -> 62
            13 -> 62
            36 -> 63
            18 -> 63
            13 -> 63
            35 -> 64
            32 -> 64
            57 -> 64
            66 -> 65
            69 -> 65
            21 -> 66
            22 -> 66
            6 -> 66
            14 -> 66
            19 -> 66
            18 -> 66
            67 -> 66
            3 -> 66
            69 -> 66
            60 -> 66
            6 -> 67
            14 -> 67
            15 -> 67
            28 -> 67
            29 -> 67
            31 -> 67
            32 -> 67
            37 -> 67
            39 -> 67
            40 -> 67
            33 -> 67
            34 -> 67
            68 -> 67
            69 -> 67
            43 -> 67
            19 -> 67
            18 -> 67
            35 -> 67
            44 -> 67
            70 -> 67
            55 -> 67
            74 -> 67
            41 -> 67
            58 -> 67
            59 -> 67
            60 -> 67
            61 -> 67
            62 -> 67
            63 -> 67
            34 -> 68
            18 -> 68
            9 -> 68
            8 -> 68
            71 -> 70
            72 -> 70
            54 -> 70
            47 -> 71
            49 -> 71
            50 -> 71
            52 -> 72
            73 -> 72
            50 -> 73
            33 -> 74
            35 -> 74
            76 -> 75
            79 -> 75
            21 -> 76
            22 -> 76
            6 -> 76
            14 -> 76
            19 -> 76
            18 -> 76
            77 -> 76
            65 -> 76
            79 -> 76
            60 -> 76
            6 -> 77
            14 -> 77
            15 -> 77
            28 -> 77
            29 -> 77
            31 -> 77
            32 -> 77
            37 -> 77
            39 -> 77
            40 -> 77
            33 -> 77
            34 -> 77
            78 -> 77
            79 -> 77
            43 -> 77
            19 -> 77
            18 -> 77
            35 -> 77
            44 -> 77
            80 -> 77
            55 -> 77
            84 -> 77
            41 -> 77
            58 -> 77
            59 -> 77
            60 -> 77
            61 -> 77
            62 -> 77
            63 -> 77
            34 -> 78
            18 -> 78
            9 -> 78
            8 -> 78
            81 -> 80
            82 -> 80
            54 -> 80
            47 -> 81
            49 -> 81
            50 -> 81
            52 -> 82
            83 -> 82
            50 -> 83
            33 -> 84
            35 -> 84
            18 -> 85
            3 -> 86
            18 -> 86
            65 -> 87
            18 -> 87
            75 -> 88
            18 -> 88
            3 -> 89
            18 -> 89
            65 -> 90
            18 -> 90
            75 -> 91
            18 -> 91
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
