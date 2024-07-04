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
    build_clustered_population_layouts                      1
    build_cop_profiles                                      1
    build_daily_heat_demand                                 1
    build_district_heat_share                               1
    build_electricity_demand                                1
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
    build_powerplants                                       1
    build_renewable_profiles                                6
    build_salt_cavern_potentials                            1
    build_shapes                                            1
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
    retrieve_eurostat_household_data                        1
    retrieve_gas_infrastructure_data                        1
    retrieve_ship_raster                                    1
    retrieve_synthetic_electricity_demand                   1
    simplify_network                                        1
    solve_sector_network                                    1
    time_aggregation                                        1
    total                                                  69

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
            0[label = "all", color = "0.28 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.60 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.30 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network", color = "0.36 0.6 0.85", style="rounded"];
            4[label = "prepare_sector_network\nsector_opts: ", color = "0.22 0.6 0.85", style="rounded"];
            5[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.20 0.6 0.85", style="rounded"];
            6[label = "base_network", color = "0.00 0.6 0.85", style="rounded"];
            7[label = "build_shapes", color = "0.25 0.6 0.85", style="rounded"];
            8[label = "retrieve_databundle", color = "0.06 0.6 0.85", style="rounded"];
            9[label = "build_ship_raster", color = "0.06 0.6 0.85", style="rounded"];
            10[label = "retrieve_ship_raster", color = "0.27 0.6 0.85", style="rounded"];
            11[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.26 0.6 0.85", style="rounded"];
            12[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.20 0.6 0.85", style="rounded"];
            13[label = "build_renewable_profiles\ntechnology: offwind-float", color = "0.20 0.6 0.85", style="rounded"];
            14[label = "cluster_gas_network", color = "0.37 0.6 0.85", style="rounded"];
            15[label = "build_gas_network", color = "0.44 0.6 0.85", style="rounded"];
            16[label = "retrieve_gas_infrastructure_data", color = "0.43 0.6 0.85", style="rounded"];
            17[label = "cluster_network\nclusters: 5", color = "0.08 0.6 0.85", style="rounded"];
            18[label = "simplify_network\nsimpl: ", color = "0.01 0.6 0.85", style="rounded"];
            19[label = "add_electricity", color = "0.53 0.6 0.85", style="rounded"];
            20[label = "build_renewable_profiles\ntechnology: solar", color = "0.20 0.6 0.85", style="rounded"];
            21[label = "build_renewable_profiles\ntechnology: solar-hsat", color = "0.20 0.6 0.85", style="rounded"];
            22[label = "build_renewable_profiles\ntechnology: onwind", color = "0.20 0.6 0.85", style="rounded"];
            23[label = "retrieve_cost_data\nyear: 2030", color = "0.11 0.6 0.85", style="rounded"];
            24[label = "build_powerplants", color = "0.62 0.6 0.85", style="rounded"];
            25[label = "build_electricity_demand", color = "0.66 0.6 0.85", style="rounded"];
            26[label = "retrieve_electricity_demand", color = "0.20 0.6 0.85", style="rounded"];
            27[label = "retrieve_synthetic_electricity_demand", color = "0.52 0.6 0.85", style="rounded"];
            28[label = "build_gas_input_locations", color = "0.21 0.6 0.85", style="rounded"];
            29[label = "time_aggregation", color = "0.58 0.6 0.85", style="rounded"];
            30[label = "prepare_network\nll: v1.5\nopts: ", color = "0.61 0.6 0.85", style="rounded"];
            31[label = "add_extra_components", color = "0.59 0.6 0.85", style="rounded"];
            32[label = "build_hourly_heat_demand", color = "0.48 0.6 0.85", style="rounded"];
            33[label = "build_daily_heat_demand\nscope: total", color = "0.12 0.6 0.85", style="rounded"];
            34[label = "build_population_layouts", color = "0.62 0.6 0.85", style="rounded"];
            35[label = "build_solar_thermal_profiles\nscope: total", color = "0.23 0.6 0.85", style="rounded"];
            36[label = "retrieve_eurostat_data", color = "0.45 0.6 0.85", style="rounded"];
            37[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.22 0.6 0.85", style="rounded"];
            38[label = "build_energy_totals", color = "0.65 0.6 0.85", style="rounded"];
            39[label = "retrieve_eurostat_household_data", color = "0.36 0.6 0.85", style="rounded"];
            40[label = "build_clustered_population_layouts", color = "0.02 0.6 0.85", style="rounded"];
            41[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.22 0.6 0.85", style="rounded"];
            42[label = "build_heat_totals", color = "0.53 0.6 0.85", style="rounded"];
            43[label = "build_shipping_demand", color = "0.17 0.6 0.85", style="rounded"];
            44[label = "build_transport_demand", color = "0.49 0.6 0.85", style="rounded"];
            45[label = "build_temperature_profiles\nscope: total", color = "0.32 0.6 0.85", style="rounded"];
            46[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.34 0.6 0.85", style="rounded"];
            47[label = "build_salt_cavern_potentials", color = "0.55 0.6 0.85", style="rounded"];
            48[label = "build_simplified_population_layouts", color = "0.46 0.6 0.85", style="rounded"];
            49[label = "build_industrial_energy_demand_per_node", color = "0.14 0.6 0.85", style="rounded"];
            50[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.27 0.6 0.85", style="rounded"];
            51[label = "build_industry_sector_ratios", color = "0.11 0.6 0.85", style="rounded"];
            52[label = "build_ammonia_production", color = "0.25 0.6 0.85", style="rounded"];
            53[label = "build_industrial_energy_demand_per_country_today", color = "0.44 0.6 0.85", style="rounded"];
            54[label = "build_industrial_production_per_country", color = "0.18 0.6 0.85", style="rounded"];
            55[label = "build_industrial_production_per_node", color = "0.41 0.6 0.85", style="rounded"];
            56[label = "build_industrial_distribution_key", color = "0.04 0.6 0.85", style="rounded"];
            57[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.09 0.6 0.85", style="rounded"];
            58[label = "build_industrial_energy_demand_per_node_today", color = "0.46 0.6 0.85", style="rounded"];
            59[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.39 0.6 0.85", style="rounded"];
            60[label = "build_temperature_profiles\nscope: rural", color = "0.32 0.6 0.85", style="rounded"];
            61[label = "build_temperature_profiles\nscope: urban", color = "0.32 0.6 0.85", style="rounded"];
            62[label = "build_cop_profiles", color = "0.55 0.6 0.85", style="rounded"];
            63[label = "build_solar_thermal_profiles\nscope: urban", color = "0.23 0.6 0.85", style="rounded"];
            64[label = "build_solar_thermal_profiles\nscope: rural", color = "0.23 0.6 0.85", style="rounded"];
            65[label = "plot_power_network_clustered", color = "0.41 0.6 0.85", style="rounded"];
            66[label = "plot_power_network", color = "0.40 0.6 0.85", style="rounded"];
            67[label = "plot_hydrogen_network", color = "0.42 0.6 0.85", style="rounded"];
            68[label = "plot_gas_network", color = "0.32 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            36 -> 1
            8 -> 1
            3 -> 2
            23 -> 2
            65 -> 2
            66 -> 2
            67 -> 2
            68 -> 2
            4 -> 3
            5 -> 4
            12 -> 4
            13 -> 4
            14 -> 4
            28 -> 4
            29 -> 4
            30 -> 4
            36 -> 4
            37 -> 4
            41 -> 4
            43 -> 4
            44 -> 4
            38 -> 4
            8 -> 4
            46 -> 4
            23 -> 4
            47 -> 4
            18 -> 4
            17 -> 4
            40 -> 4
            48 -> 4
            49 -> 4
            32 -> 4
            59 -> 4
            45 -> 4
            60 -> 4
            61 -> 4
            62 -> 4
            35 -> 4
            63 -> 4
            64 -> 4
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
            6 -> 13
            8 -> 13
            9 -> 13
            7 -> 13
            11 -> 13
            15 -> 14
            17 -> 14
            16 -> 15
            18 -> 17
            23 -> 17
            19 -> 18
            23 -> 18
            6 -> 18
            20 -> 19
            21 -> 19
            22 -> 19
            5 -> 19
            12 -> 19
            13 -> 19
            6 -> 19
            23 -> 19
            24 -> 19
            25 -> 19
            7 -> 19
            6 -> 20
            8 -> 20
            7 -> 20
            11 -> 20
            6 -> 21
            8 -> 21
            7 -> 21
            11 -> 21
            6 -> 22
            8 -> 22
            7 -> 22
            11 -> 22
            6 -> 24
            26 -> 25
            27 -> 25
            16 -> 28
            17 -> 28
            30 -> 29
            32 -> 29
            35 -> 29
            31 -> 30
            23 -> 30
            17 -> 31
            23 -> 31
            33 -> 32
            34 -> 33
            17 -> 33
            11 -> 33
            7 -> 34
            11 -> 34
            34 -> 35
            17 -> 35
            11 -> 35
            38 -> 37
            40 -> 37
            7 -> 38
            8 -> 38
            36 -> 38
            39 -> 38
            34 -> 40
            17 -> 40
            11 -> 40
            42 -> 41
            40 -> 41
            38 -> 42
            7 -> 43
            17 -> 43
            38 -> 43
            40 -> 44
            37 -> 44
            38 -> 44
            8 -> 44
            45 -> 44
            34 -> 45
            17 -> 45
            11 -> 45
            8 -> 46
            17 -> 46
            7 -> 46
            8 -> 47
            17 -> 47
            34 -> 48
            18 -> 48
            11 -> 48
            50 -> 49
            55 -> 49
            58 -> 49
            51 -> 50
            53 -> 50
            54 -> 50
            52 -> 51
            8 -> 51
            8 -> 52
            8 -> 53
            54 -> 53
            52 -> 54
            8 -> 54
            36 -> 54
            56 -> 55
            57 -> 55
            17 -> 56
            40 -> 56
            54 -> 57
            56 -> 58
            53 -> 58
            38 -> 59
            40 -> 59
            34 -> 60
            17 -> 60
            11 -> 60
            34 -> 61
            17 -> 61
            11 -> 61
            45 -> 62
            60 -> 62
            61 -> 62
            34 -> 63
            17 -> 63
            11 -> 63
            34 -> 64
            17 -> 64
            11 -> 64
            17 -> 65
            3 -> 66
            17 -> 66
            3 -> 67
            17 -> 67
            3 -> 68
            17 -> 68
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
            0[label = "all", color = "0.65 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.14 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.38 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network_myopic", color = "0.19 0.6 0.85", style="rounded"];
            4[label = "add_existing_baseyear", color = "0.32 0.6 0.85", style="rounded"];
            5[label = "prepare_sector_network\nsector_opts: ", color = "0.36 0.6 0.85", style="rounded"];
            6[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.57 0.6 0.85", style="rounded"];
            7[label = "base_network", color = "0.24 0.6 0.85", style="rounded"];
            8[label = "build_shapes", color = "0.51 0.6 0.85", style="rounded"];
            9[label = "retrieve_databundle", color = "0.31 0.6 0.85", style="rounded"];
            10[label = "build_ship_raster", color = "0.23 0.6 0.85", style="rounded"];
            11[label = "retrieve_ship_raster", color = "0.58 0.6 0.85", style="rounded"];
            12[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.21 0.6 0.85", style="rounded"];
            13[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.57 0.6 0.85", style="rounded"];
            14[label = "build_renewable_profiles\ntechnology: offwind-float", color = "0.57 0.6 0.85", style="rounded"];
            15[label = "cluster_gas_network", color = "0.06 0.6 0.85", style="rounded"];
            16[label = "build_gas_network", color = "0.29 0.6 0.85", style="rounded"];
            17[label = "retrieve_gas_infrastructure_data", color = "0.18 0.6 0.85", style="rounded"];
            18[label = "cluster_network\nclusters: 5", color = "0.59 0.6 0.85", style="rounded"];
            19[label = "simplify_network\nsimpl: ", color = "0.41 0.6 0.85", style="rounded"];
            20[label = "add_electricity", color = "0.25 0.6 0.85", style="rounded"];
            21[label = "build_renewable_profiles\ntechnology: solar", color = "0.57 0.6 0.85", style="rounded"];
            22[label = "build_renewable_profiles\ntechnology: solar-hsat", color = "0.57 0.6 0.85", style="rounded"];
            23[label = "build_renewable_profiles\ntechnology: onwind", color = "0.57 0.6 0.85", style="rounded"];
            24[label = "retrieve_cost_data\nyear: 2030", color = "0.41 0.6 0.85", style="rounded"];
            25[label = "build_powerplants", color = "0.36 0.6 0.85", style="rounded"];
            26[label = "build_electricity_demand", color = "0.12 0.6 0.85", style="rounded"];
            27[label = "retrieve_electricity_demand", color = "0.23 0.6 0.85", style="rounded"];
            28[label = "retrieve_synthetic_electricity_demand", color = "0.35 0.6 0.85", style="rounded"];
            29[label = "build_gas_input_locations", color = "0.10 0.6 0.85", style="rounded"];
            30[label = "time_aggregation", color = "0.60 0.6 0.85", style="rounded"];
            31[label = "prepare_network\nll: v1.5\nopts: ", color = "0.37 0.6 0.85", style="rounded"];
            32[label = "add_extra_components", color = "0.10 0.6 0.85", style="rounded"];
            33[label = "build_hourly_heat_demand", color = "0.11 0.6 0.85", style="rounded"];
            34[label = "build_daily_heat_demand\nscope: total", color = "0.39 0.6 0.85", style="rounded"];
            35[label = "build_population_layouts", color = "0.40 0.6 0.85", style="rounded"];
            36[label = "build_solar_thermal_profiles\nscope: total", color = "0.20 0.6 0.85", style="rounded"];
            37[label = "retrieve_eurostat_data", color = "0.45 0.6 0.85", style="rounded"];
            38[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.60 0.6 0.85", style="rounded"];
            39[label = "build_energy_totals", color = "0.48 0.6 0.85", style="rounded"];
            40[label = "retrieve_eurostat_household_data", color = "0.08 0.6 0.85", style="rounded"];
            41[label = "build_clustered_population_layouts", color = "0.07 0.6 0.85", style="rounded"];
            42[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.60 0.6 0.85", style="rounded"];
            43[label = "build_heat_totals", color = "0.52 0.6 0.85", style="rounded"];
            44[label = "build_shipping_demand", color = "0.50 0.6 0.85", style="rounded"];
            45[label = "build_transport_demand", color = "0.49 0.6 0.85", style="rounded"];
            46[label = "build_temperature_profiles\nscope: total", color = "0.63 0.6 0.85", style="rounded"];
            47[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.56 0.6 0.85", style="rounded"];
            48[label = "build_salt_cavern_potentials", color = "0.43 0.6 0.85", style="rounded"];
            49[label = "build_simplified_population_layouts", color = "0.58 0.6 0.85", style="rounded"];
            50[label = "build_industrial_energy_demand_per_node", color = "0.27 0.6 0.85", style="rounded"];
            51[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.26 0.6 0.85", style="rounded"];
            52[label = "build_industry_sector_ratios", color = "0.46 0.6 0.85", style="rounded"];
            53[label = "build_ammonia_production", color = "0.17 0.6 0.85", style="rounded"];
            54[label = "build_industrial_energy_demand_per_country_today", color = "0.02 0.6 0.85", style="rounded"];
            55[label = "build_industrial_production_per_country", color = "0.30 0.6 0.85", style="rounded"];
            56[label = "build_industrial_production_per_node", color = "0.30 0.6 0.85", style="rounded"];
            57[label = "build_industrial_distribution_key", color = "0.05 0.6 0.85", style="rounded"];
            58[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.01 0.6 0.85", style="rounded"];
            59[label = "build_industrial_energy_demand_per_node_today", color = "0.14 0.6 0.85", style="rounded"];
            60[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.56 0.6 0.85", style="rounded"];
            61[label = "build_temperature_profiles\nscope: rural", color = "0.63 0.6 0.85", style="rounded"];
            62[label = "build_temperature_profiles\nscope: urban", color = "0.63 0.6 0.85", style="rounded"];
            63[label = "build_cop_profiles", color = "0.64 0.6 0.85", style="rounded"];
            64[label = "build_solar_thermal_profiles\nscope: urban", color = "0.20 0.6 0.85", style="rounded"];
            65[label = "build_solar_thermal_profiles\nscope: rural", color = "0.20 0.6 0.85", style="rounded"];
            66[label = "build_existing_heating_distribution", color = "0.21 0.6 0.85", style="rounded"];
            67[label = "solve_sector_network_myopic", color = "0.19 0.6 0.85", style="rounded"];
            68[label = "add_brownfield", color = "0.27 0.6 0.85", style="rounded"];
            69[label = "prepare_sector_network\nsector_opts: ", color = "0.36 0.6 0.85", style="rounded"];
            70[label = "build_biomass_potentials\nplanning_horizons: 2040", color = "0.56 0.6 0.85", style="rounded"];
            71[label = "retrieve_cost_data\nyear: 2040", color = "0.41 0.6 0.85", style="rounded"];
            72[label = "build_industrial_energy_demand_per_node", color = "0.27 0.6 0.85", style="rounded"];
            73[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2040", color = "0.26 0.6 0.85", style="rounded"];
            74[label = "build_industrial_production_per_node", color = "0.30 0.6 0.85", style="rounded"];
            75[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2040", color = "0.01 0.6 0.85", style="rounded"];
            76[label = "build_district_heat_share\nplanning_horizons: 2040", color = "0.56 0.6 0.85", style="rounded"];
            77[label = "solve_sector_network_myopic", color = "0.19 0.6 0.85", style="rounded"];
            78[label = "add_brownfield", color = "0.27 0.6 0.85", style="rounded"];
            79[label = "prepare_sector_network\nsector_opts: ", color = "0.36 0.6 0.85", style="rounded"];
            80[label = "build_biomass_potentials\nplanning_horizons: 2050", color = "0.56 0.6 0.85", style="rounded"];
            81[label = "retrieve_cost_data\nyear: 2050", color = "0.41 0.6 0.85", style="rounded"];
            82[label = "build_industrial_energy_demand_per_node", color = "0.27 0.6 0.85", style="rounded"];
            83[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2050", color = "0.26 0.6 0.85", style="rounded"];
            84[label = "build_industrial_production_per_node", color = "0.30 0.6 0.85", style="rounded"];
            85[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2050", color = "0.01 0.6 0.85", style="rounded"];
            86[label = "build_district_heat_share\nplanning_horizons: 2050", color = "0.56 0.6 0.85", style="rounded"];
            87[label = "plot_power_network_clustered", color = "0.03 0.6 0.85", style="rounded"];
            88[label = "plot_power_network", color = "0.16 0.6 0.85", style="rounded"];
            89[label = "plot_power_network", color = "0.16 0.6 0.85", style="rounded"];
            90[label = "plot_power_network", color = "0.16 0.6 0.85", style="rounded"];
            91[label = "plot_hydrogen_network", color = "0.54 0.6 0.85", style="rounded"];
            92[label = "plot_hydrogen_network", color = "0.54 0.6 0.85", style="rounded"];
            93[label = "plot_hydrogen_network", color = "0.54 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            37 -> 1
            9 -> 1
            3 -> 2
            67 -> 2
            77 -> 2
            24 -> 2
            87 -> 2
            88 -> 2
            89 -> 2
            90 -> 2
            91 -> 2
            92 -> 2
            93 -> 2
            4 -> 3
            24 -> 3
            5 -> 4
            25 -> 4
            19 -> 4
            18 -> 4
            41 -> 4
            24 -> 4
            63 -> 4
            66 -> 4
            6 -> 5
            13 -> 5
            14 -> 5
            15 -> 5
            29 -> 5
            30 -> 5
            31 -> 5
            37 -> 5
            38 -> 5
            42 -> 5
            44 -> 5
            45 -> 5
            39 -> 5
            9 -> 5
            47 -> 5
            24 -> 5
            48 -> 5
            19 -> 5
            18 -> 5
            41 -> 5
            49 -> 5
            50 -> 5
            33 -> 5
            60 -> 5
            46 -> 5
            61 -> 5
            62 -> 5
            63 -> 5
            36 -> 5
            64 -> 5
            65 -> 5
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
            7 -> 14
            9 -> 14
            10 -> 14
            8 -> 14
            12 -> 14
            16 -> 15
            18 -> 15
            17 -> 16
            19 -> 18
            24 -> 18
            20 -> 19
            24 -> 19
            7 -> 19
            21 -> 20
            22 -> 20
            23 -> 20
            6 -> 20
            13 -> 20
            14 -> 20
            7 -> 20
            24 -> 20
            25 -> 20
            26 -> 20
            8 -> 20
            7 -> 21
            9 -> 21
            8 -> 21
            12 -> 21
            7 -> 22
            9 -> 22
            8 -> 22
            12 -> 22
            7 -> 23
            9 -> 23
            8 -> 23
            12 -> 23
            7 -> 25
            27 -> 26
            28 -> 26
            17 -> 29
            18 -> 29
            31 -> 30
            33 -> 30
            36 -> 30
            32 -> 31
            24 -> 31
            18 -> 32
            24 -> 32
            34 -> 33
            35 -> 34
            18 -> 34
            12 -> 34
            8 -> 35
            12 -> 35
            35 -> 36
            18 -> 36
            12 -> 36
            39 -> 38
            41 -> 38
            8 -> 39
            9 -> 39
            37 -> 39
            40 -> 39
            35 -> 41
            18 -> 41
            12 -> 41
            43 -> 42
            41 -> 42
            39 -> 43
            8 -> 44
            18 -> 44
            39 -> 44
            41 -> 45
            38 -> 45
            39 -> 45
            9 -> 45
            46 -> 45
            35 -> 46
            18 -> 46
            12 -> 46
            9 -> 47
            18 -> 47
            8 -> 47
            9 -> 48
            18 -> 48
            35 -> 49
            19 -> 49
            12 -> 49
            51 -> 50
            56 -> 50
            59 -> 50
            52 -> 51
            54 -> 51
            55 -> 51
            53 -> 52
            9 -> 52
            9 -> 53
            9 -> 54
            55 -> 54
            53 -> 55
            9 -> 55
            37 -> 55
            57 -> 56
            58 -> 56
            18 -> 57
            41 -> 57
            55 -> 58
            57 -> 59
            54 -> 59
            39 -> 60
            41 -> 60
            35 -> 61
            18 -> 61
            12 -> 61
            35 -> 62
            18 -> 62
            12 -> 62
            46 -> 63
            61 -> 63
            62 -> 63
            35 -> 64
            18 -> 64
            12 -> 64
            35 -> 65
            18 -> 65
            12 -> 65
            41 -> 66
            38 -> 66
            60 -> 66
            68 -> 67
            71 -> 67
            21 -> 68
            22 -> 68
            23 -> 68
            6 -> 68
            13 -> 68
            14 -> 68
            19 -> 68
            18 -> 68
            69 -> 68
            3 -> 68
            71 -> 68
            63 -> 68
            6 -> 69
            13 -> 69
            14 -> 69
            15 -> 69
            29 -> 69
            30 -> 69
            31 -> 69
            37 -> 69
            38 -> 69
            42 -> 69
            44 -> 69
            45 -> 69
            39 -> 69
            9 -> 69
            70 -> 69
            71 -> 69
            48 -> 69
            19 -> 69
            18 -> 69
            41 -> 69
            49 -> 69
            72 -> 69
            33 -> 69
            76 -> 69
            46 -> 69
            61 -> 69
            62 -> 69
            63 -> 69
            36 -> 69
            64 -> 69
            65 -> 69
            9 -> 70
            18 -> 70
            8 -> 70
            73 -> 72
            74 -> 72
            59 -> 72
            52 -> 73
            54 -> 73
            55 -> 73
            57 -> 74
            75 -> 74
            55 -> 75
            39 -> 76
            41 -> 76
            78 -> 77
            81 -> 77
            21 -> 78
            22 -> 78
            23 -> 78
            6 -> 78
            13 -> 78
            14 -> 78
            19 -> 78
            18 -> 78
            79 -> 78
            67 -> 78
            81 -> 78
            63 -> 78
            6 -> 79
            13 -> 79
            14 -> 79
            15 -> 79
            29 -> 79
            30 -> 79
            31 -> 79
            37 -> 79
            38 -> 79
            42 -> 79
            44 -> 79
            45 -> 79
            39 -> 79
            9 -> 79
            80 -> 79
            81 -> 79
            48 -> 79
            19 -> 79
            18 -> 79
            41 -> 79
            49 -> 79
            82 -> 79
            33 -> 79
            86 -> 79
            46 -> 79
            61 -> 79
            62 -> 79
            63 -> 79
            36 -> 79
            64 -> 79
            65 -> 79
            9 -> 80
            18 -> 80
            8 -> 80
            83 -> 82
            84 -> 82
            59 -> 82
            52 -> 83
            54 -> 83
            55 -> 83
            57 -> 84
            85 -> 84
            55 -> 85
            39 -> 86
            41 -> 86
            18 -> 87
            3 -> 88
            18 -> 88
            67 -> 89
            18 -> 89
            77 -> 90
            18 -> 90
            3 -> 91
            18 -> 91
            67 -> 92
            18 -> 92
            77 -> 93
            18 -> 93
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
