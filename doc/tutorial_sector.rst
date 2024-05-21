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
    build_renewable_profiles                                5
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
    retrieve_gas_infrastructure_data                        1
    retrieve_ship_raster                                    1
    retrieve_synthetic_electricity_demand                   1
    simplify_network                                        1
    solve_sector_network                                    1
    time_aggregation                                        1
    total                                                  67

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
            0[label = "all", color = "0.06 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.16 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.57 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network", color = "0.35 0.6 0.85", style="rounded"];
            4[label = "prepare_sector_network", color = "0.36 0.6 0.85", style="rounded"];
            5[label = "build_renewable_profiles", color = "0.12 0.6 0.85", style="rounded"];
            6[label = "base_network", color = "0.45 0.6 0.85", style="rounded"];
            7[label = "build_shapes", color = "0.66 0.6 0.85", style="rounded"];
            8[label = "retrieve_databundle", color = "0.28 0.6 0.85", style="rounded"];
            9[label = "build_ship_raster", color = "0.58 0.6 0.85", style="rounded"];
            10[label = "retrieve_ship_raster", color = "0.50 0.6 0.85", style="rounded"];
            11[label = "retrieve_cutout", color = "0.34 0.6 0.85", style="rounded"];
            12[label = "cluster_gas_network", color = "0.09 0.6 0.85", style="rounded"];
            13[label = "build_gas_network", color = "0.19 0.6 0.85", style="rounded"];
            14[label = "retrieve_gas_infrastructure_data", color = "0.52 0.6 0.85", style="rounded"];
            15[label = "cluster_network", color = "0.05 0.6 0.85", style="rounded"];
            16[label = "simplify_network", color = "0.50 0.6 0.85", style="rounded"];
            17[label = "add_electricity", color = "0.32 0.6 0.85", style="rounded"];
            18[label = "retrieve_cost_data", color = "0.39 0.6 0.85", style="rounded"];
            19[label = "build_powerplants", color = "0.33 0.6 0.85", style="rounded"];
            20[label = "build_electricity_demand", color = "0.31 0.6 0.85", style="rounded"];
            21[label = "retrieve_electricity_demand", color = "0.23 0.6 0.85", style="rounded"];
            22[label = "retrieve_synthetic_electricity_demand", color = "0.16 0.6 0.85", style="rounded"];
            23[label = "build_gas_input_locations", color = "0.62 0.6 0.85", style="rounded"];
            24[label = "time_aggregation", color = "0.28 0.6 0.85", style="rounded"];
            25[label = "prepare_network", color = "0.23 0.6 0.85", style="rounded"];
            26[label = "add_extra_components", color = "0.21 0.6 0.85", style="rounded"];
            27[label = "build_hourly_heat_demand", color = "0.21 0.6 0.85", style="rounded"];
            28[label = "build_daily_heat_demand", color = "0.24 0.6 0.85", style="rounded"];
            29[label = "build_population_layouts", color = "0.33 0.6 0.85", style="rounded"];
            30[label = "build_solar_thermal_profiles", color = "0.60 0.6 0.85", style="rounded"];
            31[label = "retrieve_eurostat_data", color = "0.03 0.6 0.85", style="rounded"];
            32[label = "build_population_weighted_energy_totals", color = "0.11 0.6 0.85", style="rounded"];
            33[label = "build_energy_totals", color = "0.62 0.6 0.85", style="rounded"];
            34[label = "build_clustered_population_layouts", color = "0.55 0.6 0.85", style="rounded"];
            35[label = "build_heat_totals", color = "0.38 0.6 0.85", style="rounded"];
            36[label = "build_shipping_demand", color = "0.37 0.6 0.85", style="rounded"];
            37[label = "build_transport_demand", color = "0.44 0.6 0.85", style="rounded"];
            38[label = "build_temperature_profiles", color = "0.54 0.6 0.85", style="rounded"];
            39[label = "build_biomass_potentials", color = "0.00 0.6 0.85", style="rounded"];
            40[label = "build_salt_cavern_potentials", color = "0.61 0.6 0.85", style="rounded"];
            41[label = "build_simplified_population_layouts", color = "0.40 0.6 0.85", style="rounded"];
            42[label = "build_industrial_energy_demand_per_node", color = "0.22 0.6 0.85", style="rounded"];
            43[label = "build_industry_sector_ratios_intermediate", color = "0.65 0.6 0.85", style="rounded"];
            44[label = "build_industry_sector_ratios", color = "0.57 0.6 0.85", style="rounded"];
            45[label = "build_ammonia_production", color = "0.01 0.6 0.85", style="rounded"];
            46[label = "build_industrial_energy_demand_per_country_today", color = "0.01 0.6 0.85", style="rounded"];
            47[label = "build_industrial_production_per_country", color = "0.59 0.6 0.85", style="rounded"];
            48[label = "build_industrial_production_per_node", color = "0.47 0.6 0.85", style="rounded"];
            49[label = "build_industrial_distribution_key", color = "0.09 0.6 0.85", style="rounded"];
            50[label = "build_industrial_production_per_country_tomorrow", color = "0.43 0.6 0.85", style="rounded"];
            51[label = "build_industrial_energy_demand_per_node_today", color = "0.64 0.6 0.85", style="rounded"];
            52[label = "build_district_heat_share", color = "0.38 0.6 0.85", style="rounded"];
            53[label = "build_cop_profiles", color = "0.14 0.6 0.85", style="rounded"];
            54[label = "plot_power_network_clustered", color = "0.25 0.6 0.85", style="rounded"];
            55[label = "plot_power_network", color = "0.65 0.6 0.85", style="rounded"];
            56[label = "plot_hydrogen_network", color = "0.30 0.6 0.85", style="rounded"];
            57[label = "plot_gas_network", color = "0.56 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            31 -> 1
            8 -> 1
            57 -> 2
            18 -> 2
            55 -> 2
            3 -> 2
            56 -> 2
            54 -> 2
            4 -> 3
            39 -> 4
            16 -> 4
            42 -> 4
            31 -> 4
            15 -> 4
            25 -> 4
            38 -> 4
            34 -> 4
            36 -> 4
            12 -> 4
            52 -> 4
            18 -> 4
            32 -> 4
            5 -> 4
            30 -> 4
            40 -> 4
            24 -> 4
            33 -> 4
            53 -> 4
            41 -> 4
            23 -> 4
            8 -> 4
            37 -> 4
            27 -> 4
            11 -> 5
            6 -> 5
            9 -> 5
            8 -> 5
            7 -> 5
            7 -> 6
            8 -> 7
            11 -> 9
            10 -> 9
            13 -> 12
            15 -> 12
            14 -> 13
            16 -> 15
            18 -> 15
            6 -> 16
            17 -> 16
            18 -> 16
            20 -> 17
            6 -> 17
            18 -> 17
            5 -> 17
            19 -> 17
            7 -> 17
            6 -> 19
            21 -> 20
            22 -> 20
            14 -> 23
            15 -> 23
            30 -> 24
            27 -> 24
            25 -> 24
            26 -> 25
            18 -> 25
            15 -> 26
            18 -> 26
            28 -> 27
            11 -> 28
            15 -> 28
            29 -> 28
            11 -> 29
            7 -> 29
            11 -> 30
            15 -> 30
            29 -> 30
            33 -> 32
            34 -> 32
            35 -> 32
            31 -> 33
            7 -> 33
            8 -> 33
            11 -> 34
            15 -> 34
            29 -> 34
            33 -> 35
            33 -> 36
            7 -> 36
            15 -> 36
            34 -> 37
            33 -> 37
            32 -> 37
            8 -> 37
            38 -> 37
            11 -> 38
            15 -> 38
            29 -> 38
            15 -> 39
            7 -> 39
            8 -> 39
            15 -> 40
            8 -> 40
            11 -> 41
            16 -> 41
            29 -> 41
            48 -> 42
            51 -> 42
            43 -> 42
            44 -> 43
            47 -> 43
            46 -> 43
            45 -> 44
            8 -> 44
            8 -> 45
            47 -> 46
            8 -> 46
            45 -> 47
            31 -> 47
            8 -> 47
            50 -> 48
            49 -> 48
            34 -> 49
            15 -> 49
            47 -> 50
            49 -> 51
            46 -> 51
            33 -> 52
            34 -> 52
            38 -> 53
            15 -> 54
            3 -> 55
            15 -> 55
            3 -> 56
            15 -> 56
            3 -> 57
            15 -> 57
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
            0[label = "all", color = "0.10 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.56 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.44 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network_myopic", color = "0.11 0.6 0.85", style="rounded"];
            4[label = "add_existing_baseyear", color = "0.24 0.6 0.85", style="rounded"];
            5[label = "prepare_sector_network", color = "0.26 0.6 0.85", style="rounded"];
            6[label = "build_renewable_profiles", color = "0.62 0.6 0.85", style="rounded"];
            7[label = "base_network", color = "0.13 0.6 0.85", style="rounded"];
            8[label = "build_shapes", color = "0.31 0.6 0.85", style="rounded"];
            9[label = "retrieve_databundle", color = "0.49 0.6 0.85", style="rounded"];
            10[label = "build_ship_raster", color = "0.39 0.6 0.85", style="rounded"];
            11[label = "retrieve_ship_raster", color = "0.61 0.6 0.85", style="rounded"];
            12[label = "retrieve_cutout", color = "0.65 0.6 0.85", style="rounded"];
            13[label = "cluster_gas_network", color = "0.05 0.6 0.85", style="rounded"];
            14[label = "build_gas_network", color = "0.18 0.6 0.85", style="rounded"];
            15[label = "retrieve_gas_infrastructure_data", color = "0.49 0.6 0.85", style="rounded"];
            16[label = "cluster_network", color = "0.45 0.6 0.85", style="rounded"];
            17[label = "simplify_network", color = "0.28 0.6 0.85", style="rounded"];
            18[label = "add_electricity", color = "0.40 0.6 0.85", style="rounded"];
            19[label = "retrieve_cost_data", color = "0.66 0.6 0.85", style="rounded"];
            20[label = "build_powerplants", color = "0.60 0.6 0.85", style="rounded"];
            21[label = "build_electricity_demand", color = "0.48 0.6 0.85", style="rounded"];
            22[label = "retrieve_electricity_demand", color = "0.60 0.6 0.85", style="rounded"];
            23[label = "retrieve_synthetic_electricity_demand", color = "0.09 0.6 0.85", style="rounded"];
            24[label = "build_gas_input_locations", color = "0.08 0.6 0.85", style="rounded"];
            25[label = "time_aggregation", color = "0.51 0.6 0.85", style="rounded"];
            26[label = "prepare_network", color = "0.58 0.6 0.85", style="rounded"];
            27[label = "add_extra_components", color = "0.25 0.6 0.85", style="rounded"];
            28[label = "build_hourly_heat_demand", color = "0.06 0.6 0.85", style="rounded"];
            29[label = "build_daily_heat_demand", color = "0.23 0.6 0.85", style="rounded"];
            30[label = "build_population_layouts", color = "0.33 0.6 0.85", style="rounded"];
            31[label = "build_solar_thermal_profiles", color = "0.62 0.6 0.85", style="rounded"];
            32[label = "retrieve_eurostat_data", color = "0.43 0.6 0.85", style="rounded"];
            33[label = "build_population_weighted_energy_totals", color = "0.12 0.6 0.85", style="rounded"];
            34[label = "build_energy_totals", color = "0.17 0.6 0.85", style="rounded"];
            35[label = "build_clustered_population_layouts", color = "0.59 0.6 0.85", style="rounded"];
            36[label = "build_heat_totals", color = "0.01 0.6 0.85", style="rounded"];
            37[label = "build_shipping_demand", color = "0.15 0.6 0.85", style="rounded"];
            38[label = "build_transport_demand", color = "0.16 0.6 0.85", style="rounded"];
            39[label = "build_temperature_profiles", color = "0.41 0.6 0.85", style="rounded"];
            40[label = "build_biomass_potentials", color = "0.53 0.6 0.85", style="rounded"];
            41[label = "build_salt_cavern_potentials", color = "0.54 0.6 0.85", style="rounded"];
            42[label = "build_simplified_population_layouts", color = "0.42 0.6 0.85", style="rounded"];
            43[label = "build_industrial_energy_demand_per_node", color = "0.28 0.6 0.85", style="rounded"];
            44[label = "build_industry_sector_ratios_intermediate", color = "0.35 0.6 0.85", style="rounded"];
            45[label = "build_industry_sector_ratios", color = "0.58 0.6 0.85", style="rounded"];
            46[label = "build_ammonia_production", color = "0.02 0.6 0.85", style="rounded"];
            47[label = "build_industrial_energy_demand_per_country_today", color = "0.07 0.6 0.85", style="rounded"];
            48[label = "build_industrial_production_per_country", color = "0.40 0.6 0.85", style="rounded"];
            49[label = "build_industrial_production_per_node", color = "0.46 0.6 0.85", style="rounded"];
            50[label = "build_industrial_distribution_key", color = "0.20 0.6 0.85", style="rounded"];
            51[label = "build_industrial_production_per_country_tomorrow", color = "0.47 0.6 0.85", style="rounded"];
            52[label = "build_industrial_energy_demand_per_node_today", color = "0.30 0.6 0.85", style="rounded"];
            53[label = "build_district_heat_share", color = "0.17 0.6 0.85", style="rounded"];
            54[label = "build_cop_profiles", color = "0.14 0.6 0.85", style="rounded"];
            55[label = "build_existing_heating_distribution", color = "0.19 0.6 0.85", style="rounded"];
            56[label = "add_brownfield", color = "0.33 0.6 0.85", style="rounded"];
            57[label = "plot_power_network_clustered", color = "0.56 0.6 0.85", style="rounded"];
            58[label = "plot_power_network", color = "0.51 0.6 0.85", style="rounded"];
            59[label = "plot_hydrogen_network", color = "0.47 0.6 0.85", style="rounded"];
            1 -> 0
            9 -> 1
            32 -> 1
            2 -> 1
            58 -> 2
            3 -> 2
            59 -> 2
            19 -> 2
            57 -> 2
            19 -> 3
            4 -> 3
            56 -> 3
            17 -> 4
            35 -> 4
            16 -> 4
            20 -> 4
            54 -> 4
            55 -> 4
            19 -> 4
            5 -> 4
            25 -> 5
            33 -> 5
            54 -> 5
            40 -> 5
            37 -> 5
            41 -> 5
            38 -> 5
            53 -> 5
            34 -> 5
            39 -> 5
            13 -> 5
            42 -> 5
            26 -> 5
            32 -> 5
            28 -> 5
            17 -> 5
            16 -> 5
            35 -> 5
            24 -> 5
            43 -> 5
            6 -> 5
            31 -> 5
            9 -> 5
            19 -> 5
            8 -> 6
            7 -> 6
            9 -> 6
            12 -> 6
            10 -> 6
            8 -> 7
            9 -> 8
            12 -> 10
            11 -> 10
            14 -> 13
            16 -> 13
            15 -> 14
            17 -> 16
            19 -> 16
            18 -> 17
            7 -> 17
            19 -> 17
            8 -> 18
            21 -> 18
            7 -> 18
            20 -> 18
            19 -> 18
            6 -> 18
            7 -> 20
            22 -> 21
            23 -> 21
            15 -> 24
            16 -> 24
            31 -> 25
            26 -> 25
            28 -> 25
            19 -> 26
            27 -> 26
            19 -> 27
            16 -> 27
            29 -> 28
            12 -> 29
            16 -> 29
            30 -> 29
            12 -> 30
            8 -> 30
            12 -> 31
            16 -> 31
            30 -> 31
            34 -> 33
            36 -> 33
            35 -> 33
            9 -> 34
            32 -> 34
            8 -> 34
            12 -> 35
            16 -> 35
            30 -> 35
            34 -> 36
            34 -> 37
            8 -> 37
            16 -> 37
            39 -> 38
            33 -> 38
            35 -> 38
            9 -> 38
            34 -> 38
            12 -> 39
            16 -> 39
            30 -> 39
            9 -> 40
            8 -> 40
            16 -> 40
            9 -> 41
            16 -> 41
            12 -> 42
            17 -> 42
            30 -> 42
            49 -> 43
            52 -> 43
            44 -> 43
            45 -> 44
            48 -> 44
            47 -> 44
            9 -> 45
            46 -> 45
            9 -> 46
            9 -> 47
            48 -> 47
            9 -> 48
            32 -> 48
            46 -> 48
            50 -> 49
            51 -> 49
            35 -> 50
            16 -> 50
            48 -> 51
            47 -> 52
            50 -> 52
            34 -> 53
            35 -> 53
            39 -> 54
            53 -> 55
            35 -> 55
            33 -> 55
            17 -> 56
            3 -> 56
            16 -> 56
            54 -> 56
            19 -> 56
            5 -> 56
            6 -> 56
            16 -> 57
            16 -> 58
            3 -> 58
            16 -> 59
            3 -> 59
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
