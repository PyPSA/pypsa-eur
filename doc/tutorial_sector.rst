.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

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
   :end-before: industry:

Documentation for all options will be added successively to :ref:`config`.

Scenarios can be defined like for electricity-only studies, but with additional configuration namespaces. Define scenario entries in ``config/scenarios.yaml`` and enable them via ``run.scenarios.enable: true`` to sweep different combinations of settings.  

Execution
---------

To run an overnight / greenfiled scenario with the specifications above, run

.. code:: console

    $ snakemake -call all --configfile config/test/config.overnight.yaml

which will result in the following jobs ``snakemake`` wants to run, some of
which were already included in the electricity-only tutorial:

.. code:: console

    job                                                 count
    ------------------------------------------------  -------
    add_transmission_projects_and_dlr                       1
    all                                                     1
    base_network                                            1
    build_ammonia_production                                1
    build_ates_potentials                                   1
    build_biomass_potentials                                1
    build_central_heating_temperature_profiles              1
    build_clustered_population_layouts                      1
    build_cop_profiles                                      1
    build_daily_heat_demand                                 1
    build_dh_areas                                          1
    build_direct_heat_source_utilisation_profiles           1
    build_district_heat_share                               1
    build_electricity_demand_base                           1
    build_energy_totals                                     1
    build_existing_heating_distribution                     1
    build_gas_input_locations                               1
    build_geothermal_heat_potential                         1
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
    build_osm_boundaries                                    4
    build_population_layouts                                1
    build_population_weighted_energy_totals                 2
    build_powerplants                                       1
    build_ptes_operations                                   1
    build_renewable_profiles                                6
    build_river_heat_potential                              1
    build_salt_cavern_potentials                            1
    build_sea_heat_potential                                1
    build_shapes                                            1
    build_ship_raster                                       1
    build_shipping_demand                                   1
    build_solar_rooftop_potentials                          1
    build_solar_thermal_profiles                            1
    build_temperature_profiles                              1
    build_transmission_projects                             1
    build_transport_demand                                  1
    chain_busmaps                                           1
    cluster_gas_network                                     1
    cluster_network                                         1
    compose_network                                         1
    determine_availability_matrix                           6
    make_summary                                            1
    plot_balance_map                                        7
    plot_balance_map_interactive                            7
    plot_balance_timeseries                                 1
    plot_base_network                                       1
    plot_clustered_network                                  1
    plot_cop_profiles                                       1
    plot_gas_network                                        1
    plot_heatmap_timeseries                                 1
    plot_hydrogen_network                                   1
    plot_power_network                                      1
    plot_summary                                            1
    process_cost_data                                       1
    retrieve_*                                             37
    simplify_network                                        1
    solve_network                                           1
    time_aggregation                                        1
    total                                                 126

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
            0[label = "all", color = "0.00 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.39 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.32 0.6 0.85", style="rounded"];
            3[label = "solve_network", color = "0.65 0.6 0.85", style="rounded"];
            4[label = "compose_network", color = "0.28 0.6 0.85", style="rounded"];
            5[label = "build_renewable_profiles", color = "0.20 0.6 0.85", style="rounded"];
            6[label = "determine_availability_matrix", color = "0.29 0.6 0.85", style="rounded"];
            7[label = "retrieve_corine", color = "0.44 0.6 0.85", style="rounded"];
            8[label = "retrieve_natura", color = "0.58 0.6 0.85", style="rounded"];
            9[label = "retrieve_luisa_land_cover", color = "0.56 0.6 0.85", style="rounded"];
            10[label = "build_shapes", color = "0.22 0.6 0.85", style="rounded"];
            11[label = "retrieve_eez", color = "0.47 0.6 0.85", style="rounded"];
            12[label = "retrieve_eu_nuts_2021", color = "0.48 0.6 0.85", style="rounded"];
            13[label = "build_osm_boundaries", color = "0.18 0.6 0.85", style="rounded"];
            14[label = "retrieve_osm_boundaries", color = "0.59 0.6 0.85", style="rounded"];
            15[label = "retrieve_jrc_ardeco", color = "0.55 0.6 0.85", style="rounded"];
            16[label = "retrieve_gdp_per_capita", color = "0.50 0.6 0.85", style="rounded"];
            17[label = "retrieve_population_count", color = "0.60 0.6 0.85", style="rounded"];
            18[label = "cluster_network", color = "0.27 0.6 0.85", style="rounded"];
            19[label = "simplify_network", color = "0.64 0.6 0.85", style="rounded"];
            20[label = "add_transmission_projects_and_dlr", color = "0.00 0.6 0.85", style="rounded"];
            21[label = "base_network", color = "0.01 0.6 0.85", style="rounded"];
            22[label = "retrieve_osm_archive", color = "0.59 0.6 0.85", style="rounded"];
            23[label = "build_transmission_projects", color = "0.25 0.6 0.85", style="rounded"];
            24[label = "build_electricity_demand_base", color = "0.09 0.6 0.85", style="rounded"];
            25[label = "build_electricity_demand", color = "0.08 0.6 0.85", style="rounded"];
            26[label = "retrieve_electricity_demand", color = "0.47 0.6 0.85", style="rounded"];
            27[label = "retrieve_synthetic_electricity_demand", color = "0.61 0.6 0.85", style="rounded"];
            28[label = "retrieve_cutout", color = "0.46 0.6 0.85", style="rounded"];
            29[label = "build_ship_raster", color = "0.23 0.6 0.85", style="rounded"];
            30[label = "retrieve_ship_raster", color = "0.61 0.6 0.85", style="rounded"];
            31[label = "process_cost_data", color = "0.40 0.6 0.85", style="rounded"];
            32[label = "retrieve_cost_data", color = "0.45 0.6 0.85", style="rounded"];
            33[label = "build_powerplants", color = "0.19 0.6 0.85", style="rounded"];
            34[label = "retrieve_powerplants", color = "0.60 0.6 0.85", style="rounded"];
            35[label = "build_monthly_prices", color = "0.17 0.6 0.85", style="rounded"];
            36[label = "retrieve_monthly_co2_prices", color = "0.57 0.6 0.85", style="rounded"];
            37[label = "retrieve_monthly_fuel_prices", color = "0.57 0.6 0.85", style="rounded"];
            38[label = "time_aggregation", color = "0.66 0.6 0.85", style="rounded"];
            39[label = "build_hourly_heat_demand", color = "0.12 0.6 0.85", style="rounded"];
            40[label = "build_daily_heat_demand", color = "0.06 0.6 0.85", style="rounded"];
            41[label = "build_population_layouts", color = "0.18 0.6 0.85", style="rounded"];
            42[label = "retrieve_worldbank_urban_population", color = "0.63 0.6 0.85", style="rounded"];
            43[label = "build_solar_thermal_profiles", color = "0.24 0.6 0.85", style="rounded"];
            44[label = "chain_busmaps", color = "0.26 0.6 0.85", style="rounded"];
            45[label = "build_solar_rooftop_potentials", color = "0.24 0.6 0.85", style="rounded"];
            46[label = "build_geothermal_heat_potential", color = "0.11 0.6 0.85", style="rounded"];
            47[label = "retrieve_geothermal_heat_utilisation_potentials", color = "0.52 0.6 0.85", style="rounded"];
            48[label = "retrieve_lau_regions", color = "0.55 0.6 0.85", style="rounded"];
            49[label = "build_river_heat_potential", color = "0.21 0.6 0.85", style="rounded"];
            50[label = "retrieve_hera_data_test_cutout", color = "0.54 0.6 0.85", style="rounded"];
            51[label = "build_dh_areas", color = "0.06 0.6 0.85", style="rounded"];
            52[label = "retrieve_dh_areas", color = "0.46 0.6 0.85", style="rounded"];
            53[label = "cluster_gas_network", color = "0.26 0.6 0.85", style="rounded"];
            54[label = "build_gas_network", color = "0.10 0.6 0.85", style="rounded"];
            55[label = "retrieve_gas_infrastructure_data", color = "0.50 0.6 0.85", style="rounded"];
            56[label = "build_gas_input_locations", color = "0.10 0.6 0.85", style="rounded"];
            57[label = "retrieve_gem_europe_gas_tracker", color = "0.51 0.6 0.85", style="rounded"];
            58[label = "retrieve_eurostat_balances", color = "0.49 0.6 0.85", style="rounded"];
            59[label = "build_population_weighted_energy_totals", color = "0.19 0.6 0.85", style="rounded"];
            60[label = "build_energy_totals", color = "0.09 0.6 0.85", style="rounded"];
            61[label = "retrieve_ghg_emissions", color = "0.52 0.6 0.85", style="rounded"];
            62[label = "retrieve_jrc_idees", color = "0.55 0.6 0.85", style="rounded"];
            63[label = "retrieve_eurostat_household_balances", color = "0.49 0.6 0.85", style="rounded"];
            64[label = "build_clustered_population_layouts", color = "0.05 0.6 0.85", style="rounded"];
            65[label = "build_heat_totals", color = "0.12 0.6 0.85", style="rounded"];
            66[label = "retrieve_country_hdd", color = "0.45 0.6 0.85", style="rounded"];
            67[label = "build_shipping_demand", color = "0.23 0.6 0.85", style="rounded"];
            68[label = "retrieve_attributed_ports", color = "0.42 0.6 0.85", style="rounded"];
            69[label = "build_transport_demand", color = "0.26 0.6 0.85", style="rounded"];
            70[label = "retrieve_mobility_profiles", color = "0.56 0.6 0.85", style="rounded"];
            71[label = "build_temperature_profiles", color = "0.25 0.6 0.85", style="rounded"];
            72[label = "build_biomass_potentials", color = "0.03 0.6 0.85", style="rounded"];
            73[label = "retrieve_enspreso_biomass", color = "0.48 0.6 0.85", style="rounded"];
            74[label = "retrieve_eu_nuts_2013", color = "0.48 0.6 0.85", style="rounded"];
            75[label = "retrieve_nuts3_population", color = "0.58 0.6 0.85", style="rounded"];
            76[label = "build_salt_cavern_potentials", color = "0.21 0.6 0.85", style="rounded"];
            77[label = "retrieve_h2_salt_caverns", color = "0.53 0.6 0.85", style="rounded"];
            78[label = "build_industrial_energy_demand_per_node", color = "0.14 0.6 0.85", style="rounded"];
            79[label = "build_industry_sector_ratios_intermediate", color = "0.16 0.6 0.85", style="rounded"];
            80[label = "build_industry_sector_ratios", color = "0.16 0.6 0.85", style="rounded"];
            81[label = "build_ammonia_production", color = "0.02 0.6 0.85", style="rounded"];
            82[label = "retrieve_nitrogen_statistics", color = "0.58 0.6 0.85", style="rounded"];
            83[label = "build_industrial_energy_demand_per_country_today", color = "0.13 0.6 0.85", style="rounded"];
            84[label = "build_industrial_production_per_country", color = "0.15 0.6 0.85", style="rounded"];
            85[label = "build_industrial_production_per_node", color = "0.16 0.6 0.85", style="rounded"];
            86[label = "build_industrial_distribution_key", color = "0.13 0.6 0.85", style="rounded"];
            87[label = "retrieve_hotmaps_industrial_sites", color = "0.54 0.6 0.85", style="rounded"];
            88[label = "retrieve_gem_steel_plant_tracker", color = "0.52 0.6 0.85", style="rounded"];
            89[label = "build_industrial_production_per_country_tomorrow", color = "0.15 0.6 0.85", style="rounded"];
            90[label = "build_industrial_energy_demand_per_node_today", color = "0.14 0.6 0.85", style="rounded"];
            91[label = "build_district_heat_share", color = "0.07 0.6 0.85", style="rounded"];
            92[label = "build_existing_heating_distribution", color = "0.10 0.6 0.85", style="rounded"];
            93[label = "build_cop_profiles", color = "0.06 0.6 0.85", style="rounded"];
            94[label = "build_sea_heat_potential", color = "0.22 0.6 0.85", style="rounded"];
            95[label = "retrieve_seawater_temperature", color = "0.61 0.6 0.85", style="rounded"];
            96[label = "build_central_heating_temperature_profiles", color = "0.04 0.6 0.85", style="rounded"];
            97[label = "build_ptes_operations", color = "0.19 0.6 0.85", style="rounded"];
            98[label = "build_direct_heat_source_utilisation_profiles", color = "0.07 0.6 0.85", style="rounded"];
            99[label = "build_ates_potentials", color = "0.02 0.6 0.85", style="rounded"];
            100[label = "retrieve_aquifer_data_bgr", color = "0.41 0.6 0.85", style="rounded"];
            101[label = "plot_base_network", color = "0.35 0.6 0.85", style="rounded"];
            102[label = "plot_clustered_network", color = "0.35 0.6 0.85", style="rounded"];
            103[label = "plot_power_network", color = "0.39 0.6 0.85", style="rounded"];
            104[label = "plot_cop_profiles", color = "0.36 0.6 0.85", style="rounded"];
            105[label = "plot_balance_timeseries", color = "0.34 0.6 0.85", style="rounded"];
            106[label = "plot_heatmap_timeseries", color = "0.37 0.6 0.85", style="rounded"];
            107[label = "plot_hydrogen_network", color = "0.38 0.6 0.85", style="rounded"];
            108[label = "plot_gas_network", color = "0.36 0.6 0.85", style="rounded"];
            109[label = "plot_balance_map", color = "0.32 0.6 0.85", style="rounded"];
            110[label = "plot_balance_map_interactive", color = "0.32 0.6 0.85", style="rounded"];
            105 -> 0
            108 -> 0
            107 -> 0
            109 -> 0
            101 -> 0
            103 -> 0
            104 -> 0
            1 -> 0
            106 -> 0
            102 -> 0
            110 -> 0
            61 -> 1
            2 -> 1
            58 -> 1
            3 -> 2
            4 -> 3
            46 -> 4
            99 -> 4
            97 -> 4
            67 -> 4
            60 -> 4
            45 -> 4
            18 -> 4
            69 -> 4
            76 -> 4
            93 -> 4
            71 -> 4
            19 -> 4
            59 -> 4
            32 -> 4
            33 -> 4
            49 -> 4
            31 -> 4
            35 -> 4
            64 -> 4
            38 -> 4
            44 -> 4
            78 -> 4
            43 -> 4
            58 -> 4
            24 -> 4
            85 -> 4
            91 -> 4
            98 -> 4
            39 -> 4
            53 -> 4
            72 -> 4
            56 -> 4
            92 -> 4
            5 -> 4
            10 -> 5
            28 -> 5
            18 -> 5
            6 -> 5
            7 -> 6
            28 -> 6
            18 -> 6
            9 -> 6
            8 -> 6
            29 -> 6
            10 -> 6
            15 -> 10
            12 -> 10
            13 -> 10
            16 -> 10
            11 -> 10
            17 -> 10
            11 -> 13
            14 -> 13
            24 -> 18
            19 -> 18
            21 -> 18
            20 -> 19
            21 -> 19
            23 -> 20
            21 -> 20
            10 -> 21
            22 -> 21
            10 -> 23
            21 -> 23
            10 -> 24
            19 -> 24
            25 -> 24
            27 -> 25
            26 -> 25
            28 -> 29
            30 -> 29
            19 -> 31
            32 -> 31
            34 -> 33
            18 -> 33
            36 -> 35
            37 -> 35
            39 -> 38
            18 -> 38
            43 -> 38
            40 -> 39
            28 -> 40
            18 -> 40
            41 -> 40
            10 -> 41
            28 -> 41
            42 -> 41
            28 -> 43
            18 -> 43
            41 -> 43
            18 -> 44
            19 -> 44
            28 -> 45
            5 -> 45
            41 -> 45
            48 -> 46
            18 -> 46
            47 -> 46
            51 -> 49
            50 -> 49
            18 -> 49
            18 -> 51
            52 -> 51
            54 -> 53
            18 -> 53
            55 -> 54
            57 -> 56
            18 -> 56
            55 -> 56
            60 -> 59
            64 -> 59
            65 -> 59
            58 -> 60
            62 -> 60
            63 -> 60
            61 -> 60
            10 -> 60
            28 -> 64
            18 -> 64
            41 -> 64
            60 -> 65
            66 -> 65
            60 -> 67
            10 -> 67
            18 -> 67
            68 -> 67
            60 -> 69
            70 -> 69
            18 -> 69
            64 -> 69
            71 -> 69
            59 -> 69
            28 -> 71
            18 -> 71
            41 -> 71
            18 -> 72
            73 -> 72
            58 -> 72
            75 -> 72
            74 -> 72
            10 -> 72
            77 -> 76
            18 -> 76
            79 -> 78
            90 -> 78
            85 -> 78
            83 -> 79
            80 -> 79
            84 -> 79
            62 -> 80
            81 -> 80
            82 -> 81
            60 -> 83
            62 -> 83
            84 -> 83
            58 -> 84
            62 -> 84
            81 -> 84
            89 -> 85
            86 -> 85
            64 -> 86
            87 -> 86
            18 -> 86
            88 -> 86
            84 -> 89
            83 -> 90
            86 -> 90
            60 -> 91
            64 -> 91
            64 -> 92
            59 -> 92
            91 -> 92
            18 -> 93
            49 -> 93
            96 -> 93
            94 -> 93
            71 -> 93
            97 -> 93
            51 -> 94
            95 -> 94
            18 -> 94
            71 -> 96
            18 -> 96
            18 -> 97
            96 -> 97
            96 -> 98
            51 -> 99
            18 -> 99
            100 -> 99
            96 -> 99
            21 -> 101
            18 -> 102
            3 -> 103
            18 -> 103
            93 -> 104
            3 -> 105
            3 -> 106
            3 -> 107
            18 -> 107
            3 -> 108
            18 -> 108
            3 -> 109
            18 -> 109
            3 -> 110
            18 -> 110
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
configuration namespaces. The myopic foresight mode relies on the top-level
``planning_horizons`` list to define the sequence of investment horizons:

.. literalinclude:: ../config/test/config.myopic.yaml
   :language: yaml
   :start-at: planning_horizons:
   :end-before: countries:

For allowed wildcard values (e.g., ``{horizon}``), refer to :ref:`wildcards`.

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

.. code:: console

    $ snakemake -call all --configfile config/test/config.myopic.yaml

which will result in additional jobs ``snakemake`` wants to run, which
translates to the following job summary. Note how the rules
``compose_network`` and ``solve_network`` are now called multiple times (once
per planning horizon), while most other rules remain the same as in the
overnight case:

.. code:: console

    job                                                 count
    ------------------------------------------------  -------
    add_transmission_projects_and_dlr                       1
    all                                                     1
    base_network                                            1
    build_ammonia_production                                1
    build_ates_potentials                                   3
    build_biomass_potentials                                3
    build_central_heating_temperature_profiles              3
    build_clustered_population_layouts                      1
    build_cop_profiles                                      3
    build_daily_heat_demand                                 1
    build_dh_areas                                          1
    build_direct_heat_source_utilisation_profiles           3
    build_district_heat_share                               3
    build_electricity_demand_base                           1
    build_energy_totals                                     1
    build_existing_heating_distribution                     3
    build_gas_input_locations                               1
    build_geothermal_heat_potential                         1
    build_heat_totals                                       1
    build_hourly_heat_demand                                1
    build_industrial_distribution_key                       1
    build_industrial_energy_demand_per_country_today        1
    build_industrial_energy_demand_per_node                 3
    build_industrial_energy_demand_per_node_today           1
    build_industrial_production_per_country                 1
    build_industrial_production_per_country_tomorrow        3
    build_industrial_production_per_node                    3
    build_industry_sector_ratios                            1
    build_industry_sector_ratios_intermediate               3
    ...
    chain_busmaps                                           1
    cluster_gas_network                                     1
    cluster_network                                         1
    compose_network                                         3
    determine_availability_matrix                           6
    make_summary                                            1
    ...
    simplify_network                                        1
    solve_network                                           3
    time_aggregation                                        1
    total                                                 197

This sequential workflow nicely outlines how the pathway optimisation with
myopic foresight is implemented: each planning horizon depends on the solved
network from the previous horizon, which allows for warm-starting the
optimisation and tracking existing infrastructure over time.

You can visualise the full DAG with::

    snakemake --dag all --configfile config/test/config.myopic.yaml | dot -Tpng -o myopic-dag.png

|


Scaling-Up
==========

If you now feel confident and want to tackle runs with larger temporal, technological and
spatial scope, clean-up the repository and after modifying the ``config/config.yaml`` file
target the collection rule ``all`` again without providing the test
configuration file.

.. code:: console

    $ snakemake -call purge
    $ snakemake -call all

.. note::

    It is good practice to perform a dry-run using the option `-n`, before you
    commit to a run:

    .. code:: console

        $ snakemake -call all -n
