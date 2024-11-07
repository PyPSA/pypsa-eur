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

.. code:: console

    $ snakemake all --configfile config/test/config.overnight.yaml

which will result in the following jobs ``snakemake`` wants to run, some of
which were already included in the electricity-only tutorial:

.. code:: console

    job                                                 count
    ------------------------------------------------  -------
    add_electricity                                         1
    add_transmission_projects_and_dlr                       1
    all                                                     1
    base_network                                            1
    build_ammonia_production                                1
    build_biomass_potentials                                1
    build_central_heating_temperature_profiles              1
    build_clustered_population_layouts                      1
    build_cop_profiles                                      1
    build_daily_heat_demand                                 1
    build_district_heat_share                               1
    build_electricity_demand                                1
    build_electricity_demand_base                           1
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
    build_solar_thermal_profiles                            1
    build_temperature_profiles                              1
    build_transmission_projects                             1
    build_transport_demand                                  1
    cluster_gas_network                                     1
    cluster_network                                         1
    determine_availability_matrix                           6
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
    retrieve_eez                                            1
    retrieve_electricity_demand                             1
    retrieve_eurostat_data                                  1
    retrieve_eurostat_household_data                        1
    retrieve_gas_infrastructure_data                        1
    retrieve_gem_europe_gas_tracker                         1
    retrieve_gem_steel_plant_tracker                        1
    retrieve_hotmaps_industrial_sites                       1
    retrieve_jrc_enspreso_biomass                           1
    retrieve_jrc_idees                                      1
    retrieve_naturalearth_countries                         1
    retrieve_osm_prebuilt                                   1
    retrieve_ship_raster                                    1
    retrieve_synthetic_electricity_demand                   1
    retrieve_usgs_ammonia_production                        1
    retrieve_worldbank_urban_population                     1
    simplify_network                                        1
    solve_sector_network                                    1
    time_aggregation                                        1
    total                                                  83

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
            0[label = "all", color = "0.41 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.60 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.62 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network", color = "0.62 0.6 0.85", style="rounded"];
            4[label = "prepare_sector_network\nsector_opts: ", color = "0.45 0.6 0.85", style="rounded"];
            5[label = "build_renewable_profiles", color = "0.20 0.6 0.85", style="rounded"];
            6[label = "determine_availability_matrix\ntechnology: offwind-ac", color = "0.24 0.6 0.85", style="rounded"];
            7[label = "retrieve_databundle", color = "0.58 0.6 0.85", style="rounded"];
            8[label = "build_ship_raster", color = "0.51 0.6 0.85", style="rounded"];
            9[label = "retrieve_ship_raster", color = "0.03 0.6 0.85", style="rounded"];
            10[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.32 0.6 0.85", style="rounded"];
            11[label = "build_shapes", color = "0.11 0.6 0.85", style="rounded"];
            12[label = "retrieve_naturalearth_countries", color = "0.63 0.6 0.85", style="rounded"];
            13[label = "retrieve_eez", color = "0.00 0.6 0.85", style="rounded"];
            14[label = "cluster_network\nclusters: 5", color = "0.64 0.6 0.85", style="rounded"];
            15[label = "simplify_network", color = "0.21 0.6 0.85", style="rounded"];
            16[label = "add_transmission_projects_and_dlr", color = "0.17 0.6 0.85", style="rounded"];
            17[label = "base_network", color = "0.53 0.6 0.85", style="rounded"];
            18[label = "retrieve_osm_prebuilt", color = "0.21 0.6 0.85", style="rounded"];
            19[label = "build_transmission_projects", color = "0.02 0.6 0.85", style="rounded"];
            20[label = "build_electricity_demand_base", color = "0.44 0.6 0.85", style="rounded"];
            21[label = "build_electricity_demand", color = "0.16 0.6 0.85", style="rounded"];
            22[label = "retrieve_electricity_demand", color = "0.06 0.6 0.85", style="rounded"];
            23[label = "retrieve_synthetic_electricity_demand", color = "0.09 0.6 0.85", style="rounded"];
            24[label = "build_renewable_profiles", color = "0.20 0.6 0.85", style="rounded"];
            25[label = "determine_availability_matrix\ntechnology: offwind-dc", color = "0.24 0.6 0.85", style="rounded"];
            26[label = "build_renewable_profiles", color = "0.20 0.6 0.85", style="rounded"];
            27[label = "determine_availability_matrix\ntechnology: offwind-float", color = "0.24 0.6 0.85", style="rounded"];
            28[label = "cluster_gas_network", color = "0.39 0.6 0.85", style="rounded"];
            29[label = "build_gas_network", color = "0.29 0.6 0.85", style="rounded"];
            30[label = "retrieve_gas_infrastructure_data", color = "0.25 0.6 0.85", style="rounded"];
            31[label = "build_gas_input_locations", color = "0.58 0.6 0.85", style="rounded"];
            32[label = "retrieve_gem_europe_gas_tracker", color = "0.05 0.6 0.85", style="rounded"];
            33[label = "time_aggregation", color = "0.66 0.6 0.85", style="rounded"];
            34[label = "prepare_network\nll: v1.5\nopts: ", color = "0.55 0.6 0.85", style="rounded"];
            35[label = "add_electricity", color = "0.36 0.6 0.85", style="rounded"];
            36[label = "build_renewable_profiles", color = "0.20 0.6 0.85", style="rounded"];
            37[label = "determine_availability_matrix\ntechnology: solar", color = "0.24 0.6 0.85", style="rounded"];
            38[label = "build_renewable_profiles", color = "0.20 0.6 0.85", style="rounded"];
            39[label = "determine_availability_matrix\ntechnology: solar-hsat", color = "0.24 0.6 0.85", style="rounded"];
            40[label = "build_renewable_profiles", color = "0.20 0.6 0.85", style="rounded"];
            41[label = "determine_availability_matrix\ntechnology: onwind", color = "0.24 0.6 0.85", style="rounded"];
            42[label = "retrieve_cost_data\nyear: 2030", color = "0.55 0.6 0.85", style="rounded"];
            43[label = "build_powerplants", color = "0.18 0.6 0.85", style="rounded"];
            44[label = "build_hourly_heat_demand", color = "0.29 0.6 0.85", style="rounded"];
            45[label = "build_daily_heat_demand", color = "0.40 0.6 0.85", style="rounded"];
            46[label = "build_population_layouts", color = "0.27 0.6 0.85", style="rounded"];
            47[label = "retrieve_worldbank_urban_population", color = "0.30 0.6 0.85", style="rounded"];
            48[label = "build_solar_thermal_profiles", color = "0.27 0.6 0.85", style="rounded"];
            49[label = "retrieve_eurostat_data", color = "0.13 0.6 0.85", style="rounded"];
            50[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.24 0.6 0.85", style="rounded"];
            51[label = "build_energy_totals", color = "0.26 0.6 0.85", style="rounded"];
            52[label = "retrieve_jrc_idees", color = "0.48 0.6 0.85", style="rounded"];
            53[label = "retrieve_eurostat_household_data", color = "0.12 0.6 0.85", style="rounded"];
            54[label = "build_clustered_population_layouts", color = "0.35 0.6 0.85", style="rounded"];
            55[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.24 0.6 0.85", style="rounded"];
            56[label = "build_heat_totals", color = "0.01 0.6 0.85", style="rounded"];
            57[label = "build_shipping_demand", color = "0.60 0.6 0.85", style="rounded"];
            58[label = "build_transport_demand", color = "0.50 0.6 0.85", style="rounded"];
            59[label = "build_temperature_profiles", color = "0.54 0.6 0.85", style="rounded"];
            60[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.45 0.6 0.85", style="rounded"];
            61[label = "retrieve_jrc_enspreso_biomass", color = "0.07 0.6 0.85", style="rounded"];
            62[label = "build_salt_cavern_potentials", color = "0.18 0.6 0.85", style="rounded"];
            63[label = "build_industrial_energy_demand_per_node", color = "0.65 0.6 0.85", style="rounded"];
            64[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.64 0.6 0.85", style="rounded"];
            65[label = "build_industry_sector_ratios", color = "0.51 0.6 0.85", style="rounded"];
            66[label = "build_ammonia_production", color = "0.15 0.6 0.85", style="rounded"];
            67[label = "retrieve_usgs_ammonia_production", color = "0.38 0.6 0.85", style="rounded"];
            68[label = "build_industrial_energy_demand_per_country_today", color = "0.65 0.6 0.85", style="rounded"];
            69[label = "build_industrial_production_per_country", color = "0.11 0.6 0.85", style="rounded"];
            70[label = "build_industrial_production_per_node", color = "0.07 0.6 0.85", style="rounded"];
            71[label = "build_industrial_distribution_key", color = "0.48 0.6 0.85", style="rounded"];
            72[label = "retrieve_hotmaps_industrial_sites", color = "0.20 0.6 0.85", style="rounded"];
            73[label = "retrieve_gem_steel_plant_tracker", color = "0.10 0.6 0.85", style="rounded"];
            74[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.34 0.6 0.85", style="rounded"];
            75[label = "build_industrial_energy_demand_per_node_today", color = "0.28 0.6 0.85", style="rounded"];
            76[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.57 0.6 0.85", style="rounded"];
            77[label = "build_cop_profiles", color = "0.02 0.6 0.85", style="rounded"];
            78[label = "build_central_heating_temperature_profiles", color = "0.15 0.6 0.85", style="rounded"];
            79[label = "plot_power_network_clustered", color = "0.43 0.6 0.85", style="rounded"];
            80[label = "plot_power_network", color = "0.05 0.6 0.85", style="rounded"];
            81[label = "plot_hydrogen_network", color = "0.52 0.6 0.85", style="rounded"];
            82[label = "plot_gas_network", color = "0.46 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            49 -> 1
            7 -> 1
            3 -> 2
            42 -> 2
            79 -> 2
            80 -> 2
            81 -> 2
            82 -> 2
            4 -> 3
            5 -> 4
            24 -> 4
            26 -> 4
            28 -> 4
            31 -> 4
            33 -> 4
            34 -> 4
            49 -> 4
            50 -> 4
            55 -> 4
            57 -> 4
            58 -> 4
            51 -> 4
            7 -> 4
            60 -> 4
            42 -> 4
            62 -> 4
            15 -> 4
            14 -> 4
            54 -> 4
            63 -> 4
            44 -> 4
            70 -> 4
            76 -> 4
            59 -> 4
            77 -> 4
            48 -> 4
            6 -> 5
            11 -> 5
            14 -> 5
            10 -> 5
            7 -> 6
            8 -> 6
            11 -> 6
            14 -> 6
            10 -> 6
            9 -> 8
            10 -> 8
            12 -> 11
            13 -> 11
            7 -> 11
            15 -> 14
            20 -> 14
            16 -> 15
            17 -> 15
            17 -> 16
            19 -> 16
            18 -> 17
            11 -> 17
            17 -> 19
            11 -> 19
            15 -> 20
            11 -> 20
            21 -> 20
            22 -> 21
            23 -> 21
            25 -> 24
            11 -> 24
            14 -> 24
            10 -> 24
            7 -> 25
            8 -> 25
            11 -> 25
            14 -> 25
            10 -> 25
            27 -> 26
            11 -> 26
            14 -> 26
            10 -> 26
            7 -> 27
            8 -> 27
            11 -> 27
            14 -> 27
            10 -> 27
            29 -> 28
            14 -> 28
            30 -> 29
            32 -> 31
            30 -> 31
            14 -> 31
            34 -> 33
            44 -> 33
            48 -> 33
            35 -> 34
            42 -> 34
            36 -> 35
            38 -> 35
            40 -> 35
            5 -> 35
            24 -> 35
            26 -> 35
            14 -> 35
            42 -> 35
            43 -> 35
            20 -> 35
            37 -> 36
            11 -> 36
            14 -> 36
            10 -> 36
            7 -> 37
            11 -> 37
            14 -> 37
            10 -> 37
            39 -> 38
            11 -> 38
            14 -> 38
            10 -> 38
            7 -> 39
            11 -> 39
            14 -> 39
            10 -> 39
            41 -> 40
            11 -> 40
            14 -> 40
            10 -> 40
            7 -> 41
            11 -> 41
            14 -> 41
            10 -> 41
            14 -> 43
            45 -> 44
            46 -> 45
            14 -> 45
            10 -> 45
            11 -> 46
            47 -> 46
            10 -> 46
            46 -> 48
            14 -> 48
            10 -> 48
            51 -> 50
            54 -> 50
            11 -> 51
            7 -> 51
            52 -> 51
            49 -> 51
            53 -> 51
            46 -> 54
            14 -> 54
            10 -> 54
            56 -> 55
            54 -> 55
            51 -> 56
            11 -> 57
            14 -> 57
            51 -> 57
            54 -> 58
            50 -> 58
            51 -> 58
            7 -> 58
            59 -> 58
            46 -> 59
            14 -> 59
            10 -> 59
            61 -> 60
            49 -> 60
            7 -> 60
            14 -> 60
            11 -> 60
            7 -> 62
            14 -> 62
            64 -> 63
            70 -> 63
            75 -> 63
            65 -> 64
            68 -> 64
            69 -> 64
            66 -> 65
            52 -> 65
            67 -> 66
            51 -> 68
            52 -> 68
            69 -> 68
            66 -> 69
            52 -> 69
            49 -> 69
            71 -> 70
            74 -> 70
            14 -> 71
            54 -> 71
            72 -> 71
            73 -> 71
            69 -> 74
            71 -> 75
            68 -> 75
            51 -> 76
            54 -> 76
            78 -> 77
            59 -> 77
            14 -> 77
            59 -> 78
            14 -> 78
            14 -> 79
            3 -> 80
            14 -> 80
            3 -> 81
            14 -> 81
            3 -> 82
            14 -> 82
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

.. code:: console

    $ snakemake all --configfile config/test/config.myopic.yaml

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
            0[label = "all", color = "0.58 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.14 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.16 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network_myopic", color = "0.04 0.6 0.85", style="rounded"];
            4[label = "add_existing_baseyear", color = "0.29 0.6 0.85", style="rounded"];
            5[label = "prepare_sector_network\nsector_opts: ", color = "0.01 0.6 0.85", style="rounded"];
            6[label = "build_renewable_profiles", color = "0.29 0.6 0.85", style="rounded"];
            7[label = "determine_availability_matrix\ntechnology: offwind-ac", color = "0.48 0.6 0.85", style="rounded"];
            8[label = "retrieve_databundle", color = "0.25 0.6 0.85", style="rounded"];
            9[label = "build_ship_raster", color = "0.35 0.6 0.85", style="rounded"];
            10[label = "retrieve_ship_raster", color = "0.36 0.6 0.85", style="rounded"];
            11[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.37 0.6 0.85", style="rounded"];
            12[label = "build_shapes", color = "0.64 0.6 0.85", style="rounded"];
            13[label = "retrieve_naturalearth_countries", color = "0.39 0.6 0.85", style="rounded"];
            14[label = "retrieve_eez", color = "0.43 0.6 0.85", style="rounded"];
            15[label = "cluster_network\nclusters: 5", color = "0.36 0.6 0.85", style="rounded"];
            16[label = "simplify_network", color = "0.13 0.6 0.85", style="rounded"];
            17[label = "add_transmission_projects_and_dlr", color = "0.05 0.6 0.85", style="rounded"];
            18[label = "base_network", color = "0.34 0.6 0.85", style="rounded"];
            19[label = "retrieve_osm_prebuilt", color = "0.39 0.6 0.85", style="rounded"];
            20[label = "build_transmission_projects", color = "0.17 0.6 0.85", style="rounded"];
            21[label = "build_electricity_demand_base", color = "0.41 0.6 0.85", style="rounded"];
            22[label = "build_electricity_demand", color = "0.26 0.6 0.85", style="rounded"];
            23[label = "retrieve_electricity_demand", color = "0.32 0.6 0.85", style="rounded"];
            24[label = "retrieve_synthetic_electricity_demand", color = "0.60 0.6 0.85", style="rounded"];
            25[label = "build_renewable_profiles", color = "0.29 0.6 0.85", style="rounded"];
            26[label = "determine_availability_matrix\ntechnology: offwind-dc", color = "0.48 0.6 0.85", style="rounded"];
            27[label = "build_renewable_profiles", color = "0.29 0.6 0.85", style="rounded"];
            28[label = "determine_availability_matrix\ntechnology: offwind-float", color = "0.48 0.6 0.85", style="rounded"];
            29[label = "cluster_gas_network", color = "0.50 0.6 0.85", style="rounded"];
            30[label = "build_gas_network", color = "0.12 0.6 0.85", style="rounded"];
            31[label = "retrieve_gas_infrastructure_data", color = "0.09 0.6 0.85", style="rounded"];
            32[label = "build_gas_input_locations", color = "0.06 0.6 0.85", style="rounded"];
            33[label = "retrieve_gem_europe_gas_tracker", color = "0.11 0.6 0.85", style="rounded"];
            34[label = "time_aggregation", color = "0.64 0.6 0.85", style="rounded"];
            35[label = "prepare_network\nll: v1.5\nopts: ", color = "0.25 0.6 0.85", style="rounded"];
            36[label = "add_electricity", color = "0.30 0.6 0.85", style="rounded"];
            37[label = "build_renewable_profiles", color = "0.29 0.6 0.85", style="rounded"];
            38[label = "determine_availability_matrix\ntechnology: solar", color = "0.48 0.6 0.85", style="rounded"];
            39[label = "build_renewable_profiles", color = "0.29 0.6 0.85", style="rounded"];
            40[label = "determine_availability_matrix\ntechnology: solar-hsat", color = "0.48 0.6 0.85", style="rounded"];
            41[label = "build_renewable_profiles", color = "0.29 0.6 0.85", style="rounded"];
            42[label = "determine_availability_matrix\ntechnology: onwind", color = "0.48 0.6 0.85", style="rounded"];
            43[label = "retrieve_cost_data\nyear: 2030", color = "0.61 0.6 0.85", style="rounded"];
            44[label = "build_powerplants", color = "0.51 0.6 0.85", style="rounded"];
            45[label = "build_hourly_heat_demand", color = "0.07 0.6 0.85", style="rounded"];
            46[label = "build_daily_heat_demand", color = "0.12 0.6 0.85", style="rounded"];
            47[label = "build_population_layouts", color = "0.40 0.6 0.85", style="rounded"];
            48[label = "retrieve_worldbank_urban_population", color = "0.65 0.6 0.85", style="rounded"];
            49[label = "build_solar_thermal_profiles", color = "0.40 0.6 0.85", style="rounded"];
            50[label = "retrieve_eurostat_data", color = "0.45 0.6 0.85", style="rounded"];
            51[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.02 0.6 0.85", style="rounded"];
            52[label = "build_energy_totals", color = "0.23 0.6 0.85", style="rounded"];
            53[label = "retrieve_jrc_idees", color = "0.35 0.6 0.85", style="rounded"];
            54[label = "retrieve_eurostat_household_data", color = "0.19 0.6 0.85", style="rounded"];
            55[label = "build_clustered_population_layouts", color = "0.24 0.6 0.85", style="rounded"];
            56[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.02 0.6 0.85", style="rounded"];
            57[label = "build_heat_totals", color = "0.66 0.6 0.85", style="rounded"];
            58[label = "build_shipping_demand", color = "0.59 0.6 0.85", style="rounded"];
            59[label = "build_transport_demand", color = "0.19 0.6 0.85", style="rounded"];
            60[label = "build_temperature_profiles", color = "0.27 0.6 0.85", style="rounded"];
            61[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.08 0.6 0.85", style="rounded"];
            62[label = "retrieve_jrc_enspreso_biomass", color = "0.18 0.6 0.85", style="rounded"];
            63[label = "build_salt_cavern_potentials", color = "0.57 0.6 0.85", style="rounded"];
            64[label = "build_industrial_energy_demand_per_node", color = "0.13 0.6 0.85", style="rounded"];
            65[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.05 0.6 0.85", style="rounded"];
            66[label = "build_industry_sector_ratios", color = "0.28 0.6 0.85", style="rounded"];
            67[label = "build_ammonia_production", color = "0.22 0.6 0.85", style="rounded"];
            68[label = "retrieve_usgs_ammonia_production", color = "0.49 0.6 0.85", style="rounded"];
            69[label = "build_industrial_energy_demand_per_country_today", color = "0.20 0.6 0.85", style="rounded"];
            70[label = "build_industrial_production_per_country", color = "0.18 0.6 0.85", style="rounded"];
            71[label = "build_industrial_production_per_node", color = "0.32 0.6 0.85", style="rounded"];
            72[label = "build_industrial_distribution_key", color = "0.55 0.6 0.85", style="rounded"];
            73[label = "retrieve_hotmaps_industrial_sites", color = "0.16 0.6 0.85", style="rounded"];
            74[label = "retrieve_gem_steel_plant_tracker", color = "0.47 0.6 0.85", style="rounded"];
            75[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.21 0.6 0.85", style="rounded"];
            76[label = "build_industrial_energy_demand_per_node_today", color = "0.00 0.6 0.85", style="rounded"];
            77[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.08 0.6 0.85", style="rounded"];
            78[label = "build_cop_profiles", color = "0.44 0.6 0.85", style="rounded"];
            79[label = "build_central_heating_temperature_profiles", color = "0.42 0.6 0.85", style="rounded"];
            80[label = "build_existing_heating_distribution", color = "0.42 0.6 0.85", style="rounded"];
            81[label = "solve_sector_network_myopic", color = "0.04 0.6 0.85", style="rounded"];
            82[label = "add_brownfield", color = "0.37 0.6 0.85", style="rounded"];
            83[label = "prepare_sector_network\nsector_opts: ", color = "0.01 0.6 0.85", style="rounded"];
            84[label = "build_biomass_potentials\nplanning_horizons: 2040", color = "0.08 0.6 0.85", style="rounded"];
            85[label = "retrieve_cost_data\nyear: 2040", color = "0.61 0.6 0.85", style="rounded"];
            86[label = "build_industrial_energy_demand_per_node", color = "0.13 0.6 0.85", style="rounded"];
            87[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2040", color = "0.05 0.6 0.85", style="rounded"];
            88[label = "build_industrial_production_per_node", color = "0.32 0.6 0.85", style="rounded"];
            89[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2040", color = "0.21 0.6 0.85", style="rounded"];
            90[label = "build_district_heat_share\nplanning_horizons: 2040", color = "0.08 0.6 0.85", style="rounded"];
            91[label = "solve_sector_network_myopic", color = "0.04 0.6 0.85", style="rounded"];
            92[label = "add_brownfield", color = "0.37 0.6 0.85", style="rounded"];
            93[label = "prepare_sector_network\nsector_opts: ", color = "0.01 0.6 0.85", style="rounded"];
            94[label = "build_biomass_potentials\nplanning_horizons: 2050", color = "0.08 0.6 0.85", style="rounded"];
            95[label = "retrieve_cost_data\nyear: 2050", color = "0.61 0.6 0.85", style="rounded"];
            96[label = "build_industrial_energy_demand_per_node", color = "0.13 0.6 0.85", style="rounded"];
            97[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2050", color = "0.05 0.6 0.85", style="rounded"];
            98[label = "build_industrial_production_per_node", color = "0.32 0.6 0.85", style="rounded"];
            99[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2050", color = "0.21 0.6 0.85", style="rounded"];
            100[label = "build_district_heat_share\nplanning_horizons: 2050", color = "0.08 0.6 0.85", style="rounded"];
            101[label = "plot_power_network_clustered", color = "0.27 0.6 0.85", style="rounded"];
            102[label = "plot_power_network", color = "0.54 0.6 0.85", style="rounded"];
            103[label = "plot_power_network", color = "0.54 0.6 0.85", style="rounded"];
            104[label = "plot_power_network", color = "0.54 0.6 0.85", style="rounded"];
            105[label = "plot_hydrogen_network", color = "0.02 0.6 0.85", style="rounded"];
            106[label = "plot_hydrogen_network", color = "0.02 0.6 0.85", style="rounded"];
            107[label = "plot_hydrogen_network", color = "0.02 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            50 -> 1
            8 -> 1
            3 -> 2
            81 -> 2
            91 -> 2
            43 -> 2
            101 -> 2
            102 -> 2
            103 -> 2
            104 -> 2
            105 -> 2
            106 -> 2
            107 -> 2
            4 -> 3
            43 -> 3
            5 -> 4
            44 -> 4
            16 -> 4
            15 -> 4
            55 -> 4
            43 -> 4
            78 -> 4
            80 -> 4
            52 -> 4
            6 -> 5
            25 -> 5
            27 -> 5
            29 -> 5
            32 -> 5
            34 -> 5
            35 -> 5
            50 -> 5
            51 -> 5
            56 -> 5
            58 -> 5
            59 -> 5
            52 -> 5
            8 -> 5
            61 -> 5
            43 -> 5
            63 -> 5
            16 -> 5
            15 -> 5
            55 -> 5
            64 -> 5
            45 -> 5
            71 -> 5
            77 -> 5
            60 -> 5
            78 -> 5
            49 -> 5
            7 -> 6
            12 -> 6
            15 -> 6
            11 -> 6
            8 -> 7
            9 -> 7
            12 -> 7
            15 -> 7
            11 -> 7
            10 -> 9
            11 -> 9
            13 -> 12
            14 -> 12
            8 -> 12
            16 -> 15
            21 -> 15
            17 -> 16
            18 -> 16
            18 -> 17
            20 -> 17
            19 -> 18
            12 -> 18
            18 -> 20
            12 -> 20
            16 -> 21
            12 -> 21
            22 -> 21
            23 -> 22
            24 -> 22
            26 -> 25
            12 -> 25
            15 -> 25
            11 -> 25
            8 -> 26
            9 -> 26
            12 -> 26
            15 -> 26
            11 -> 26
            28 -> 27
            12 -> 27
            15 -> 27
            11 -> 27
            8 -> 28
            9 -> 28
            12 -> 28
            15 -> 28
            11 -> 28
            30 -> 29
            15 -> 29
            31 -> 30
            33 -> 32
            31 -> 32
            15 -> 32
            35 -> 34
            45 -> 34
            49 -> 34
            36 -> 35
            43 -> 35
            37 -> 36
            39 -> 36
            41 -> 36
            6 -> 36
            25 -> 36
            27 -> 36
            15 -> 36
            43 -> 36
            44 -> 36
            21 -> 36
            38 -> 37
            12 -> 37
            15 -> 37
            11 -> 37
            8 -> 38
            12 -> 38
            15 -> 38
            11 -> 38
            40 -> 39
            12 -> 39
            15 -> 39
            11 -> 39
            8 -> 40
            12 -> 40
            15 -> 40
            11 -> 40
            42 -> 41
            12 -> 41
            15 -> 41
            11 -> 41
            8 -> 42
            12 -> 42
            15 -> 42
            11 -> 42
            15 -> 44
            46 -> 45
            47 -> 46
            15 -> 46
            11 -> 46
            12 -> 47
            48 -> 47
            11 -> 47
            47 -> 49
            15 -> 49
            11 -> 49
            52 -> 51
            55 -> 51
            12 -> 52
            8 -> 52
            53 -> 52
            50 -> 52
            54 -> 52
            47 -> 55
            15 -> 55
            11 -> 55
            57 -> 56
            55 -> 56
            52 -> 57
            12 -> 58
            15 -> 58
            52 -> 58
            55 -> 59
            51 -> 59
            52 -> 59
            8 -> 59
            60 -> 59
            47 -> 60
            15 -> 60
            11 -> 60
            62 -> 61
            50 -> 61
            8 -> 61
            15 -> 61
            12 -> 61
            8 -> 63
            15 -> 63
            65 -> 64
            71 -> 64
            76 -> 64
            66 -> 65
            69 -> 65
            70 -> 65
            67 -> 66
            53 -> 66
            68 -> 67
            52 -> 69
            53 -> 69
            70 -> 69
            67 -> 70
            53 -> 70
            50 -> 70
            72 -> 71
            75 -> 71
            15 -> 72
            55 -> 72
            73 -> 72
            74 -> 72
            70 -> 75
            72 -> 76
            69 -> 76
            52 -> 77
            55 -> 77
            79 -> 78
            60 -> 78
            15 -> 78
            60 -> 79
            15 -> 79
            55 -> 80
            51 -> 80
            77 -> 80
            82 -> 81
            85 -> 81
            37 -> 82
            39 -> 82
            41 -> 82
            6 -> 82
            25 -> 82
            27 -> 82
            16 -> 82
            15 -> 82
            83 -> 82
            3 -> 82
            85 -> 82
            78 -> 82
            6 -> 83
            25 -> 83
            27 -> 83
            29 -> 83
            32 -> 83
            34 -> 83
            35 -> 83
            50 -> 83
            51 -> 83
            56 -> 83
            58 -> 83
            59 -> 83
            52 -> 83
            8 -> 83
            84 -> 83
            85 -> 83
            63 -> 83
            16 -> 83
            15 -> 83
            55 -> 83
            86 -> 83
            45 -> 83
            88 -> 83
            90 -> 83
            60 -> 83
            78 -> 83
            49 -> 83
            62 -> 84
            50 -> 84
            8 -> 84
            15 -> 84
            12 -> 84
            87 -> 86
            88 -> 86
            76 -> 86
            66 -> 87
            69 -> 87
            70 -> 87
            72 -> 88
            89 -> 88
            70 -> 89
            52 -> 90
            55 -> 90
            92 -> 91
            95 -> 91
            37 -> 92
            39 -> 92
            41 -> 92
            6 -> 92
            25 -> 92
            27 -> 92
            16 -> 92
            15 -> 92
            93 -> 92
            81 -> 92
            95 -> 92
            78 -> 92
            6 -> 93
            25 -> 93
            27 -> 93
            29 -> 93
            32 -> 93
            34 -> 93
            35 -> 93
            50 -> 93
            51 -> 93
            56 -> 93
            58 -> 93
            59 -> 93
            52 -> 93
            8 -> 93
            94 -> 93
            95 -> 93
            63 -> 93
            16 -> 93
            15 -> 93
            55 -> 93
            96 -> 93
            45 -> 93
            98 -> 93
            100 -> 93
            60 -> 93
            78 -> 93
            49 -> 93
            62 -> 94
            50 -> 94
            8 -> 94
            15 -> 94
            12 -> 94
            97 -> 96
            98 -> 96
            76 -> 96
            66 -> 97
            69 -> 97
            70 -> 97
            72 -> 98
            99 -> 98
            70 -> 99
            52 -> 100
            55 -> 100
            15 -> 101
            3 -> 102
            15 -> 102
            81 -> 103
            15 -> 103
            91 -> 104
            15 -> 104
            3 -> 105
            15 -> 105
            81 -> 106
            15 -> 106
            91 -> 107
            15 -> 107
    }

|


Scaling-Up
==========

If you now feel confident and want to tackle runs with larger temporal, technological and
spatial scope, clean-up the repository and after modifying the ``config/config.yaml`` file
target the collection rule ``all`` again without providing the test
configuration file.

.. code:: console

    $ snakemake purge
    $ snakemake all

.. note::

    It is good practice to perform a dry-run using the option `-n`, before you
    commit to a run:

    .. code:: console

        $ snakemake all -n
