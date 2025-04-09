..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

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
   :end-before: industry:

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
    build_direct_heat_source_utilisation_profiles           1
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
    build_osm_boundaries                                    4
    build_population_layouts                                1
    build_population_weighted_energy_totals                 2
    build_powerplants                                       1
    build_renewable_profiles                                6
    build_salt_cavern_potentials                            1
    build_shapes                                            1
    build_ship_raster                                       1
    build_shipping_demand                                   1
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
    retrieve_databundle                                     1
    retrieve_eez                                            1
    retrieve_electricity_demand                             1
    retrieve_eurostat_data                                  1
    retrieve_eurostat_household_data                        1
    retrieve_gas_infrastructure_data                        1
    retrieve_gem_europe_gas_tracker                         1
    retrieve_gem_steel_plant_tracker                        1
    retrieve_hotmaps_industrial_sites                       1
    retrieve_jrc_ardeco                                     1
    retrieve_jrc_enspreso_biomass                           1
    retrieve_jrc_idees                                      1
    retrieve_nuts_2013_shapes                               1
    retrieve_nuts_2021_shapes                               1
    retrieve_osm_boundaries                                 4
    retrieve_osm_prebuilt                                   1
    retrieve_ship_raster                                    1
    retrieve_synthetic_electricity_demand                   1
    retrieve_usgs_ammonia_production                        1
    retrieve_worldbank_urban_population                     1
    simplify_network                                        1
    solve_sector_network                                    1
    time_aggregation                                        1
    total                                                  92

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
            0[label = "all", color = "0.25 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.29 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.30 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network", color = "0.34 0.6 0.85", style="rounded"];
            4[label = "prepare_sector_network", color = "0.32 0.6 0.85", style="rounded"];
            5[label = "build_renewable_profiles", color = "0.43 0.6 0.85", style="rounded"];
            6[label = "determine_availability_matrix\ntechnology: offwind-ac", color = "0.09 0.6 0.85", style="rounded"];
            7[label = "retrieve_databundle", color = "0.45 0.6 0.85", style="rounded"];
            8[label = "build_ship_raster", color = "0.11 0.6 0.85", style="rounded"];
            9[label = "retrieve_ship_raster", color = "0.28 0.6 0.85", style="rounded"];
            10[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.65 0.6 0.85", style="rounded"];
            11[label = "build_shapes", color = "0.31 0.6 0.85", style="rounded"];
            12[label = "retrieve_eez", color = "0.02 0.6 0.85", style="rounded"];
            13[label = "retrieve_nuts_2021_shapes", color = "0.38 0.6 0.85", style="rounded"];
            14[label = "build_osm_boundaries", color = "0.08 0.6 0.85", style="rounded"];
            15[label = "retrieve_osm_boundaries\ncountry: BA", color = "0.35 0.6 0.85", style="rounded"];
            16[label = "build_osm_boundaries", color = "0.08 0.6 0.85", style="rounded"];
            17[label = "retrieve_osm_boundaries\ncountry: MD", color = "0.35 0.6 0.85", style="rounded"];
            18[label = "build_osm_boundaries", color = "0.08 0.6 0.85", style="rounded"];
            19[label = "retrieve_osm_boundaries\ncountry: UA", color = "0.35 0.6 0.85", style="rounded"];
            20[label = "build_osm_boundaries", color = "0.08 0.6 0.85", style="rounded"];
            21[label = "retrieve_osm_boundaries\ncountry: XK", color = "0.35 0.6 0.85", style="rounded"];
            22[label = "retrieve_jrc_ardeco", color = "0.22 0.6 0.85", style="rounded"];
            23[label = "cluster_network\nclusters: 5", color = "0.39 0.6 0.85", style="rounded"];
            24[label = "simplify_network", color = "0.13 0.6 0.85", style="rounded"];
            25[label = "add_transmission_projects_and_dlr", color = "0.44 0.6 0.85", style="rounded"];
            26[label = "base_network", color = "0.62 0.6 0.85", style="rounded"];
            27[label = "retrieve_osm_prebuilt", color = "0.52 0.6 0.85", style="rounded"];
            28[label = "build_transmission_projects", color = "0.32 0.6 0.85", style="rounded"];
            29[label = "build_electricity_demand_base", color = "0.47 0.6 0.85", style="rounded"];
            30[label = "build_electricity_demand", color = "0.34 0.6 0.85", style="rounded"];
            31[label = "retrieve_electricity_demand", color = "0.61 0.6 0.85", style="rounded"];
            32[label = "retrieve_synthetic_electricity_demand", color = "0.56 0.6 0.85", style="rounded"];
            33[label = "build_renewable_profiles", color = "0.43 0.6 0.85", style="rounded"];
            34[label = "determine_availability_matrix\ntechnology: offwind-dc", color = "0.09 0.6 0.85", style="rounded"];
            35[label = "build_renewable_profiles", color = "0.43 0.6 0.85", style="rounded"];
            36[label = "determine_availability_matrix\ntechnology: offwind-float", color = "0.09 0.6 0.85", style="rounded"];
            37[label = "cluster_gas_network", color = "0.03 0.6 0.85", style="rounded"];
            38[label = "build_gas_network", color = "0.50 0.6 0.85", style="rounded"];
            39[label = "retrieve_gas_infrastructure_data", color = "0.49 0.6 0.85", style="rounded"];
            40[label = "build_gas_input_locations", color = "0.23 0.6 0.85", style="rounded"];
            41[label = "retrieve_gem_europe_gas_tracker", color = "0.33 0.6 0.85", style="rounded"];
            42[label = "time_aggregation\nsector_opts: ", color = "0.06 0.6 0.85", style="rounded"];
            43[label = "prepare_network\nopts: ", color = "0.00 0.6 0.85", style="rounded"];
            44[label = "add_electricity", color = "0.15 0.6 0.85", style="rounded"];
            45[label = "build_renewable_profiles", color = "0.43 0.6 0.85", style="rounded"];
            46[label = "determine_availability_matrix\ntechnology: solar", color = "0.09 0.6 0.85", style="rounded"];
            47[label = "build_renewable_profiles", color = "0.43 0.6 0.85", style="rounded"];
            48[label = "determine_availability_matrix\ntechnology: solar-hsat", color = "0.09 0.6 0.85", style="rounded"];
            49[label = "build_renewable_profiles", color = "0.43 0.6 0.85", style="rounded"];
            50[label = "determine_availability_matrix\ntechnology: onwind", color = "0.09 0.6 0.85", style="rounded"];
            51[label = "retrieve_cost_data\nyear: 2040", color = "0.10 0.6 0.85", style="rounded"];
            52[label = "build_powerplants", color = "0.42 0.6 0.85", style="rounded"];
            53[label = "build_hourly_heat_demand", color = "0.41 0.6 0.85", style="rounded"];
            54[label = "build_daily_heat_demand", color = "0.66 0.6 0.85", style="rounded"];
            55[label = "build_population_layouts", color = "0.16 0.6 0.85", style="rounded"];
            56[label = "retrieve_worldbank_urban_population", color = "0.12 0.6 0.85", style="rounded"];
            57[label = "retrieve_eurostat_data", color = "0.38 0.6 0.85", style="rounded"];
            58[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.37 0.6 0.85", style="rounded"];
            59[label = "build_energy_totals", color = "0.43 0.6 0.85", style="rounded"];
            60[label = "retrieve_jrc_idees", color = "0.36 0.6 0.85", style="rounded"];
            61[label = "retrieve_eurostat_household_data", color = "0.46 0.6 0.85", style="rounded"];
            62[label = "build_clustered_population_layouts", color = "0.05 0.6 0.85", style="rounded"];
            63[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.37 0.6 0.85", style="rounded"];
            64[label = "build_heat_totals", color = "0.04 0.6 0.85", style="rounded"];
            65[label = "build_shipping_demand", color = "0.16 0.6 0.85", style="rounded"];
            66[label = "build_transport_demand", color = "0.03 0.6 0.85", style="rounded"];
            67[label = "build_temperature_profiles", color = "0.10 0.6 0.85", style="rounded"];
            68[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.20 0.6 0.85", style="rounded"];
            69[label = "retrieve_jrc_enspreso_biomass", color = "0.49 0.6 0.85", style="rounded"];
            70[label = "retrieve_nuts_2013_shapes", color = "0.57 0.6 0.85", style="rounded"];
            71[label = "build_salt_cavern_potentials", color = "0.64 0.6 0.85", style="rounded"];
            72[label = "build_industrial_energy_demand_per_node", color = "0.26 0.6 0.85", style="rounded"];
            73[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.24 0.6 0.85", style="rounded"];
            74[label = "build_industry_sector_ratios", color = "0.27 0.6 0.85", style="rounded"];
            75[label = "build_ammonia_production", color = "0.54 0.6 0.85", style="rounded"];
            76[label = "retrieve_usgs_ammonia_production", color = "0.53 0.6 0.85", style="rounded"];
            77[label = "build_industrial_energy_demand_per_country_today", color = "0.17 0.6 0.85", style="rounded"];
            78[label = "build_industrial_production_per_country", color = "0.26 0.6 0.85", style="rounded"];
            79[label = "build_industrial_production_per_node", color = "0.58 0.6 0.85", style="rounded"];
            80[label = "build_industrial_distribution_key", color = "0.11 0.6 0.85", style="rounded"];
            81[label = "retrieve_hotmaps_industrial_sites", color = "0.21 0.6 0.85", style="rounded"];
            82[label = "retrieve_gem_steel_plant_tracker", color = "0.12 0.6 0.85", style="rounded"];
            83[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.54 0.6 0.85", style="rounded"];
            84[label = "build_industrial_energy_demand_per_node_today", color = "0.55 0.6 0.85", style="rounded"];
            85[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.31 0.6 0.85", style="rounded"];
            86[label = "build_cop_profiles", color = "0.19 0.6 0.85", style="rounded"];
            87[label = "build_central_heating_temperature_profiles\nplanning_horizons: 2030", color = "0.17 0.6 0.85", style="rounded"];
            88[label = "build_direct_heat_source_utilisation_profiles", color = "0.41 0.6 0.85", style="rounded"];
            89[label = "plot_power_network_clustered", color = "0.58 0.6 0.85", style="rounded"];
            90[label = "plot_power_network", color = "0.13 0.6 0.85", style="rounded"];
            91[label = "plot_hydrogen_network", color = "0.24 0.6 0.85", style="rounded"];
            92[label = "plot_gas_network", color = "0.04 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            57 -> 1
            7 -> 1
            3 -> 2
            51 -> 2
            89 -> 2
            90 -> 2
            91 -> 2
            92 -> 2
            4 -> 3
            5 -> 4
            33 -> 4
            35 -> 4
            37 -> 4
            40 -> 4
            42 -> 4
            43 -> 4
            57 -> 4
            58 -> 4
            63 -> 4
            65 -> 4
            66 -> 4
            59 -> 4
            7 -> 4
            68 -> 4
            51 -> 4
            71 -> 4
            24 -> 4
            23 -> 4
            62 -> 4
            72 -> 4
            53 -> 4
            79 -> 4
            85 -> 4
            67 -> 4
            86 -> 4
            88 -> 4
            6 -> 5
            11 -> 5
            23 -> 5
            10 -> 5
            7 -> 6
            8 -> 6
            11 -> 6
            23 -> 6
            10 -> 6
            9 -> 8
            10 -> 8
            12 -> 11
            13 -> 11
            14 -> 11
            16 -> 11
            18 -> 11
            20 -> 11
            22 -> 11
            7 -> 11
            15 -> 14
            12 -> 14
            17 -> 16
            12 -> 16
            19 -> 18
            12 -> 18
            21 -> 20
            12 -> 20
            24 -> 23
            29 -> 23
            25 -> 24
            26 -> 24
            26 -> 25
            28 -> 25
            27 -> 26
            11 -> 26
            26 -> 28
            11 -> 28
            24 -> 29
            11 -> 29
            30 -> 29
            31 -> 30
            32 -> 30
            34 -> 33
            11 -> 33
            23 -> 33
            10 -> 33
            7 -> 34
            8 -> 34
            11 -> 34
            23 -> 34
            10 -> 34
            36 -> 35
            11 -> 35
            23 -> 35
            10 -> 35
            7 -> 36
            8 -> 36
            11 -> 36
            23 -> 36
            10 -> 36
            38 -> 37
            23 -> 37
            39 -> 38
            41 -> 40
            39 -> 40
            23 -> 40
            43 -> 42
            53 -> 42
            44 -> 43
            51 -> 43
            45 -> 44
            47 -> 44
            49 -> 44
            5 -> 44
            33 -> 44
            35 -> 44
            23 -> 44
            51 -> 44
            52 -> 44
            29 -> 44
            46 -> 45
            11 -> 45
            23 -> 45
            10 -> 45
            7 -> 46
            11 -> 46
            23 -> 46
            10 -> 46
            48 -> 47
            11 -> 47
            23 -> 47
            10 -> 47
            7 -> 48
            11 -> 48
            23 -> 48
            10 -> 48
            50 -> 49
            11 -> 49
            23 -> 49
            10 -> 49
            7 -> 50
            11 -> 50
            23 -> 50
            10 -> 50
            23 -> 52
            54 -> 53
            55 -> 54
            23 -> 54
            10 -> 54
            11 -> 55
            56 -> 55
            10 -> 55
            59 -> 58
            62 -> 58
            11 -> 59
            7 -> 59
            60 -> 59
            57 -> 59
            61 -> 59
            55 -> 62
            23 -> 62
            10 -> 62
            64 -> 63
            62 -> 63
            59 -> 64
            11 -> 65
            23 -> 65
            59 -> 65
            62 -> 66
            58 -> 66
            59 -> 66
            7 -> 66
            67 -> 66
            55 -> 67
            23 -> 67
            10 -> 67
            69 -> 68
            57 -> 68
            70 -> 68
            23 -> 68
            7 -> 68
            11 -> 68
            7 -> 71
            23 -> 71
            73 -> 72
            79 -> 72
            84 -> 72
            74 -> 73
            77 -> 73
            78 -> 73
            75 -> 74
            60 -> 74
            76 -> 75
            59 -> 77
            60 -> 77
            78 -> 77
            75 -> 78
            60 -> 78
            57 -> 78
            80 -> 79
            83 -> 79
            23 -> 80
            62 -> 80
            81 -> 80
            82 -> 80
            78 -> 83
            80 -> 84
            77 -> 84
            59 -> 85
            62 -> 85
            87 -> 86
            67 -> 86
            23 -> 86
            67 -> 87
            23 -> 87
            87 -> 88
            23 -> 89
            3 -> 90
            23 -> 90
            3 -> 91
            23 -> 91
            3 -> 92
            23 -> 92
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
            0[label = "all", color = "0.50 0.6 0.85", style="rounded"];
            1[label = "plot_summary", color = "0.26 0.6 0.85", style="rounded"];
            2[label = "make_summary", color = "0.51 0.6 0.85", style="rounded"];
            3[label = "solve_sector_network_myopic", color = "0.09 0.6 0.85", style="rounded"];
            4[label = "add_existing_baseyear", color = "0.26 0.6 0.85", style="rounded"];
            5[label = "prepare_sector_network", color = "0.25 0.6 0.85", style="rounded"];
            6[label = "build_renewable_profiles", color = "0.53 0.6 0.85", style="rounded"];
            7[label = "determine_availability_matrix\ntechnology: offwind-ac", color = "0.58 0.6 0.85", style="rounded"];
            8[label = "retrieve_databundle", color = "0.28 0.6 0.85", style="rounded"];
            9[label = "build_ship_raster", color = "0.07 0.6 0.85", style="rounded"];
            10[label = "retrieve_ship_raster", color = "0.37 0.6 0.85", style="rounded"];
            11[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.55 0.6 0.85", style="rounded,dashed"];
            12[label = "build_shapes", color = "0.00 0.6 0.85", style="rounded"];
            13[label = "retrieve_eez", color = "0.66 0.6 0.85", style="rounded"];
            14[label = "retrieve_nuts_2021_shapes", color = "0.60 0.6 0.85", style="rounded"];
            15[label = "build_osm_boundaries", color = "0.21 0.6 0.85", style="rounded"];
            16[label = "retrieve_osm_boundaries\ncountry: BA", color = "0.23 0.6 0.85", style="rounded"];
            17[label = "build_osm_boundaries", color = "0.21 0.6 0.85", style="rounded"];
            18[label = "retrieve_osm_boundaries\ncountry: MD", color = "0.23 0.6 0.85", style="rounded"];
            19[label = "build_osm_boundaries", color = "0.21 0.6 0.85", style="rounded"];
            20[label = "retrieve_osm_boundaries\ncountry: UA", color = "0.23 0.6 0.85", style="rounded"];
            21[label = "build_osm_boundaries", color = "0.21 0.6 0.85", style="rounded"];
            22[label = "retrieve_osm_boundaries\ncountry: XK", color = "0.23 0.6 0.85", style="rounded"];
            23[label = "retrieve_jrc_ardeco", color = "0.29 0.6 0.85", style="rounded"];
            24[label = "cluster_network\nclusters: 5", color = "0.30 0.6 0.85", style="rounded"];
            25[label = "simplify_network", color = "0.02 0.6 0.85", style="rounded"];
            26[label = "add_transmission_projects_and_dlr", color = "0.60 0.6 0.85", style="rounded"];
            27[label = "base_network", color = "0.30 0.6 0.85", style="rounded"];
            28[label = "retrieve_osm_prebuilt", color = "0.43 0.6 0.85", style="rounded"];
            29[label = "build_transmission_projects", color = "0.51 0.6 0.85", style="rounded"];
            30[label = "build_electricity_demand_base", color = "0.64 0.6 0.85", style="rounded"];
            31[label = "build_electricity_demand", color = "0.16 0.6 0.85", style="rounded"];
            32[label = "retrieve_electricity_demand", color = "0.33 0.6 0.85", style="rounded"];
            33[label = "retrieve_synthetic_electricity_demand", color = "0.56 0.6 0.85", style="rounded"];
            34[label = "build_renewable_profiles", color = "0.53 0.6 0.85", style="rounded"];
            35[label = "determine_availability_matrix\ntechnology: offwind-dc", color = "0.58 0.6 0.85", style="rounded"];
            36[label = "build_renewable_profiles", color = "0.53 0.6 0.85", style="rounded"];
            37[label = "determine_availability_matrix\ntechnology: offwind-float", color = "0.58 0.6 0.85", style="rounded"];
            38[label = "cluster_gas_network", color = "0.13 0.6 0.85", style="rounded"];
            39[label = "build_gas_network", color = "0.22 0.6 0.85", style="rounded"];
            40[label = "retrieve_gas_infrastructure_data", color = "0.33 0.6 0.85", style="rounded"];
            41[label = "build_gas_input_locations", color = "0.12 0.6 0.85", style="rounded"];
            42[label = "retrieve_gem_europe_gas_tracker", color = "0.17 0.6 0.85", style="rounded"];
            43[label = "time_aggregation\nsector_opts: ", color = "0.40 0.6 0.85", style="rounded"];
            44[label = "prepare_network\nopts: ", color = "0.09 0.6 0.85", style="rounded"];
            45[label = "add_electricity", color = "0.03 0.6 0.85", style="rounded"];
            46[label = "build_renewable_profiles", color = "0.53 0.6 0.85", style="rounded"];
            47[label = "determine_availability_matrix\ntechnology: solar", color = "0.58 0.6 0.85", style="rounded"];
            48[label = "build_renewable_profiles", color = "0.53 0.6 0.85", style="rounded"];
            49[label = "determine_availability_matrix\ntechnology: solar-hsat", color = "0.58 0.6 0.85", style="rounded"];
            50[label = "build_renewable_profiles", color = "0.53 0.6 0.85", style="rounded"];
            51[label = "determine_availability_matrix\ntechnology: onwind", color = "0.58 0.6 0.85", style="rounded"];
            52[label = "retrieve_cost_data\nyear: 2040", color = "0.08 0.6 0.85", style="rounded"];
            53[label = "build_powerplants", color = "0.63 0.6 0.85", style="rounded"];
            54[label = "build_hourly_heat_demand", color = "0.27 0.6 0.85", style="rounded"];
            55[label = "build_daily_heat_demand", color = "0.02 0.6 0.85", style="rounded"];
            56[label = "build_population_layouts", color = "0.47 0.6 0.85", style="rounded"];
            57[label = "retrieve_worldbank_urban_population", color = "0.36 0.6 0.85", style="rounded"];
            58[label = "retrieve_eurostat_data", color = "0.39 0.6 0.85", style="rounded"];
            59[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.55 0.6 0.85", style="rounded"];
            60[label = "build_energy_totals", color = "0.57 0.6 0.85", style="rounded"];
            61[label = "retrieve_jrc_idees", color = "0.01 0.6 0.85", style="rounded"];
            62[label = "retrieve_eurostat_household_data", color = "0.28 0.6 0.85", style="rounded"];
            63[label = "build_clustered_population_layouts", color = "0.61 0.6 0.85", style="rounded"];
            64[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.55 0.6 0.85", style="rounded"];
            65[label = "build_heat_totals", color = "0.44 0.6 0.85", style="rounded"];
            66[label = "build_shipping_demand", color = "0.14 0.6 0.85", style="rounded"];
            67[label = "build_transport_demand", color = "0.34 0.6 0.85", style="rounded"];
            68[label = "build_temperature_profiles", color = "0.31 0.6 0.85", style="rounded"];
            69[label = "build_biomass_potentials\nplanning_horizons: 2030", color = "0.15 0.6 0.85", style="rounded"];
            70[label = "retrieve_jrc_enspreso_biomass", color = "0.52 0.6 0.85", style="rounded"];
            71[label = "retrieve_nuts_2013_shapes", color = "0.24 0.6 0.85", style="rounded"];
            72[label = "retrieve_cost_data\nyear: 2030", color = "0.08 0.6 0.85", style="rounded"];
            73[label = "build_salt_cavern_potentials", color = "0.05 0.6 0.85", style="rounded"];
            74[label = "build_industrial_energy_demand_per_node", color = "0.38 0.6 0.85", style="rounded"];
            75[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2030", color = "0.49 0.6 0.85", style="rounded"];
            76[label = "build_industry_sector_ratios", color = "0.39 0.6 0.85", style="rounded"];
            77[label = "build_ammonia_production", color = "0.11 0.6 0.85", style="rounded"];
            78[label = "retrieve_usgs_ammonia_production", color = "0.36 0.6 0.85", style="rounded"];
            79[label = "build_industrial_energy_demand_per_country_today", color = "0.41 0.6 0.85", style="rounded"];
            80[label = "build_industrial_production_per_country", color = "0.08 0.6 0.85", style="rounded"];
            81[label = "build_industrial_production_per_node", color = "0.12 0.6 0.85", style="rounded"];
            82[label = "build_industrial_distribution_key", color = "0.34 0.6 0.85", style="rounded"];
            83[label = "retrieve_hotmaps_industrial_sites", color = "0.24 0.6 0.85", style="rounded"];
            84[label = "retrieve_gem_steel_plant_tracker", color = "0.10 0.6 0.85", style="rounded"];
            85[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2030", color = "0.07 0.6 0.85", style="rounded"];
            86[label = "build_industrial_energy_demand_per_node_today", color = "0.37 0.6 0.85", style="rounded"];
            87[label = "build_district_heat_share\nplanning_horizons: 2030", color = "0.20 0.6 0.85", style="rounded"];
            88[label = "build_cop_profiles", color = "0.18 0.6 0.85", style="rounded"];
            89[label = "build_central_heating_temperature_profiles\nplanning_horizons: 2030", color = "0.06 0.6 0.85", style="rounded"];
            90[label = "build_direct_heat_source_utilisation_profiles", color = "0.52 0.6 0.85", style="rounded"];
            91[label = "build_existing_heating_distribution", color = "0.46 0.6 0.85", style="rounded"];
            92[label = "solve_sector_network_myopic", color = "0.09 0.6 0.85", style="rounded"];
            93[label = "add_brownfield", color = "0.54 0.6 0.85", style="rounded"];
            94[label = "prepare_sector_network", color = "0.25 0.6 0.85", style="rounded"];
            95[label = "build_biomass_potentials\nplanning_horizons: 2040", color = "0.15 0.6 0.85", style="rounded"];
            96[label = "build_industrial_energy_demand_per_node", color = "0.38 0.6 0.85", style="rounded"];
            97[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2040", color = "0.49 0.6 0.85", style="rounded"];
            98[label = "build_industrial_production_per_node", color = "0.12 0.6 0.85", style="rounded"];
            99[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2040", color = "0.07 0.6 0.85", style="rounded"];
            100[label = "build_district_heat_share\nplanning_horizons: 2040", color = "0.20 0.6 0.85", style="rounded"];
            101[label = "build_cop_profiles", color = "0.18 0.6 0.85", style="rounded"];
            102[label = "build_central_heating_temperature_profiles\nplanning_horizons: 2040", color = "0.06 0.6 0.85", style="rounded"];
            103[label = "build_direct_heat_source_utilisation_profiles", color = "0.52 0.6 0.85", style="rounded"];
            104[label = "solve_sector_network_myopic", color = "0.09 0.6 0.85", style="rounded"];
            105[label = "add_brownfield", color = "0.54 0.6 0.85", style="rounded"];
            106[label = "prepare_sector_network", color = "0.25 0.6 0.85", style="rounded"];
            107[label = "build_biomass_potentials\nplanning_horizons: 2050", color = "0.15 0.6 0.85", style="rounded"];
            108[label = "retrieve_cost_data\nyear: 2050", color = "0.08 0.6 0.85", style="rounded"];
            109[label = "build_industrial_energy_demand_per_node", color = "0.38 0.6 0.85", style="rounded"];
            110[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2050", color = "0.49 0.6 0.85", style="rounded"];
            111[label = "build_industrial_production_per_node", color = "0.12 0.6 0.85", style="rounded"];
            112[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2050", color = "0.07 0.6 0.85", style="rounded"];
            113[label = "build_district_heat_share\nplanning_horizons: 2050", color = "0.20 0.6 0.85", style="rounded"];
            114[label = "build_cop_profiles", color = "0.18 0.6 0.85", style="rounded"];
            115[label = "build_central_heating_temperature_profiles\nplanning_horizons: 2050", color = "0.06 0.6 0.85", style="rounded"];
            116[label = "build_direct_heat_source_utilisation_profiles", color = "0.52 0.6 0.85", style="rounded"];
            117[label = "plot_power_network_clustered", color = "0.40 0.6 0.85", style="rounded"];
            118[label = "plot_power_network", color = "0.53 0.6 0.85", style="rounded"];
            119[label = "plot_power_network", color = "0.53 0.6 0.85", style="rounded"];
            120[label = "plot_power_network", color = "0.53 0.6 0.85", style="rounded"];
            121[label = "plot_hydrogen_network", color = "0.32 0.6 0.85", style="rounded"];
            122[label = "plot_hydrogen_network", color = "0.32 0.6 0.85", style="rounded"];
            123[label = "plot_hydrogen_network", color = "0.32 0.6 0.85", style="rounded"];
            124[label = "plot_gas_network", color = "0.29 0.6 0.85", style="rounded"];
            125[label = "plot_gas_network", color = "0.29 0.6 0.85", style="rounded"];
            126[label = "plot_gas_network", color = "0.29 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            58 -> 1
            8 -> 1
            3 -> 2
            92 -> 2
            104 -> 2
            72 -> 2
            117 -> 2
            118 -> 2
            119 -> 2
            120 -> 2
            121 -> 2
            122 -> 2
            123 -> 2
            124 -> 2
            125 -> 2
            126 -> 2
            4 -> 3
            72 -> 3
            5 -> 4
            53 -> 4
            25 -> 4
            24 -> 4
            63 -> 4
            72 -> 4
            88 -> 4
            91 -> 4
            60 -> 4
            6 -> 5
            34 -> 5
            36 -> 5
            38 -> 5
            41 -> 5
            43 -> 5
            44 -> 5
            58 -> 5
            59 -> 5
            64 -> 5
            66 -> 5
            67 -> 5
            60 -> 5
            8 -> 5
            69 -> 5
            72 -> 5
            73 -> 5
            25 -> 5
            24 -> 5
            63 -> 5
            74 -> 5
            54 -> 5
            81 -> 5
            87 -> 5
            68 -> 5
            88 -> 5
            90 -> 5
            7 -> 6
            12 -> 6
            24 -> 6
            11 -> 6
            8 -> 7
            9 -> 7
            12 -> 7
            24 -> 7
            11 -> 7
            10 -> 9
            11 -> 9
            13 -> 12
            14 -> 12
            15 -> 12
            17 -> 12
            19 -> 12
            21 -> 12
            23 -> 12
            8 -> 12
            16 -> 15
            13 -> 15
            18 -> 17
            13 -> 17
            20 -> 19
            13 -> 19
            22 -> 21
            13 -> 21
            25 -> 24
            30 -> 24
            26 -> 25
            27 -> 25
            27 -> 26
            29 -> 26
            28 -> 27
            12 -> 27
            27 -> 29
            12 -> 29
            25 -> 30
            12 -> 30
            31 -> 30
            32 -> 31
            33 -> 31
            35 -> 34
            12 -> 34
            24 -> 34
            11 -> 34
            8 -> 35
            9 -> 35
            12 -> 35
            24 -> 35
            11 -> 35
            37 -> 36
            12 -> 36
            24 -> 36
            11 -> 36
            8 -> 37
            9 -> 37
            12 -> 37
            24 -> 37
            11 -> 37
            39 -> 38
            24 -> 38
            40 -> 39
            42 -> 41
            40 -> 41
            24 -> 41
            44 -> 43
            54 -> 43
            45 -> 44
            52 -> 44
            46 -> 45
            48 -> 45
            50 -> 45
            6 -> 45
            34 -> 45
            36 -> 45
            24 -> 45
            52 -> 45
            53 -> 45
            30 -> 45
            47 -> 46
            12 -> 46
            24 -> 46
            11 -> 46
            8 -> 47
            12 -> 47
            24 -> 47
            11 -> 47
            49 -> 48
            12 -> 48
            24 -> 48
            11 -> 48
            8 -> 49
            12 -> 49
            24 -> 49
            11 -> 49
            51 -> 50
            12 -> 50
            24 -> 50
            11 -> 50
            8 -> 51
            12 -> 51
            24 -> 51
            11 -> 51
            24 -> 53
            55 -> 54
            56 -> 55
            24 -> 55
            11 -> 55
            12 -> 56
            57 -> 56
            11 -> 56
            60 -> 59
            63 -> 59
            12 -> 60
            8 -> 60
            61 -> 60
            58 -> 60
            62 -> 60
            56 -> 63
            24 -> 63
            11 -> 63
            65 -> 64
            63 -> 64
            60 -> 65
            12 -> 66
            24 -> 66
            60 -> 66
            63 -> 67
            59 -> 67
            60 -> 67
            8 -> 67
            68 -> 67
            56 -> 68
            24 -> 68
            11 -> 68
            70 -> 69
            58 -> 69
            71 -> 69
            24 -> 69
            8 -> 69
            12 -> 69
            8 -> 73
            24 -> 73
            75 -> 74
            81 -> 74
            86 -> 74
            76 -> 75
            79 -> 75
            80 -> 75
            77 -> 76
            61 -> 76
            78 -> 77
            60 -> 79
            61 -> 79
            80 -> 79
            77 -> 80
            61 -> 80
            58 -> 80
            82 -> 81
            85 -> 81
            24 -> 82
            63 -> 82
            83 -> 82
            84 -> 82
            80 -> 85
            82 -> 86
            79 -> 86
            60 -> 87
            63 -> 87
            89 -> 88
            68 -> 88
            24 -> 88
            68 -> 89
            24 -> 89
            89 -> 90
            63 -> 91
            59 -> 91
            87 -> 91
            93 -> 92
            52 -> 92
            46 -> 93
            48 -> 93
            50 -> 93
            6 -> 93
            34 -> 93
            36 -> 93
            25 -> 93
            24 -> 93
            94 -> 93
            3 -> 93
            52 -> 93
            101 -> 93
            6 -> 94
            34 -> 94
            36 -> 94
            38 -> 94
            41 -> 94
            43 -> 94
            44 -> 94
            58 -> 94
            59 -> 94
            64 -> 94
            66 -> 94
            67 -> 94
            60 -> 94
            8 -> 94
            95 -> 94
            52 -> 94
            73 -> 94
            25 -> 94
            24 -> 94
            63 -> 94
            96 -> 94
            54 -> 94
            98 -> 94
            100 -> 94
            68 -> 94
            101 -> 94
            103 -> 94
            70 -> 95
            58 -> 95
            71 -> 95
            24 -> 95
            8 -> 95
            12 -> 95
            97 -> 96
            98 -> 96
            86 -> 96
            76 -> 97
            79 -> 97
            80 -> 97
            82 -> 98
            99 -> 98
            80 -> 99
            60 -> 100
            63 -> 100
            102 -> 101
            68 -> 101
            24 -> 101
            68 -> 102
            24 -> 102
            102 -> 103
            105 -> 104
            108 -> 104
            46 -> 105
            48 -> 105
            50 -> 105
            6 -> 105
            34 -> 105
            36 -> 105
            25 -> 105
            24 -> 105
            106 -> 105
            92 -> 105
            108 -> 105
            114 -> 105
            6 -> 106
            34 -> 106
            36 -> 106
            38 -> 106
            41 -> 106
            43 -> 106
            44 -> 106
            58 -> 106
            59 -> 106
            64 -> 106
            66 -> 106
            67 -> 106
            60 -> 106
            8 -> 106
            107 -> 106
            108 -> 106
            73 -> 106
            25 -> 106
            24 -> 106
            63 -> 106
            109 -> 106
            54 -> 106
            111 -> 106
            113 -> 106
            68 -> 106
            114 -> 106
            116 -> 106
            70 -> 107
            58 -> 107
            71 -> 107
            24 -> 107
            8 -> 107
            12 -> 107
            110 -> 109
            111 -> 109
            86 -> 109
            76 -> 110
            79 -> 110
            80 -> 110
            82 -> 111
            112 -> 111
            80 -> 112
            60 -> 113
            63 -> 113
            115 -> 114
            68 -> 114
            24 -> 114
            68 -> 115
            24 -> 115
            115 -> 116
            24 -> 117
            3 -> 118
            24 -> 118
            92 -> 119
            24 -> 119
            104 -> 120
            24 -> 120
            3 -> 121
            24 -> 121
            92 -> 122
            24 -> 122
            104 -> 123
            24 -> 123
            3 -> 124
            24 -> 124
            92 -> 125
            24 -> 125
            104 -> 126
            24 -> 126
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
