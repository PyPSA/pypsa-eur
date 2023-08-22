..
  SPDX-FileCopyrightText: 2023 The PyPSA-Eur Authors

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

which will result in the following *additional* jobs ``snakemake`` wants to run
on top of those already included in the electricity-only tutorial:

.. code:: bash

    job                                                 count    min threads    max threads
    ------------------------------------------------  -------  -------------  -------------
    all                                                     1              1              1
    build_ammonia_production                                1              1              1
    build_biomass_potentials                                1              1              1
    build_clustered_population_layouts                      1              1              1
    build_cop_profiles                                      1              1              1
    build_gas_input_locations                               1              1              1
    build_gas_network                                       1              1              1
    build_heat_demands                                      3              1              1
    build_industrial_distribution_key                       1              1              1
    build_industrial_energy_demand_per_country_today        1              1              1
    build_industrial_energy_demand_per_node                 1              1              1
    build_industrial_energy_demand_per_node_today           1              1              1
    build_industrial_production_per_country                 1              1              1
    build_industrial_production_per_country_tomorrow        1              1              1
    build_industrial_production_per_node                    1              1              1
    build_industry_sector_ratios                            1              1              1
    build_population_weighted_energy_totals                 1              1              1
    build_salt_cavern_potentials                            1              1              1
    build_shipping_demand                                   1              1              1
    build_simplified_population_layouts                     1              1              1
    build_solar_thermal_profiles                            3              1              1
    build_temperature_profiles                              3              1              1
    build_transport_demand                                  1              1              1
    cluster_gas_network                                     1              1              1
    cluster_network                                         1              1              1
    copy_config                                             1              1              1
    make_summary                                            1              1              1
    plot_network                                            1              1              1
    plot_summary                                            1              1              1
    prepare_sector_network                                  1              1              1
    retrieve_gas_infrastructure_data                        1              1              1
    retrieve_sector_databundle                              1              1              1
    solve_sector_network                                    1              1              1

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
        0[label = "all", color = "0.51 0.6 0.85", style="rounded"];
        1[label = "plot_summary", color = "0.54 0.6 0.85", style="rounded"];
        2[label = "make_summary", color = "0.44 0.6 0.85", style="rounded"];
        3[label = "solve_sector_network", color = "0.46 0.6 0.85", style="rounded"];
        4[label = "prepare_sector_network", color = "0.09 0.6 0.85", style="rounded"];
        5[label = "cluster_gas_network", color = "0.38 0.6 0.85", style="rounded"];
        6[label = "build_gas_network", color = "0.00 0.6 0.85", style="rounded"];
        7[label = "retrieve_gas_infrastructure_data", color = "0.33 0.6 0.85", style="rounded"];
        8[label = "cluster_network", color = "0.26 0.6 0.85", style="rounded"];
        9[label = "simplify_network", color = "0.03 0.6 0.85", style="rounded"];
        10[label = "add_electricity", color = "0.25 0.6 0.85", style="rounded"];
        11[label = "build_renewable_profiles", color = "0.07 0.6 0.85", style="rounded"];
        12[label = "base_network", color = "0.16 0.6 0.85", style="rounded"];
        13[label = "build_shapes", color = "0.65 0.6 0.85", style="rounded"];
        14[label = "retrieve_databundle", color = "0.20 0.6 0.85", style="rounded"];
        15[label = "retrieve_natura_raster", color = "0.10 0.6 0.85", style="rounded"];
        16[label = "build_bus_regions", color = "0.11 0.6 0.85", style="rounded"];
        17[label = "build_ship_raster", color = "0.56 0.6 0.85", style="rounded"];
        18[label = "retrieve_ship_raster", color = "0.15 0.6 0.85", style="rounded"];
        19[label = "retrieve_cost_data", color = "0.50 0.6 0.85", style="rounded"];
        20[label = "build_powerplants", color = "0.49 0.6 0.85", style="rounded"];
        21[label = "build_electricity_demand", color = "0.39 0.6 0.85", style="rounded"];
        22[label = "retrieve_electricity_demand", color = "0.05 0.6 0.85", style="rounded"];
        23[label = "build_gas_input_locations", color = "0.45 0.6 0.85", style="rounded"];
        24[label = "prepare_network", color = "0.31 0.6 0.85", style="rounded"];
        25[label = "add_extra_components", color = "0.23 0.6 0.85", style="rounded"];
        26[label = "build_energy_totals", color = "0.19 0.6 0.85", style="rounded"];
        27[label = "build_population_weighted_energy_totals", color = "0.27 0.6 0.85", style="rounded"];
        28[label = "build_clustered_population_layouts", color = "0.64 0.6 0.85", style="rounded"];
        29[label = "build_population_layouts", color = "0.43 0.6 0.85", style="rounded"];
        30[label = "build_shipping_demand", color = "0.57 0.6 0.85", style="rounded"];
        31[label = "build_transport_demand", color = "0.53 0.6 0.85", style="rounded"];
        32[label = "build_temperature_profiles", color = "0.58 0.6 0.85", style="rounded"];
        33[label = "build_biomass_potentials", color = "0.30 0.6 0.85", style="rounded"];
        34[label = "build_salt_cavern_potentials", color = "0.47 0.6 0.85", style="rounded"];
        35[label = "build_simplified_population_layouts", color = "0.32 0.6 0.85", style="rounded"];
        36[label = "build_industrial_energy_demand_per_node", color = "0.14 0.6 0.85", style="rounded"];
        37[label = "build_industry_sector_ratios", color = "0.18 0.6 0.85", style="rounded"];
        38[label = "build_ammonia_production", color = "0.48 0.6 0.85", style="rounded"];
        39[label = "build_industrial_production_per_node", color = "0.12 0.6 0.85", style="rounded"];
        40[label = "build_industrial_distribution_key", color = "0.61 0.6 0.85", style="rounded"];
        41[label = "build_industrial_production_per_country_tomorrow", color = "0.22 0.6 0.85", style="rounded"];
        42[label = "build_industrial_production_per_country", color = "0.59 0.6 0.85", style="rounded"];
        43[label = "build_industrial_energy_demand_per_node_today", color = "0.62 0.6 0.85", style="rounded"];
        44[label = "build_industrial_energy_demand_per_country_today", color = "0.41 0.6 0.85", style="rounded"];
        45[label = "build_heat_demands", color = "0.08 0.6 0.85", style="rounded"];
        46[label = "build_cop_profiles", color = "0.52 0.6 0.85", style="rounded"];
        47[label = "build_solar_thermal_profiles", color = "0.17 0.6 0.85", style="rounded"];
        48[label = "copy_config", color = "0.40 0.6 0.85", style="rounded"];
        49[label = "plot_network", color = "0.60 0.6 0.85", style="rounded"];
        1 -> 0
        2 -> 1
        49 -> 2
        19 -> 2
        3 -> 2
        48 -> 3
        4 -> 3
        19 -> 3
        9 -> 4
        11 -> 4
        45 -> 4
        36 -> 4
        47 -> 4
        26 -> 4
        27 -> 4
        8 -> 4
        33 -> 4
        24 -> 4
        35 -> 4
        5 -> 4
        23 -> 4
        34 -> 4
        19 -> 4
        31 -> 4
        46 -> 4
        30 -> 4
        32 -> 4
        28 -> 4
        6 -> 5
        8 -> 5
        7 -> 6
        19 -> 8
        9 -> 8
        19 -> 9
        10 -> 9
        16 -> 9
        14 -> 10
        21 -> 10
        20 -> 10
        19 -> 10
        11 -> 10
        16 -> 10
        13 -> 10
        12 -> 10
        14 -> 11
        17 -> 11
        15 -> 11
        16 -> 11
        12 -> 11
        13 -> 11
        13 -> 12
        14 -> 13
        12 -> 16
        13 -> 16
        18 -> 17
        12 -> 20
        22 -> 21
        8 -> 23
        7 -> 23
        25 -> 24
        19 -> 24
        19 -> 25
        8 -> 25
        13 -> 26
        28 -> 27
        26 -> 27
        8 -> 28
        29 -> 28
        13 -> 29
        13 -> 30
        8 -> 30
        26 -> 30
        32 -> 31
        28 -> 31
        27 -> 31
        26 -> 31
        8 -> 32
        29 -> 32
        13 -> 33
        14 -> 33
        8 -> 33
        8 -> 34
        9 -> 35
        29 -> 35
        37 -> 36
        39 -> 36
        43 -> 36
        38 -> 37
        41 -> 39
        40 -> 39
        28 -> 40
        8 -> 40
        42 -> 41
        38 -> 42
        44 -> 43
        40 -> 43
        38 -> 44
        42 -> 44
        8 -> 45
        29 -> 45
        32 -> 46
        8 -> 47
        29 -> 47
        8 -> 49
        3 -> 49
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

which will result in the following *additional* jobs ``snakemake`` wants to run:

.. code:: bash

    job                                                 count    min threads    max threads
    ------------------------------------------------  -------  -------------  -------------
    all                                                     1              1              1
    add_brownfield                                          2              1              1
    add_existing_baseyear                                   1              1              1
    plot_network                                            3              1              1
    plot_summary                                            1              1              1
    prepare_sector_network                                  3              1              1
    solve_sector_network_myopic                             3              1              1

which translates to the following workflow diagram which nicely outlines
how the sequential pathway optimisation with myopic foresight is
implemented in the workflow:

.. graphviz::
    :class: full-width
    :align: center

    digraph snakemake_dag {
        graph[bgcolor=white, margin=0];
        node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
        edge[penwidth=2, color=grey];
        0[label = "all", color = "0.38 0.6 0.85", style="rounded"];
        1[label = "plot_summary", color = "0.61 0.6 0.85", style="rounded"];
        2[label = "make_summary", color = "0.51 0.6 0.85", style="rounded"];
        3[label = "solve_sector_network_myopic", color = "0.32 0.6 0.85", style="rounded"];
        4[label = "add_existing_baseyear", color = "0.20 0.6 0.85", style="rounded"];
        5[label = "prepare_sector_network", color = "0.14 0.6 0.85", style="rounded"];
        6[label = "prepare_network", color = "0.06 0.6 0.85", style="rounded"];
        7[label = "add_extra_components", color = "0.00 0.6 0.85", style="rounded"];
        8[label = "cluster_network", color = "0.18 0.6 0.85", style="rounded"];
        9[label = "simplify_network", color = "0.30 0.6 0.85", style="rounded"];
        10[label = "add_electricity", color = "0.24 0.6 0.85", style="rounded"];
        11[label = "build_renewable_profiles", color = "0.40 0.6 0.85", style="rounded"];
        12[label = "base_network", color = "0.11 0.6 0.85", style="rounded"];
        13[label = "build_shapes", color = "0.29 0.6 0.85", style="rounded"];
        14[label = "retrieve_databundle", color = "0.58 0.6 0.85", style="rounded"];
        15[label = "retrieve_natura_raster", color = "0.39 0.6 0.85", style="rounded"];
        16[label = "build_bus_regions", color = "0.60 0.6 0.85", style="rounded"];
        17[label = "build_ship_raster", color = "0.65 0.6 0.85", style="rounded"];
        18[label = "retrieve_ship_raster", color = "0.09 0.6 0.85", style="rounded"];
        19[label = "retrieve_cost_data", color = "0.04 0.6 0.85", style="rounded"];
        20[label = "build_powerplants", color = "0.28 0.6 0.85", style="rounded"];
        21[label = "build_electricity_demand", color = "0.46 0.6 0.85", style="rounded"];
        22[label = "retrieve_electricity_demand", color = "0.44 0.6 0.85", style="rounded"];
        23[label = "build_energy_totals", color = "0.53 0.6 0.85", style="rounded"];
        24[label = "build_population_weighted_energy_totals", color = "0.03 0.6 0.85", style="rounded"];
        25[label = "build_clustered_population_layouts", color = "0.34 0.6 0.85", style="rounded"];
        26[label = "build_population_layouts", color = "0.63 0.6 0.85", style="rounded"];
        27[label = "build_shipping_demand", color = "0.05 0.6 0.85", style="rounded"];
        28[label = "build_transport_demand", color = "0.52 0.6 0.85", style="rounded"];
        29[label = "build_temperature_profiles", color = "0.16 0.6 0.85", style="rounded"];
        30[label = "build_biomass_potentials", color = "0.47 0.6 0.85", style="rounded"];
        31[label = "build_salt_cavern_potentials", color = "0.48 0.6 0.85", style="rounded"];
        32[label = "build_simplified_population_layouts", color = "0.08 0.6 0.85", style="rounded"];
        33[label = "build_industrial_energy_demand_per_node", color = "0.22 0.6 0.85", style="rounded"];
        34[label = "build_industry_sector_ratios", color = "0.56 0.6 0.85", style="rounded"];
        35[label = "build_ammonia_production", color = "0.57 0.6 0.85", style="rounded"];
        36[label = "build_industrial_production_per_node", color = "0.66 0.6 0.85", style="rounded"];
        37[label = "build_industrial_distribution_key", color = "0.41 0.6 0.85", style="rounded"];
        38[label = "build_industrial_production_per_country_tomorrow", color = "0.54 0.6 0.85", style="rounded"];
        39[label = "build_industrial_production_per_country", color = "0.10 0.6 0.85", style="rounded"];
        40[label = "build_industrial_energy_demand_per_node_today", color = "0.55 0.6 0.85", style="rounded"];
        41[label = "build_industrial_energy_demand_per_country_today", color = "0.35 0.6 0.85", style="rounded"];
        42[label = "build_heat_demands", color = "0.49 0.6 0.85", style="rounded"];
        43[label = "build_cop_profiles", color = "0.01 0.6 0.85", style="rounded"];
        44[label = "build_solar_thermal_profiles", color = "0.45 0.6 0.85", style="rounded"];
        45[label = "copy_config", color = "0.33 0.6 0.85", style="rounded"];
        46[label = "add_brownfield", color = "0.59 0.6 0.85", style="rounded"];
        47[label = "plot_network", color = "0.15 0.6 0.85", style="rounded"];
        1 -> 0
        2 -> 1
        3 -> 2
        19 -> 2
        47 -> 2
        46 -> 3
        19 -> 3
        4 -> 3
        45 -> 3
        43 -> 4
        19 -> 4
        20 -> 4
        9 -> 4
        5 -> 4
        25 -> 4
        8 -> 4
        28 -> 5
        23 -> 5
        11 -> 5
        33 -> 5
        24 -> 5
        43 -> 5
        19 -> 5
        27 -> 5
        6 -> 5
        31 -> 5
        32 -> 5
        44 -> 5
        9 -> 5
        30 -> 5
        25 -> 5
        29 -> 5
        42 -> 5
        8 -> 5
        7 -> 6
        19 -> 6
        19 -> 7
        8 -> 7
        9 -> 8
        19 -> 8
        10 -> 9
        19 -> 9
        16 -> 9
        11 -> 10
        19 -> 10
        14 -> 10
        20 -> 10
        12 -> 10
        21 -> 10
        16 -> 10
        13 -> 10
        15 -> 11
        14 -> 11
        13 -> 11
        12 -> 11
        16 -> 11
        17 -> 11
        13 -> 12
        14 -> 13
        13 -> 16
        12 -> 16
        18 -> 17
        12 -> 20
        22 -> 21
        13 -> 23
        25 -> 24
        23 -> 24
        8 -> 25
        26 -> 25
        13 -> 26
        13 -> 27
        23 -> 27
        8 -> 27
        24 -> 28
        25 -> 28
        29 -> 28
        23 -> 28
        8 -> 29
        26 -> 29
        13 -> 30
        14 -> 30
        8 -> 30
        8 -> 31
        9 -> 32
        26 -> 32
        34 -> 33
        36 -> 33
        40 -> 33
        35 -> 34
        37 -> 36
        38 -> 36
        25 -> 37
        8 -> 37
        39 -> 38
        35 -> 39
        41 -> 40
        37 -> 40
        39 -> 41
        35 -> 41
        8 -> 42
        26 -> 42
        29 -> 43
        8 -> 44
        26 -> 44
        3 -> 46
        19 -> 46
        5 -> 46
        43 -> 46
        3 -> 47
        8 -> 47
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
