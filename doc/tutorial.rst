..
  SPDX-FileCopyrightText: 2019-2024 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _tutorial:

###############################
Tutorial: Electricity-Only
###############################

.. raw:: html

    <iframe width="832" height="468" src="https://www.youtube.com/embed/mAwhQnNRIvs" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

.. note::
    If you have not done it yet, follow the :ref:`installation` steps first.

In this tutorial, we will build a heavily simplified power system model for
Belgium. But before getting started with **PyPSA-Eur** it makes sense to be familiar
with its general modelling framework `PyPSA <https://pypsa.readthedocs.io>`__.

Running the tutorial requires limited computational resources compared to the
full model, which allows the user to explore most of its functionalities on a
local machine. The tutorial will cover examples on how to configure and
customise the PyPSA-Eur model and run the ``snakemake`` workflow step by step
from network creation to the solved network. The configuration for the tutorial
is located at ``config/test/config.electricity.yaml``. It includes parts deviating from
the default config file ``config/config.default.yaml``. To run the tutorial with this
configuration, execute

.. code:: console
    :class: full-width

    $ snakemake results/test-elec/networks/base_s_6_elec_lcopt_.nc --configfile config/test/config.electricity.yaml

This configuration is set to download a reduced cutout via the rule :mod:`retrieve_cutout`.
For more information on the data dependencies of PyPSA-Eur, continue reading :ref:`data`.

How to configure runs?
===========================

The model can be adapted to only include selected countries (e.g. Belgium) instead of all European countries to limit the spatial scope.

.. literalinclude:: ../config/test/config.electricity.yaml
   :language: yaml
   :start-at: countries:
   :end-before: snapshots:

Likewise, the example's temporal scope can be restricted (e.g. to a single week).

.. literalinclude:: ../config/test/config.electricity.yaml
   :language: yaml
   :start-at: snapshots:
   :end-before: electricity:

It is also possible to allow less or more carbon-dioxide emissions. Here, we limit the emissions of Belgium to 100 Mt per year.

.. literalinclude:: ../config/test/config.electricity.yaml
   :language: yaml
   :start-at: electricity:
   :end-before: extendable_carriers:

PyPSA-Eur also includes a database of existing conventional powerplants.
We can select which types of existing powerplants we like to be extendable:

.. literalinclude:: ../config/test/config.electricity.yaml
   :language: yaml
   :start-at: extendable_carriers:
   :end-before: renewable_carriers:

To accurately model the temporal and spatial availability of renewables such as
wind and solar energy, we rely on historical weather data. It is advisable to
adapt the required range of coordinates to the selection of countries.

.. literalinclude:: ../config/test/config.electricity.yaml
   :language: yaml
   :start-at: atlite:
   :end-before: renewable:

We can also decide which weather data source should be used to calculate
potentials and capacity factor time-series for each carrier. For example, we may
want to use the ERA-5 dataset for solar and not the default SARAH-3 dataset.

.. literalinclude:: ../config/test/config.electricity.yaml
   :language: yaml
   :start-at: solar:
   :end-at: cutout:

Finally, it is possible to pick a solver. For instance, this tutorial uses the
open-source solver GLPK.

.. literalinclude:: ../config/test/config.electricity.yaml
   :language: yaml
   :start-at: solver:
   :end-before: plotting:

Note, that ``config/test/config.electricity.yaml`` only includes changes relative to
the default configuration. There are many more configuration options, which are
documented at :ref:`config`.


How to use ``snakemake`` rules?
===================================

Open a terminal, go into the PyPSA-Eur directory, and activate the ``pypsa-eur`` environment with

.. code:: console

    $ mamba activate pypsa-eur

Let's say based on the modifications above we would like to solve a very simplified model
clustered down to 6 buses and every 24 hours aggregated to one snapshot. The command

.. code:: console

    $ snakemake results/test-elec/networks/base_s_6_elec_lcopt_.nc --configfile config/test/config.electricity.yaml

orders ``snakemake`` to run the rule :mod:`solve_network` that produces the solved network and stores it in ``results/networks`` with the name ``base_s_6_elec_lcopt_.nc``:

.. literalinclude:: ../rules/solve_electricity.smk
   :start-at: rule solve_network:
   :end-before: rule solve_operations_network:

This triggers a workflow of multiple preceding jobs that depend on each rule's inputs and outputs:

.. graphviz::
    :class: full-width
    :align: center

    digraph snakemake_dag {
        graph[bgcolor=white, margin=0];
        node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
        edge[penwidth=2, color=grey];
            0[label = "solve_network", color = "0.19 0.6 0.85", style="rounded"];
            1[label = "prepare_network\nll: copt\nopts: ", color = "0.24 0.6 0.85", style="rounded"];
            2[label = "add_electricity", color = "0.35 0.6 0.85", style="rounded"];
            3[label = "build_renewable_profiles", color = "0.15 0.6 0.85", style="rounded"];
            4[label = "determine_availability_matrix\ntechnology: solar", color = "0.39 0.6 0.85", style="rounded"];
            5[label = "retrieve_databundle", color = "0.65 0.6 0.85", style="rounded"];
            6[label = "build_shapes", color = "0.45 0.6 0.85", style="rounded"];
            7[label = "retrieve_naturalearth_countries", color = "0.03 0.6 0.85", style="rounded"];
            8[label = "retrieve_eez", color = "0.17 0.6 0.85", style="rounded"];
            9[label = "cluster_network\nclusters: 6", color = "0.38 0.6 0.85", style="rounded"];
            10[label = "simplify_network", color = "0.14 0.6 0.85", style="rounded"];
            11[label = "add_transmission_projects_and_dlr", color = "0.61 0.6 0.85", style="rounded"];
            12[label = "base_network", color = "0.36 0.6 0.85", style="rounded"];
            13[label = "retrieve_osm_prebuilt", color = "0.22 0.6 0.85", style="rounded"];
            14[label = "build_line_rating", color = "0.50 0.6 0.85", style="rounded"];
            15[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.02 0.6 0.85", style="rounded"];
            16[label = "build_transmission_projects", color = "0.08 0.6 0.85", style="rounded"];
            17[label = "build_electricity_demand_base", color = "0.11 0.6 0.85", style="rounded"];
            18[label = "build_electricity_demand", color = "0.60 0.6 0.85", style="rounded"];
            19[label = "retrieve_electricity_demand", color = "0.60 0.6 0.85", style="rounded"];
            20[label = "retrieve_synthetic_electricity_demand", color = "0.32 0.6 0.85", style="rounded"];
            21[label = "build_renewable_profiles", color = "0.15 0.6 0.85", style="rounded"];
            22[label = "determine_availability_matrix\ntechnology: solar-hsat", color = "0.39 0.6 0.85", style="rounded"];
            23[label = "build_renewable_profiles", color = "0.15 0.6 0.85", style="rounded"];
            24[label = "determine_availability_matrix\ntechnology: onwind", color = "0.39 0.6 0.85", style="rounded"];
            25[label = "build_renewable_profiles", color = "0.15 0.6 0.85", style="rounded"];
            26[label = "determine_availability_matrix\ntechnology: offwind-ac", color = "0.39 0.6 0.85", style="rounded"];
            27[label = "build_ship_raster", color = "0.12 0.6 0.85", style="rounded"];
            28[label = "retrieve_ship_raster", color = "0.44 0.6 0.85", style="rounded"];
            29[label = "build_renewable_profiles", color = "0.15 0.6 0.85", style="rounded"];
            30[label = "determine_availability_matrix\ntechnology: offwind-dc", color = "0.39 0.6 0.85", style="rounded"];
            31[label = "build_renewable_profiles", color = "0.15 0.6 0.85", style="rounded"];
            32[label = "determine_availability_matrix\ntechnology: offwind-float", color = "0.39 0.6 0.85", style="rounded"];
            33[label = "retrieve_cost_data\nyear: 2030", color = "0.01 0.6 0.85", style="rounded"];
            34[label = "build_powerplants", color = "0.52 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            33 -> 1
            3 -> 2
            21 -> 2
            23 -> 2
            25 -> 2
            29 -> 2
            31 -> 2
            9 -> 2
            33 -> 2
            34 -> 2
            17 -> 2
            4 -> 3
            6 -> 3
            9 -> 3
            15 -> 3
            5 -> 4
            6 -> 4
            9 -> 4
            15 -> 4
            7 -> 6
            8 -> 6
            5 -> 6
            10 -> 9
            17 -> 9
            11 -> 10
            12 -> 10
            12 -> 11
            14 -> 11
            16 -> 11
            13 -> 12
            6 -> 12
            12 -> 14
            15 -> 14
            12 -> 16
            6 -> 16
            10 -> 17
            6 -> 17
            18 -> 17
            19 -> 18
            20 -> 18
            22 -> 21
            6 -> 21
            9 -> 21
            15 -> 21
            5 -> 22
            6 -> 22
            9 -> 22
            15 -> 22
            24 -> 23
            6 -> 23
            9 -> 23
            15 -> 23
            5 -> 24
            6 -> 24
            9 -> 24
            15 -> 24
            26 -> 25
            6 -> 25
            9 -> 25
            15 -> 25
            5 -> 26
            27 -> 26
            6 -> 26
            9 -> 26
            15 -> 26
            28 -> 27
            15 -> 27
            30 -> 29
            6 -> 29
            9 -> 29
            15 -> 29
            5 -> 30
            27 -> 30
            6 -> 30
            9 -> 30
            15 -> 30
            32 -> 31
            6 -> 31
            9 -> 31
            15 -> 31
            5 -> 32
            27 -> 32
            6 -> 32
            9 -> 32
            15 -> 32
            9 -> 34
    }

|

In the terminal, this will show up as a list of jobs to be run:

.. code:: console

    Building DAG of jobs...
    Job stats:
    job                                      count
    -------------------------------------  -------
    add_electricity                              1
    add_transmission_projects_and_dlr            1
    base_network                                 1
    build_electricity_demand                     1
    build_electricity_demand_base                1
    build_line_rating                            1
    build_powerplants                            1
    build_renewable_profiles                     6
    build_shapes                                 1
    build_ship_raster                            1
    build_transmission_projects                  1
    cluster_network                              1
    determine_availability_matrix                6
    prepare_network                              1
    retrieve_cost_data                           1
    retrieve_cutout                              1
    retrieve_databundle                          1
    retrieve_eez                                 1
    retrieve_electricity_demand                  1
    retrieve_naturalearth_countries              1
    retrieve_osm_prebuilt                        1
    retrieve_ship_raster                         1
    retrieve_synthetic_electricity_demand        1
    simplify_network                             1
    solve_network                                1
    total                                       35


``snakemake`` then runs these jobs in the correct order.

A job (here ``simplify_network``) will display its attributes and normally some logs below this block:

.. code:: console

    rule simplify_network:
        input: resources/test/networks/base_extended.nc, resources/test/regions_onshore.geojson, resources/test/regions_offshore.geojson
        output: resources/test/networks/base_s.nc, resources/test/regions_onshore_base_s.geojson, resources/test/regions_offshore_base_s.geojson, resources/test/busmap_base_s.csv
        log: logs/test/simplify_network.log
        jobid: 10
        benchmark: benchmarks/test/simplify_network_b
        reason: Forced execution
        resources: tmpdir=<TBD>, mem_mb=12000, mem_mib=11445

Once the whole worktree is finished, it should state so in the terminal.

You will notice that many intermediate stages are saved, namely the outputs of each individual ``snakemake`` rule.

You can produce any output file occurring in the ``Snakefile`` by running

.. code:: console

    $ snakemake <output file>

For example, you can explore the evolution of the PyPSA networks by running

#. ``snakemake resources/networks/base.nc --configfile config/test/config.electricity.yaml``
#. ``snakemake resources/networks/base_s.nc --configfile config/test/config.electricity.yaml``
#. ``snakemake resources/networks/base_s_6.nc --configfile config/test/config.electricity.yaml``
#. ``snakemake resources/networks/base_s_6_elec_lcopt_.nc --configfile config/test/config.electricity.yaml``

To run all combinations of wildcard values provided in the ``config/config.yaml`` under ``scenario:``,
you can use the collection rule ``solve_elec_networks``.

.. code:: console

    $ snakemake solve_elec_networks --configfile config/test/config.electricity.yaml

If you now feel confident and want to tackle runs with larger temporal and
spatial scope, clean-up the repository and after modifying the ``config/config.yaml`` file
target the collection rule ``solve_elec_networks`` again without providing the test
configuration file.

.. code:: console

    $ snakemake purge
    $ snakemake solve_elec_networks

.. note::

    It is good practice to perform a dry-run using the option `-n`, before you
    commit to a run:

    .. code:: console

        $ snakemake solve_elec_networks -n

How to analyse results?
===============================

The solved networks can be analysed just like any other PyPSA network (e.g. in
Jupyter Notebooks).

.. code:: python

    import pypsa

    n = pypsa.Network("results/networks/base_s_6_elec_lcopt_.nc")

For inspiration, read the `examples section in the PyPSA documentation <https://pypsa.readthedocs.io/en/latest/examples-basic.html>`__.
