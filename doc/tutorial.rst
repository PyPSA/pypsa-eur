..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

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

    $ snakemake results/test-elec/networks/base_s_6_elec_.nc --configfile config/test/config.electricity.yaml

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
open-source solver HiGHS.

.. literalinclude:: ../config/test/config.electricity.yaml
   :language: yaml
   :start-at: solver:
   :end-before: check_objective:

Note, that ``config/test/config.electricity.yaml`` only includes changes relative to
the default configuration. There are many more configuration options, which are
documented at :ref:`config`.


How to use ``snakemake`` rules?
===================================

Open a terminal, go into the PyPSA-Eur directory, and activate the ``pypsa-eur`` environment with

.. code:: console

    $ conda activate pypsa-eur

Let's say based on the modifications above we would like to solve a very simplified model
clustered down to 6 buses and every 24 hours aggregated to one snapshot. The command

.. code:: console

    $ snakemake results/test-elec/networks/base_s_6_elec_.nc --configfile config/test/config.electricity.yaml

orders ``snakemake`` to run the rule :mod:`solve_network` that produces the solved network and stores it in ``results/networks`` with the name ``base_s_6_elec_.nc``:

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
            0[label = "solve_network", color = "0.15 0.6 0.85", style="rounded"];
            1[label = "prepare_network\nopts: ", color = "0.07 0.6 0.85", style="rounded"];
            2[label = "add_electricity", color = "0.49 0.6 0.85", style="rounded"];
            3[label = "build_renewable_profiles", color = "0.25 0.6 0.85", style="rounded"];
            4[label = "determine_availability_matrix\ntechnology: solar", color = "0.10 0.6 0.85", style="rounded"];
            5[label = "retrieve_databundle", color = "0.02 0.6 0.85", style="rounded"];
            6[label = "build_shapes", color = "0.05 0.6 0.85", style="rounded"];
            7[label = "retrieve_eez", color = "0.10 0.6 0.85", style="rounded"];
            8[label = "retrieve_nuts_2021_shapes", color = "0.49 0.6 0.85", style="rounded"];
            9[label = "build_osm_boundaries", color = "0.43 0.6 0.85", style="rounded"];
            10[label = "retrieve_osm_boundaries\ncountry: BA", color = "0.52 0.6 0.85", style="rounded"];
            11[label = "build_osm_boundaries", color = "0.43 0.6 0.85", style="rounded"];
            12[label = "retrieve_osm_boundaries\ncountry: MD", color = "0.52 0.6 0.85", style="rounded"];
            13[label = "build_osm_boundaries", color = "0.43 0.6 0.85", style="rounded"];
            14[label = "retrieve_osm_boundaries\ncountry: UA", color = "0.52 0.6 0.85", style="rounded"];
            15[label = "build_osm_boundaries", color = "0.43 0.6 0.85", style="rounded"];
            16[label = "retrieve_osm_boundaries\ncountry: XK", color = "0.52 0.6 0.85", style="rounded"];
            17[label = "retrieve_jrc_ardeco", color = "0.06 0.6 0.85", style="rounded"];
            18[label = "cluster_network\nclusters: 6", color = "0.12 0.6 0.85", style="rounded"];
            19[label = "simplify_network", color = "0.33 0.6 0.85", style="rounded"];
            20[label = "add_transmission_projects_and_dlr", color = "0.46 0.6 0.85", style="rounded"];
            21[label = "base_network", color = "0.00 0.6 0.85", style="rounded"];
            22[label = "retrieve_osm_prebuilt", color = "0.41 0.6 0.85", style="rounded"];
            23[label = "build_line_rating", color = "0.57 0.6 0.85", style="rounded"];
            24[label = "retrieve_cutout\ncutout: be-03-2013-era5", color = "0.18 0.6 0.85", style="rounded"];
            25[label = "build_transmission_projects", color = "0.04 0.6 0.85", style="rounded"];
            26[label = "build_electricity_demand_base", color = "0.61 0.6 0.85", style="rounded"];
            27[label = "build_electricity_demand", color = "0.16 0.6 0.85", style="rounded"];
            28[label = "retrieve_electricity_demand", color = "0.21 0.6 0.85", style="rounded"];
            29[label = "retrieve_synthetic_electricity_demand", color = "0.13 0.6 0.85", style="rounded"];
            30[label = "build_renewable_profiles", color = "0.25 0.6 0.85", style="rounded"];
            31[label = "determine_availability_matrix\ntechnology: solar-hsat", color = "0.10 0.6 0.85", style="rounded"];
            32[label = "build_renewable_profiles", color = "0.25 0.6 0.85", style="rounded"];
            33[label = "determine_availability_matrix\ntechnology: onwind", color = "0.10 0.6 0.85", style="rounded"];
            34[label = "build_renewable_profiles", color = "0.25 0.6 0.85", style="rounded"];
            35[label = "determine_availability_matrix\ntechnology: offwind-ac", color = "0.10 0.6 0.85", style="rounded"];
            36[label = "build_ship_raster", color = "0.36 0.6 0.85", style="rounded"];
            37[label = "retrieve_ship_raster", color = "0.39 0.6 0.85", style="rounded"];
            38[label = "build_renewable_profiles", color = "0.25 0.6 0.85", style="rounded"];
            39[label = "determine_availability_matrix\ntechnology: offwind-dc", color = "0.10 0.6 0.85", style="rounded"];
            40[label = "build_renewable_profiles", color = "0.25 0.6 0.85", style="rounded"];
            41[label = "determine_availability_matrix\ntechnology: offwind-float", color = "0.10 0.6 0.85", style="rounded"];
            42[label = "retrieve_cost_data\nyear: 2040", color = "0.64 0.6 0.85", style="rounded"];
            43[label = "build_powerplants", color = "0.39 0.6 0.85", style="rounded"];
            1 -> 0
            2 -> 1
            42 -> 1
            3 -> 2
            30 -> 2
            32 -> 2
            34 -> 2
            38 -> 2
            40 -> 2
            18 -> 2
            42 -> 2
            43 -> 2
            26 -> 2
            4 -> 3
            6 -> 3
            18 -> 3
            24 -> 3
            5 -> 4
            6 -> 4
            18 -> 4
            24 -> 4
            7 -> 6
            8 -> 6
            9 -> 6
            11 -> 6
            13 -> 6
            15 -> 6
            17 -> 6
            5 -> 6
            10 -> 9
            7 -> 9
            12 -> 11
            7 -> 11
            14 -> 13
            7 -> 13
            16 -> 15
            7 -> 15
            19 -> 18
            26 -> 18
            20 -> 19
            21 -> 19
            21 -> 20
            23 -> 20
            25 -> 20
            22 -> 21
            6 -> 21
            21 -> 23
            24 -> 23
            21 -> 25
            6 -> 25
            19 -> 26
            6 -> 26
            27 -> 26
            28 -> 27
            29 -> 27
            31 -> 30
            6 -> 30
            18 -> 30
            24 -> 30
            5 -> 31
            6 -> 31
            18 -> 31
            24 -> 31
            33 -> 32
            6 -> 32
            18 -> 32
            24 -> 32
            5 -> 33
            6 -> 33
            18 -> 33
            24 -> 33
            35 -> 34
            6 -> 34
            18 -> 34
            24 -> 34
            5 -> 35
            36 -> 35
            6 -> 35
            18 -> 35
            24 -> 35
            37 -> 36
            24 -> 36
            39 -> 38
            6 -> 38
            18 -> 38
            24 -> 38
            5 -> 39
            36 -> 39
            6 -> 39
            18 -> 39
            24 -> 39
            41 -> 40
            6 -> 40
            18 -> 40
            24 -> 40
            5 -> 41
            36 -> 41
            6 -> 41
            18 -> 41
            24 -> 41
            18 -> 43
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
    build_osm_boundaries                         4
    build_powerplants                            1
    build_renewable_profiles                     6
    build_shapes                                 1
    build_ship_raster                            1
    build_transmission_projects                  1
    cluster_network                              1
    determine_availability_matrix                6
    prepare_network                              1
    retrieve_cost_data                           1
    retrieve_databundle                          1
    retrieve_eez                                 1
    retrieve_electricity_demand                  1
    retrieve_jrc_ardeco                          1
    retrieve_nuts_2021_shapes                    1
    retrieve_osm_boundaries                      4
    retrieve_osm_prebuilt                        1
    retrieve_ship_raster                         1
    retrieve_synthetic_electricity_demand        1
    simplify_network                             1
    solve_network                                1
    total                                       43


``snakemake`` then runs these jobs in the correct order.

A job (here ``build_powerplants``) will display its attributes and normally some logs below this block:

.. code:: console

    rule build_powerplants:
        input: resources/test/networks/base_s_6.nc, data/custom_powerplants.csv
        output: resources/test/powerplants_s_6.csv
        log: logs/test/build_powerplants_s_6.log
        jobid: 43
        benchmark: benchmarks/test/build_powerplants_s_6
        reason: Missing output files: resources/test/powerplants_s_6.csv; Input files updated by another job: resources/test/networks/base_s_6.nc
        wildcards: clusters=6
        resources: tmpdir=<TBD>, mem_mb=7000, mem_mib=6676

Once the whole worktree is finished, it should state so in the terminal.

You will notice that many intermediate stages are saved, namely the outputs of each individual ``snakemake`` rule.

You can produce any output file occurring in the ``Snakefile`` by running

.. code:: console

    $ snakemake <output file>

For example, you can explore the evolution of the PyPSA networks by running

#. ``snakemake resources/networks/base.nc --configfile config/test/config.electricity.yaml``
#. ``snakemake resources/networks/base_s.nc --configfile config/test/config.electricity.yaml``
#. ``snakemake resources/networks/base_s_6.nc --configfile config/test/config.electricity.yaml``
#. ``snakemake resources/networks/base_s_6_elec_.nc --configfile config/test/config.electricity.yaml``

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

    n = pypsa.Network("results/networks/base_s_6_elec_.nc")

For inspiration, read the `examples section in the PyPSA documentation <https://pypsa.readthedocs.io/en/latest/examples-basic.html>`__.
