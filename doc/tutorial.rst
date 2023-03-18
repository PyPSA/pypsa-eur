..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

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
is located at ``test/config.electricity.yaml``. It includes parts deviating from
the default config file ``config.default.yaml``. To run the tutorial with this
configuration, execute

.. code:: bash
    :class: full-width

    snakemake -call --configfile test/config.electricity.yaml results/networks/elec_s_6_ec_lcopt_Co2L-24H.nc

This configuration is set to download a reduced data set via the rules :mod:`retrieve_databundle`,
:mod:`retrieve_natura_raster`, :mod:`retrieve_cutout`.
For more information on the data dependencies of PyPSA-Eur, continue reading :ref:`data`.

How to configure runs?
===========================

The model can be adapted to only include selected countries (e.g. Belgium) instead of all European countries to limit the spatial scope.

.. literalinclude:: ../test/config.electricity.yaml
   :language: yaml
   :start-at: countries:
   :end-before: snapshots:

Likewise, the example's temporal scope can be restricted (e.g. to a single week).

.. literalinclude:: ../test/config.electricity.yaml
   :language: yaml
   :start-at: snapshots:
   :end-before: electricity:

It is also possible to allow less or more carbon-dioxide emissions. Here, we limit the emissions of Belgium to 100 Mt per year.

.. literalinclude:: ../test/config.electricity.yaml
   :language: yaml
   :start-at: electricity:
   :end-before: extendable_carriers:

PyPSA-Eur also includes a database of existing conventional powerplants.
We can select which types of existing powerplants we like to be extendable:

.. literalinclude:: ../test/config.electricity.yaml
   :language: yaml
   :start-at: extendable_carriers:
   :end-before: renewable_carriers:

To accurately model the temporal and spatial availability of renewables such as
wind and solar energy, we rely on historical weather data. It is advisable to
adapt the required range of coordinates to the selection of countries.

.. literalinclude:: ../test/config.electricity.yaml
   :language: yaml
   :start-at: atlite:
   :end-before: renewable:

We can also decide which weather data source should be used to calculate
potentials and capacity factor time-series for each carrier. For example, we may
want to use the ERA-5 dataset for solar and not the default SARAH-2 dataset.

.. literalinclude:: ../test/config.electricity.yaml
   :language: yaml
   :start-at: solar:
   :end-at: cutout:

Finally, it is possible to pick a solver. For instance, this tutorial uses the
open-source solver GLPK.

.. literalinclude:: ../test/config.electricity.yaml
   :language: yaml
   :start-at: solver:
   :end-before: plotting:

Note, that ``test/config.electricity.yaml`` only includes changes relative to
the default configuration. There are many more configuration options, which are
documented at :ref:`config`.


How to use ``snakemake`` rules?
===================================

Open a terminal, go into the PyPSA-Eur directory, and activate the ``pypsa-eur`` environment with

.. code:: bash

    mamba activate pypsa-eur

Let's say based on the modifications above we would like to solve a very simplified model
clustered down to 6 buses and every 24 hours aggregated to one snapshot. The command

.. code:: bash

    snakemake -call --configfile test/config.electricity.yaml results/networks/elec_s_6_ec_lcopt_Co2L-24H.nc

orders ``snakemake`` to run the rule :mod:`solve_network` that produces the solved network and stores it in ``results/networks`` with the name ``elec_s_6_ec_lcopt_Co2L-24H.nc``:

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
        0[label = "solve_network", color = "0.21 0.6 0.85", style="rounded"];
        1[label = "prepare_network\nll: copt\nopts: Co2L-24H", color = "0.02 0.6 0.85", style="rounded"];
        2[label = "add_extra_components", color = "0.37 0.6 0.85", style="rounded"];
        3[label = "cluster_network\nclusters: 6", color = "0.39 0.6 0.85", style="rounded"];
        4[label = "simplify_network\nsimpl: ", color = "0.11 0.6 0.85", style="rounded"];
        5[label = "add_electricity", color = "0.23 0.6 0.85", style="rounded"];
        6[label = "build_renewable_profiles\ntechnology: onwind", color = "0.57 0.6 0.85", style="rounded"];
        7[label = "base_network", color = "0.09 0.6 0.85", style="rounded"];
        8[label = "build_shapes", color = "0.41 0.6 0.85", style="rounded"];
        9[label = "retrieve_databundle", color = "0.28 0.6 0.85", style="rounded"];
        10[label = "retrieve_natura_raster", color = "0.62 0.6 0.85", style="rounded"];
        11[label = "build_bus_regions", color = "0.53 0.6 0.85", style="rounded"];
        12[label = "retrieve_cutout\ncutout: europe-2013-era5", color = "0.05 0.6 0.85", style="rounded,dashed"];
        13[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.57 0.6 0.85", style="rounded"];
        14[label = "build_ship_raster", color = "0.64 0.6 0.85", style="rounded"];
        15[label = "retrieve_ship_raster", color = "0.07 0.6 0.85", style="rounded,dashed"];
        16[label = "retrieve_cutout\ncutout: europe-2013-sarah", color = "0.05 0.6 0.85", style="rounded,dashed"];
        17[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.57 0.6 0.85", style="rounded"];
        18[label = "build_renewable_profiles\ntechnology: solar", color = "0.57 0.6 0.85", style="rounded"];
        19[label = "build_hydro_profile", color = "0.44 0.6 0.85", style="rounded"];
        20[label = "retrieve_cost_data", color = "0.30 0.6 0.85", style="rounded"];
        21[label = "build_powerplants", color = "0.16 0.6 0.85", style="rounded"];
        22[label = "build_load_data", color = "0.00 0.6 0.85", style="rounded"];
        23[label = "retrieve_load_data", color = "0.34 0.6 0.85", style="rounded,dashed"];
        1 -> 0
        2 -> 1
        20 -> 1
        3 -> 2
        20 -> 2
        4 -> 3
        20 -> 3
        5 -> 4
        20 -> 4
        11 -> 4
        6 -> 5
        13 -> 5
        17 -> 5
        18 -> 5
        19 -> 5
        7 -> 5
        20 -> 5
        11 -> 5
        21 -> 5
        9 -> 5
        22 -> 5
        8 -> 5
        7 -> 6
        9 -> 6
        10 -> 6
        8 -> 6
        11 -> 6
        12 -> 6
        8 -> 7
        9 -> 8
        8 -> 11
        7 -> 11
        7 -> 13
        9 -> 13
        10 -> 13
        14 -> 13
        8 -> 13
        11 -> 13
        12 -> 13
        15 -> 14
        12 -> 14
        16 -> 14
        7 -> 17
        9 -> 17
        10 -> 17
        14 -> 17
        8 -> 17
        11 -> 17
        12 -> 17
        7 -> 18
        9 -> 18
        10 -> 18
        8 -> 18
        11 -> 18
        16 -> 18
        8 -> 19
        12 -> 19
        7 -> 21
        23 -> 22
    }

|

In the terminal, this will show up as a list of jobs to be run:

.. code:: bash

    Building DAG of jobs...
    job                         count    min threads    max threads
    ------------------------  -------  -------------  -------------
    add_electricity                 1              1              1
    add_extra_components            1              1              1
    base_network                    1              1              1
    build_bus_regions               1              1              1
    build_hydro_profile             1              1              1
    build_load_data                 1              1              1
    build_powerplants               1              1              1
    build_renewable_profiles        4              1              1
    build_shapes                    1              1              1
    build_ship_raster               1              1              1
    cluster_network                 1              1              1
    prepare_network                 1              1              1
    retrieve_cost_data              1              1              1
    retrieve_databundle             1              1              1
    retrieve_natura_raster          1              1              1
    simplify_network                1              1              1
    solve_network                   1              1              1
    total                          20              1              1


``snakemake`` then runs these jobs in the correct order.

A job (here ``simplify_network``) will display its attributes and normally some logs below this block:

.. code:: bash

    [Mon Jan 1 00:00:00 2023]
    rule simplify_network:
        input: networks/elec.nc, resources/costs.csv, resources/regions_onshore.geojson, resources/regions_offshore.geojson
        output: networks/elec_s.nc, resources/regions_onshore_elec_s.geojson, resources/regions_offshore_elec_s.geojson, resources/busmap_elec_s.csv, resources/connection_costs_s.csv
        log: logs/simplify_network/elec_s.log
        jobid: 4
        benchmark: benchmarks/simplify_network/elec_s
        reason: Missing output files: resources/busmap_elec_s.csv, resources/regions_onshore_elec_s.geojson, networks/elec_s.nc, resources/regions_offshore_elec_s.geojson; Input files updated by another job: resources/regions_offshore.geojson, resources/regions_onshore.geojson, resources/costs.csv, networks/elec.nc
        wildcards: simpl=
        resources: tmpdir=/tmp, mem_mb=4000, mem_mib=3815

Once the whole worktree is finished, it should state so in the terminal.

You will notice that many intermediate stages are saved, namely the outputs of each individual ``snakemake`` rule.

You can produce any output file occurring in the ``Snakefile`` by running

.. code:: bash

    snakemake -call <output file>

For example, you can explore the evolution of the PyPSA networks by running

#. ``snakemake -call --configfile test/config.electricity.yaml resources/networks/base.nc``
#. ``snakemake -call --configfile test/config.electricity.yaml resources/networks/elec.nc``
#. ``snakemake -call --configfile test/config.electricity.yaml resources/networks/elec_s.nc``
#. ``snakemake -call --configfile test/config.electricity.yaml resources/networks/elec_s_6.nc``
#. ``snakemake -call --configfile test/config.electricity.yaml resources/networks/elec_s_6_ec_lcopt_Co2L-24H.nc``

To run all combinations of wildcard values provided in the ``config.yaml`` under ``scenario:``,
you can use the collection rule ``solve_elec_networks``.

.. code:: bash

    snakemake -call --configfile test/config.electricity.yaml solve_elec_networks

If you now feel confident and want to tackle runs with larger temporal and
spatial scope, clean-up the repository and after modifying the ``config.yaml`` file
target the collection rule ``solve_elec_networks`` again without providing the test
configuration file.

.. code:: bash

    snakemake -call purge
    snakemake -call solve_elec_networks

.. note::

    It is good practice to perform a dry-run using the option `-n`, before you
    commit to a run:

    .. code:: bash

        snakemake -call solve_elec_networks -n

How to analyse results?
===============================

The solved networks can be analysed just like any other PyPSA network (e.g. in
Jupyter Notebooks).

.. code:: python

    import pypsa

    n = pypsa.Network("results/networks/elec_s_6_ec_lcopt_Co2L-24H.nc")

For inspiration, read the `examples section in the PyPSA documentation <https://pypsa.readthedocs.io/en/latest/examples-basic.html>`_.
