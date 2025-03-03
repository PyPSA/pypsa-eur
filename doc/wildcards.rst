..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

  SPDX-License-Identifier: CC-BY-4.0

.. _wildcards:

#########
Wildcards
#########

It is easy to run PyPSA-Eur for multiple scenarios using the wildcards feature of ``snakemake``.
Wildcards allow to generalise a rule to produce all files that follow a regular expression pattern
which e.g. defines one particular scenario. One can think of a wildcard as a parameter that shows
up in the input/output file names of the ``Snakefile`` and thereby determines which rules to run,
what data to retrieve and what files to produce.

.. note::
    Detailed explanations of how wildcards work in ``snakemake`` can be found in the
    `relevant section of the documentation <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards>`__.

.. _cutout_wc:

The ``{cutout}`` wildcard
=========================

The ``{cutout}`` wildcard facilitates running the rule :mod:`build_cutout`
for all cutout configurations specified under ``atlite: cutouts:``.
These cutouts will be stored in a folder specified by ``{cutout}``.

.. _technology:

The ``{technology}`` wildcard
=============================

The ``{technology}`` wildcard specifies for which renewable energy technology to produce availability time
series and potentials using the rule :mod:`build_renewable_profiles`.
It can take the values ``onwind``, ``offwind-ac``, ``offwind-dc``, ``offwind-float``, and ``solar`` but **not** ``hydro``
(since hydroelectric plant profiles are created by a different rule)``

.. _clusters:

The ``{clusters}`` wildcard
===========================

The ``{clusters}`` wildcard specifies the number of buses a detailed
network model should be reduced to in the rule :mod:`cluster_network`.
The number of clusters must be lower than the total number of nodes
and higher than the number of countries. However, a country counts twice if
it has two asynchronous subnetworks (e.g. Denmark or Italy).

.. _opts:

The ``{opts}`` wildcard
=======================

The ``{opts}`` wildcard is used for electricity-only studies. It triggers
optional constraints, which are activated in either :mod:`prepare_network` or
the :mod:`solve_network` step. It may hold multiple triggers separated by ``-``,
i.e. ``Co2L-3h`` contains the ``Co2L`` trigger and the ``3h`` switch. There are
currently:


.. csv-table::
   :header-rows: 1
   :widths: 10,20,10,10
   :file: configtables/opts.csv

.. _sector_opts:

The ``{sector_opts}`` wildcard
==============================

.. warning::
    More comprehensive documentation for this wildcard will be added soon.
    To really understand the options here, look in scripts/prepare_sector_network.py

The ``{sector_opts}`` wildcard is only used for sector-coupling studies.

.. csv-table::
   :header-rows: 1
   :widths: 10,20,10,10
   :file: configtables/sector-opts.csv

.. _planning_horizons:

The ``{planning_horizons}`` wildcard
====================================

.. warning::
    More comprehensive documentation for this wildcard will be added soon.

The ``{planning_horizons}`` wildcard is only used for sector-coupling studies.
It takes years as values, e.g. 2020, 2030, 2040, 2050.
