.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _wildcards:

#########
Wildcards
#########

.. warning::

   **MAJOR CHANGES:** The wildcard-based scenario system relies on direct
   configuration.

   Key elements:

   - ``{horizon}`` for planning horizons
   - Direct configuration in ``config.yaml`` instead of wildcard combinations
   - Simplified file naming (e.g., ``clustered.nc``)

   .. note::

      Earlier releases exposed ``{clusters}``, ``{opts}``, and ``{sector_opts}``
      wildcards. See :doc:`migration` and :doc:`release_notes` for guidance on
      translating legacy targets.

.. note::

   A few common filename transitions are summarised below; see :doc:`migration`
   for the complete table.

   .. list-table::
      :header-rows: 1

      * - Legacy target
        - New target
        - Stage
      * - ``networks/base_s_{clusters}.nc``
        - ``networks/simplified.nc``
        - Simplification
      * - ``busmap_base_s_{clusters}.csv`` / ``linemap_base_s_{clusters}.csv``
        - ``busmap.csv`` / ``linemap.csv``
        - Clustering
      * - ``networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc``
        - ``networks/composed_{horizon}.nc``
        - Composition
      * - ``RESULTS/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc``
        - ``RESULTS/networks/solved_{horizon}.nc``
        - Solving / results

PyPSA-Eur uses the wildcards feature of ``snakemake`` to run multiple scenarios.
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

The ``{clusters}`` wildcard (DEPRECATED)
=========================================

.. note::

   **DEPRECATED:** This wildcard has been removed. Use
   ``clustering.cluster_network.n_clusters`` in the config file instead.
   Legacy behaviour: the ``{clusters}`` wildcard specified the number of buses
   a detailed network model should be reduced to in the rule
   :mod:`cluster_network`.

.. _opts:

The ``{opts}`` wildcard (DEPRECATED)
=====================================

.. note::

   **DEPRECATED:** This wildcard has been removed. Configure options directly
   in their respective config sections (e.g., ``electricity``, ``clustering``).
   Legacy behaviour: the ``{opts}`` wildcard triggered optional constraints in
   the electricity-only workflow and passed combinations to
   ``prepare_network``/``solve_network``. Configure these options directly so
   :mod:`compose_network` can apply them before :mod:`solve_network` runs.

.. _sector_opts:

The ``{sector_opts}`` wildcard (DEPRECATED)
============================================

.. note::

   **DEPRECATED:** This wildcard has been removed. Configure options directly
   in the ``sector`` config section. Legacy behaviour: the
   ``{sector_opts}`` wildcard controlled sector-coupling studies and is fully
   represented by the dedicated config keys.

.. _planning_horizons:

The ``{horizon}`` wildcard
===========================

The ``{horizon}`` wildcard is used for multi-period optimization studies. It
takes years as values, e.g. 2030, 2040, 2050.

The planning horizons are configured in the top-level ``planning_horizons``
config option:

.. code:: yaml

   planning_horizons: [2030, 2040, 2050]  # or single value: 2050

This wildcard appears in composed and solved network filenames:

- ``resources/{run}/networks/composed_{horizon}.nc``
- ``results/{run}/networks/solved_{horizon}.nc``

.. note::

   Legacy documentation used the ``{planning_horizons}`` wildcard name. Update
   custom scripts and configs to the ``{horizon}`` form.
