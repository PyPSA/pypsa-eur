.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
Migration Guide
##########################################

This guide covers the workflow refactor introduced in version XXX. The refactoring
aims to simplify and streamline the workflow by following the order
``base → simplified → clustered → composed → solved`` while only keeping the
``{run}`` and the ``{horizon}`` wildcard. All foresight modes, ``overnight``, ``myopic`` and ``perfect``,
follow this principle and only differ in potential additional loops over the rules: ``myopic`` iterates
``compose`` and ``solve`` for each horizon, ``perfect`` iterates over ``compose`` for each horizon and
then solves the network.

The previous wildcards ``{clusters}``, ``{opts}``, ``{sector_opts}`` have been removed. Their configuration is not handled by the ``config.yaml`` file only. The wildcard `{planning_horizons}` was renamed to `{horizon}`. 


Major Changes
=============

Rule Order Changes
------------------

The previous workflow used separate rules for network preparation: ``add_electricity``,
``prepare_network``, ``prepare_sector_network``, ``add_existing_baseyear``, ``add_brownfield``,
and ``prepare_perfect_foresight``. These have been consolidated into a single ``compose_network``
rule that handles all composition steps. Similarly, the separate
``solve_network`` (electricity-only) and ``solve_sector_network`` (sector-coupled) rules
are unified into a single ``solve_network`` rule.

The rule execution order depends on the foresight mode:

- **Overnight**: Composes and solves a single horizon:
  ``clustered.nc`` → ``composed_{horizon}.nc`` → ``solved_{horizon}.nc``

- **Myopic**: After composing and solving the first horizon, cycles over ``compose_network`` and ``solve_network`` 
  where `compose_network` reads the previous horizon's *solved* network to incorporate brownfield capacities.

- **Perfect foresight**: First composes all horizons sequentially, where each
  ``compose_network`` reads the previous horizon's *composed* network. 
  Then solves all horizons together in a single optimization.

All foresight modes and sector configurations (electricity-only and sector-coupled) support
the summary rules (``make_summary`` + plotting rules) and the ``all`` target rule.


File name changes
-----------------

With the removal of the wildcards, intermediate and final outputs were renamed as depicted in the following table.
The ``{run}`` prefix is omitted for simplicity.

.. list-table::
   :header-rows: 1

   * - Legacy
     - New
   * - ``resources/networks/base_s.nc``
     - ``resources/networks/simplified.nc``
   * - ``resources/networks/base_s_{clusters}.nc``
     - ``resources/networks/clustered.nc``
   * - ``resources/busmap_base_s_{clusters}.csv``
     - ``resources/busmap.csv``
   * - ``resources/powerplants_s_{clusters}.csv``
     - ``resources/powerplants.csv``
   * - ``resources/regions_onshore_base_s_{clusters}.geojson``
     - ``resources/onshore_regions.geojson``
   * - ``resources/regions_offshore_base_s_{clusters}.geojson``
     - ``resources/offshore_regions.geojson``
   * - ``resources/pop_layout_base_s_{clusters}.csv``
     - ``resources/pop_layout.csv``
   * - ``resources/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc``
     - ``resources/networks/composed_{horizon}.nc``
   * - ``results/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc``
     - ``results/networks/solved_{horizon}.nc``
   * - ``results/maps/base_s_*-costs-all_{planning_horizons}.pdf``
     - ``results/maps/power_network_{horizon}.pdf``
     
Configuration Changes
---------------------

- **Scenario management**: The ``scenario`` configuration section, which previously defined wildcard values for
  ``{clusters}``, ``{opts}``, and ``{sector_opts}``, has been removed. These settings are now
  specified directly in the relevant configuration sections (e.g., ``clustering.cluster_network.n_clusters``).
  To run multiple scenarios, use ``run.scenarios`` with separate scenario files
  (see :doc:`configuration` and :doc:`tutorial`).

- **Planning horizons**: The ``scenario.planning_horizons`` setting has moved to ``planning_horizons`` at the top level.
  This list defines the years for which the model is optimized and directly controls the ``{horizon}`` wildcard.

- **CO₂ handling**: The ``electricity.co2limit_enable``, ``electricity.co2limit``, and
  ``electricity.co2base`` settings are removed. CO₂ constraints are now exclusively
  configured via ``co2_budget``, which has a new structure with ``emissions_scope``,
  ``relative``, ``upper``, and ``lower`` keys.

- **Temporal resolution**: ``clustering.temporal.resolution_elec`` and ``resolution_sector``
  are merged into a single ``clustering.temporal.resolution`` setting. A new
  ``time_segmentation`` subsection provides ``enable``, ``resolution``, and ``segments``.

- **Line/link extensions**: ``lines.max_extension`` is renamed to ``lines.s_nom_max_extension``
  and ``links.max_extension`` to ``links.p_nom_max_extension``.

- **Sector toggle**: A new ``sector.enabled`` flag controls whether sector coupling is active.


**Before** (legacy)::

    scenario:
      clusters: [37]
      opts: [Co2L-24h]
      sector_opts: [Co2L0-1H-T-H-B-I]
      planning_horizons: [2030, 2040]

    clustering:
      temporal:
        resolution_elec: 24h
        resolution_sector: 1h

**Now** (new)::

    planning_horizons: [2030, 2040]

    clustering:
      cluster_network:
        n_clusters: 37
      temporal:
        resolution: 24h

    co2_budget:
      emissions_scope: CO2
      relative: true
      upper:
        2030: 0.45
        2040: 0.1

    sector:
      enabled: true

      
Snakemake Call Changes
----------------------

The streamlined workflow removes the sector-specific collection rules and replaces them with unified rules.

**Collection rules mapping:**

.. list-table::
   :header-rows: 1

   * - Legacy Rule
     - New Rule
   * - ``cluster_networks``
     - ``cluster_networks`` (targets ``clustered.nc`` instead of ``base_s_{clusters}.nc``)
   * - ``prepare_elec_networks``
     - ``compose_networks``
   * - ``prepare_sector_networks``
     - ``compose_networks``
   * - ``solve_elec_networks``
     - ``solve_networks``
   * - ``solve_sector_networks``
     - ``solve_networks``
   * - ``solve_sector_networks_perfect``
     - ``solve_networks``

Note all sector and foresight modes support the ``all`` rule. That is, you can **always** call::

    snakemake --configfile <configfile>

to trigger the whole workflow independently of the mode.


**Direct file targeting:**

When targeting individual files via Snakemake, use the new simplified paths:

.. list-table::
   :header-rows: 1

   * - Legacy
     - New
   * - ``snakemake resources/.../networks/base_s_37.nc``
     - ``snakemake resources/.../networks/clustered.nc``
   * - ``snakemake resources/.../networks/base_s_37_elec_Co2L-24h.nc``
     - ``snakemake resources/.../networks/composed_2030.nc``
   * - ``snakemake results/.../networks/base_s_37_elec_Co2L-24h.nc``
     - ``snakemake results/.../networks/solved_2030.nc``
   * - ``snakemake results/.../networks/base_s_37_Co2L-24h_Co2L0-1H-T-H-B-I_2030.nc``
     - ``snakemake results/.../networks/solved_2030.nc``

Note that the ``**config["scenario"]`` expansion pattern is no longer used. Collection rules now explicitly use ``horizon=config["planning_horizons"]`` for the ``{horizon}`` wildcard.

      
.. _fork-migration:

Migrating Forked Repositories
=============================

This section guides developers of forked PyPSA-Eur repositories through the migration process.
The streamlined workflow consolidates multiple Snakemake rules and restructures Python scripts,
which will likely cause merge conflicts when updating forks.


Summary of Structural Changes
-----------------------------

**Removed Snakemake rule files and their rules:**

- ``rules/solve_electricity.smk``: ``solve_network``, ``solve_operations_network``
- ``rules/solve_myopic.smk``: ``add_existing_baseyear``, ``add_brownfield``, ``solve_sector_network_myopic``
- ``rules/solve_overnight.smk``: ``solve_sector_network``
- ``rules/solve_perfect.smk``: ``add_existing_baseyear``, ``prepare_perfect_foresight``, ``solve_sector_network_perfect``, ``make_summary_perfect``

**Removed Snakemake rules** (from existing files):

- ``add_electricity``
- ``prepare_network``
- ``prepare_sector_network``

**New Snakemake rule files:**

- ``rules/compose.smk`` - Composition rule ``compose_network`` for all sector and foresight modes
- ``rules/solve.smk`` - Solve rule ``solve_network`` for all sector and foresight modes

**Script changes:**

The Python scripts ``add_electricity.py``, ``add_existing_baseyear.py``, ``add_brownfield.py``,
``prepare_network.py``, ``prepare_sector_network.py``, and ``prepare_perfect_foresight.py`` are
retained but refactored. Their main execution blocks are replaced by ``main()`` functions that
are imported and called by ``compose_network.py``.


Resolving Merge Conflicts in Python Scripts
-------------------------------------------

When merging the new structure into your fork, you will likely encounter conflicts in the
``add_*`` and ``prepare_*`` scripts. The key change is that the ``if __name__ == "__main__"``
block has been replaced by a ``main()`` function, looking like this:

.. code-block:: python

   def main(
       n: pypsa.Network,
       inputs,
       params,
       costs: pd.DataFrame,
   ) -> None:

(Some functions have additional parameters like ``nyears`` or ``current_horizon`` depending
on their requirements.)

In these main functions the previous ``snakemake.`` directives are replaced by the ``inputs``  and ``params`` arguments. 
To port your changes to the previous main section, update it to take into account the new accessor as shown in the following table  

.. list-table::
   :header-rows: 1

   * - Legacy Pattern
     - New Pattern
   * - ``snakemake.input.network``
     - ``inputs.network`` or ``inputs["network"]``
   * - ``snakemake.input["powerplants"]``
     - ``inputs["powerplants"]``
   * - ``snakemake.params.costs``
     - ``params.costs``
   * - ``snakemake.config["sector"]["enabled"]``
     - ``params.sector["enabled"]``

     
Note that if your custom code accesses configuration via ``snakemake.config["key"]``, you must add
the corresponding entry to the ``params`` section of the ``compose_network`` rule in
``rules/compose.smk``:

.. code-block:: python

   rule compose_network:
       params:
           your_custom_param=config_provider("path", "to", "config"),


Resolving Merge Conflicts in Snakemake Rules
--------------------------------------------

If your fork modified any of the removed rules, keep the rules deleted and migrate the changes to ``rules/compose.smk``:

1. Move custom ``input:`` entries to the ``get_compose_inputs()`` function. 
2. Move custom ``params:`` entries to the ``compose_network`` rule's ``params:`` section
3. Move custom logic from the associated script into its ``main()`` function

**File naming changes:**

Update any hardcoded file references to use the new naming convention. Remove all wildcard references except for ``{planning_horizons}`` which should be replaced with ``{horizon}``. See the examples in the following table:

.. list-table::
   :header-rows: 1

   * - Legacy
     - New
   * - ``profile_{clusters}_{technology}.nc``
     - ``profile_{technology}.nc``
   * - ``powerplants_s_{clusters}.csv``
     - ``powerplants.csv``
   * - ``costs_{planning_horizons}_processed.csv``
     - ``costs_{horizon}_processed.csv``
   * - ``regions_onshore_base_s_{clusters}.geojson``
     - ``onshore_regions.geojson``
   * - ``regions_offshore_base_s_{clusters}.geojson``
     - ``offshore_regions.geojson``

Note the special case of ``regions_onshore_*.geojson``/``regions_offshore_*.geojson`` which are renamed to ``onshore_regions``/``offshore_regions`` (reverted order of words).
