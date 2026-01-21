.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
Migration Guide
##########################################

This guide summarises the changes introduced with the PR1838 workflow
refactor and serves as the canonical reference for upgrading custom rules,
scripts, and configs.

Overview
========

- The Snakemake pipeline now runs ``base → simplified → clustered → composed →
  solved`` with one ``compose_network``/``solve_network`` pair instead of the
  scattered ``add_*``/``prepare_*`` sequence.
- File names now encode scenario information via configuration settings and the
  ``{horizon}`` wildcard rather than embedding ``{clusters}``, ``{opts}``, and
  ``{planning_horizons}`` inside every artefact.
- Configuration options such as ``planning_horizons`` and CO₂ budgets moved to
  the top level, so scenario sweeps reference config files directly instead of
  wildcard combinations.


Workflow changes
================

The table below maps legacy targets to their equivalents in the refactored
workflow.

.. list-table::
   :header-rows: 1

   * - Legacy target
     - New target
     - Notes
   * - ``networks/base_s.nc``
     - ``networks/simplified.nc``
     - Simplification outputs no longer embed cluster counts.
   * - ``busmap_base_s_{clusters}.csv`` / ``linemap_base_s_{clusters}.csv``
     - ``busmap.csv`` / ``linemap.csv``
     - Cluster counts come from ``clustering.cluster_network.n_clusters``.
   * - ``powerplants_s_{clusters}.csv``
     - ``powerplants_s.csv``
     - Produced once per clustered network.
   * - ``networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc``
     - ``networks/composed_{horizon}.nc``
     - Single entry point handled by ``scripts/compose_network.py``.
   * - ``RESULTS/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc``
     - ``RESULTS/networks/solved_{horizon}.nc``
     - Solver logs live under ``results/{run}/logs/solve_network``.
   * - ``RESULTS/maps/base_s_{clusters}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf``
     - ``RESULTS/maps/power_network_{horizon}.pdf`` (and matching carrier maps)
     - Same naming for overnight, myopic, and perfect modes.

``compose_network`` now orchestrates all per-horizon assembly steps, from
merging simplified assets to applying brownfield capacities and warm-start data.
Legacy Snakemake workflows split these responsibilities across multiple
``add_*`` and ``prepare_*`` rules, so custom injections should extend the
dedicated sections inside ``scripts/compose_network.py`` instead of re-adding
bespoke steps.

Configuration changes
=====================

1. Define ``planning_horizons`` at the top level. Provide a single value for
   overnight studies or a list for myopic/perfect runs. The ``{horizon}``
   wildcard is derived from this list.
2. CO₂ budgets share a unified schema with ``values`` (``fraction`` vs
   ``absolute``) and per-period ``upper``/``lower`` entries. Optional ``total``
   + ``distribution`` fields distribute an aggregate cap across horizons.
3. Transmission expansion limits use the ``transmission_limit: <metric><cap>``
   syntax (e.g., ``vopt`` or ``c1.25``).
4. ``existing_capacities.enabled`` gates brownfield injections; leave it
   ``false`` for greenfield runs.
5. Snakemake rules should access configuration through
   ``config_provider(...)`` or ``get_config(w)`` to ensure merged scenario files
   and wildcards are respected. Prior releases collected many of these knobs via
   the ``scenario`` wildcard block, which is now ignored; encode them directly
   in configuration files and rely on ``run.scenarios`` for sweeps.
6. Temporal resolution settings ``clustering.temporal.resolution_elec`` and
   ``clustering.temporal.resolution_sector`` have been unified into a single
   ``clustering.temporal.resolution`` setting.

Foresight modes
=======================

- **Overnight**: Require a single horizon; ``compose_network`` never looks for
  previous outputs.
- **Myopic**: ``compose_network`` imports
  ``results/{run}/networks/solved_{previous}.nc`` via the helper
  ``solved_previous_horizon``. Ensure horizons are sorted ascendingly.
- **Perfect**: ``compose_network`` uses ``resources/{run}/networks/composed_{previous}.nc`` as  
  the brownfield baseline to build the full multi-period optimisation.


