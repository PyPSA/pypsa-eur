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


Migrating forked repositories
=============================

If you maintain a fork of PyPSA-Eur with custom scripts or rules, the
streamlined workflow consolidates previously scattered logic into
``scripts/compose_network.py``. This section explains how to port your
customizations.

Customization points mapping
----------------------------

The table below shows where legacy customization points now reside:

.. list-table::
   :header-rows: 1
   :widths: 35 35 30

   * - Legacy location
     - New location
     - Section marker in compose_network.py
   * - ``add_electricity.py`` main section
     - ``add_electricity_components()`` function
     - ``ELECTRICITY COMPONENTS (from add_electricity.py)``
   * - ``prepare_sector_network.py`` main section
     - ``add_sector_components()`` function
     - ``SECTOR COMPONENTS (from prepare_sector_network.py)``
   * - ``add_existing_baseyear.py`` main section
     - ``add_existing_capacities()`` function
     - ``EXISTING CAPACITIES (from add_existing_baseyear.py)``
   * - ``add_brownfield.py`` main section
     - Inline brownfield block
     - ``BROWNFIELD FOR MYOPIC (from add_brownfield.py)``
   * - ``prepare_network.py`` main section
     - ``prepare_network_for_solving()`` function
     - ``NETWORK PREPARATION (from prepare_network.py)``
   * - ``prepare_perfect_foresight.py`` main section
     - ``concatenate_network_with_previous()`` + inline block
     - ``PERFECT FORESIGHT CONCATENATION (from prepare_perfect_foresight.py)``

Porting custom functions
------------------------

**Scenario 1: You added a function to a legacy script**

If you added a helper function (e.g., ``attach_custom_generators()``) to
``add_electricity.py``, keep the function in that file. The streamlined
workflow imports functions from legacy scripts, so your function remains
available::

    # In your fork's add_electricity.py, add your function:
    def attach_custom_generators(n, costs, params):
        """Add custom generator type."""
        ...

    # In compose_network.py, add the import and call:
    from scripts.add_electricity import (
        ...
        attach_custom_generators,  # Add your import
    )

    # Then call it in the ELECTRICITY COMPONENTS section:
    attach_custom_generators(n, costs, params)

**Scenario 2: You modified the main execution flow**

If you inserted logic between existing steps (e.g., filtering buses after
load attachment), locate the corresponding section in ``compose_network.py``
and insert your code there. Each section is clearly marked with comments like::

    # ========== ELECTRICITY COMPONENTS (from add_electricity.py) ==========

**Scenario 3: You added a custom Snakemake rule**

If you had a custom rule (e.g., ``add_custom_data``) that ran between
``add_electricity`` and ``prepare_network``, you have two options:

1. **Integrate into compose_network.py**: Add your logic as a function call
   in the appropriate section. This is preferred for logic that modifies
   the network.

2. **Keep as separate rule**: If your rule produces independent data files,
   keep it separate and add its outputs to ``get_compose_inputs()`` in
   ``rules/compose.smk``.

Adding custom inputs to compose_network
---------------------------------------

To add custom data files as inputs to the compose rule, extend
``get_compose_inputs()`` in ``rules/compose.smk``::

    def get_compose_inputs(w):
        cfg = get_config(w)
        ...
        inputs = dict(
            ...existing inputs...
        )

        # Add your custom inputs
        if cfg.get("my_custom_feature", {}).get("enabled", False):
            inputs["custom_data"] = resources("my_custom_data.csv")
            inputs["custom_profiles"] = resources("my_custom_profiles_{horizon}.nc")

        return inputs

Then access them in ``compose_network.py`` via ``snakemake.input.custom_data``.

Example: Porting a custom carrier
---------------------------------

Suppose your fork adds a "green ammonia" carrier in ``prepare_sector_network.py``.
Here's how to port it:

1. **Keep the function in the legacy script**::

       # scripts/prepare_sector_network.py
       def add_green_ammonia(n, costs, spatial, options):
           """Add green ammonia production and storage."""
           ...

2. **Add the import in compose_network.py**::

       from scripts.prepare_sector_network import (
           ...
           add_green_ammonia,
       )

3. **Call it in the SECTOR COMPONENTS section**::

       # In compose_network.py, after existing sector components:
       if sector_opts.get("green_ammonia", False):
           add_green_ammonia(n, costs, spatial, sector_opts)

4. **Add any required inputs to get_compose_inputs()**::

       # In rules/compose.smk
       if cfg["sector"].get("green_ammonia", False):
           inputs["green_ammonia_potentials"] = resources("green_ammonia_potentials.csv")

Best practices for fork maintenance
-----------------------------------

1. **Minimize changes to compose_network.py**: Keep custom functions in the
   legacy script files and only add imports and calls in ``compose_network.py``.
   This reduces merge conflicts when pulling upstream changes.

2. **Use config flags for custom features**: Gate your customizations with
   configuration options so they can be disabled when testing against upstream.

3. **Document your integration points**: Add comments indicating where your
   custom code integrates, referencing any related issues or documentation.

4. **Test incrementally**: After merging upstream changes, run the test suite
   with your customizations disabled first, then enable them one by one.

5. **Watch for function signature changes**: When upstream modifies a function
   you import, check if the signature changed. The deprecation of main sections
   means function signatures are now the stable API.


