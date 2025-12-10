.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _foresight:

#####################
Foresight Options
#####################

Planning horizons and warm starts
=================================

All foresight modes share the same ``planning_horizons`` list at the top of the
configuration file. The workflow iterates over this list and composes/solves a
network for each year while persisting two canonical artefacts:

- ``resources/{run}/networks/composed_{horizon}.nc`` — output of
  :mod:`scripts/compose_network`
- ``results/{run}/networks/solved_{horizon}.nc`` — output of
  :mod:`scripts/solve_network`

Warm-start behaviour depends on the foresight mode:

- ``overnight`` expects a single value and does **not** reuse previous years.
- ``myopic`` requires at least two horizons and feeds the solved network from
  the previous year ``RESULTS/networks/solved_{prev}.nc`` into the next
  ``compose_network`` call.
- ``perfect`` also iterates sequentially but consumes
  ``networks/composed_{prev}.nc`` to build the full multi-period optimisation.

When in doubt, inspect :mod:`rules/compose.smk` to see how the helpers
``solved_previous_horizon`` and ``compose_previous_horizon`` assemble these
paths.

.. _overnight:

Overnight (greenfield) scenarios
================================

The default is to calculate a rebuilding of the energy system to meet demand, a
so-called overnight or greenfield approach.

In this case, the ``planning_horizons`` parameter specifies the reference year
for exogenously given transition paths (e.g. the level of steel recycling). It
does not affect the year for cost and technology assumptions, which is set
separately in the config.

.. code:: yaml

  planning_horizons: 2050

  costs:
    year: 2030

For running overnight scenarios, use in the ``config/config.yaml``:

.. code:: yaml

  foresight: overnight

.. _perfect:

Perfect foresight scenarios
===========================

.. warning::

  Perfect foresight is currently implemented as an experimental test version.

For running perfect foresight scenarios, you can adjust the
 ``config/config.perfect.yaml``:

.. code:: yaml

  foresight: perfect


.. _myopic:

Myopic foresight scenarios
=============================

The myopic code can be used to investigate progressive changes in a network, for
instance, those taking place throughout a transition path. The capacities
installed in a certain time step are maintained in the network until their
operational lifetime expires.

The myopic approach was initially developed and used in the paper `Early
decarbonisation of the European Energy system pays off (2020)
<https://www.nature.com/articles/s41467-020-20015-4>`__ and later further
extended in `Speed of technological transformations required in Europe to
achieve different climate goals (2022)
<https://doi.org/10.1016/j.joule.2022.04.016>`__. The current implementation
complies with the PyPSA-Eur-Sec standard working flow and is compatible with
using the higher resolution electricity transmission model `PyPSA-Eur
<https://github.com/PyPSA/pypsa-eur>`__ rather than a one-node-per-country
model.

The current code applies the myopic approach to generators, storage technologies
and links in the power sector. It furthermore applies it to the space and water
heating sector (e.g., the share of district heating and reduced space heat
demand), industry processes (e.g., steel, direct reduced iron, and aluminum
production via primary route), the share of fuel cell and battery electric
vehicles in land transport, and the hydrogen share in shipping (see
:doc:`supply_demand` for further information).

The following subjects within the land transport and biomass currently do not
evolve with the myopic approach:

- The percentage of electric vehicles that allow demand-side management and
  vehicle-to-grid services.

- The annual biomass potential (default year and scenario for which potential is
  taken is 2030, as defined in config)

.. literalinclude:: ../config/test/config.myopic.yaml
   :language: yaml
   :start-at: biomass:
   :end-at: year:


Configuration
--------------

For running myopic foresight transition scenarios, set in ``config/config.yaml``:

.. code:: yaml

  foresight: myopic

The following options included in the ``config/config.yaml`` file  are relevant for the
myopic code.

The ``{horizon}`` wildcard indicates the year in which the network is optimized.
For a myopic optimization, this is equivalent to the investment year. The set
of values is pulled directly from the top-level ``planning_horizons`` list. To
set the investment years which are sequentially simulated for the myopic
investment planning, select for example:

.. literalinclude:: ../config/test/config.myopic.yaml
   :language: yaml
   :start-at:   planning_horizons:
   :end-before: countries:


**existing capacities**

Grouping years indicates the bins limits for grouping the existing capacities of
different technologies. Note that separate bins are defined for the power and
heating plants due to different data sources.

``grouping_years_power: [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020,
2025, 2030]``

``grouping_years_heat: [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2019]``





**threshold capacity**

If for a technology, node, and grouping bin, the capacity is lower than
threshold_capacity, it is ignored.

``threshold_capacity: 10``




**conventional carriers**

Conventional carriers indicate carriers used in the existing conventional
technologies.

    conventional_carriers:

    \- lignite

    \- coal

    \- oil

    \- uranium




Options
--------------

Carbon Budget Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The carbon budget for the entire transition path can be specified in the
``co2_budget`` section in ``config/config.yaml`` using a unified configuration
structure that works consistently across all three foresight modes (overnight,
myopic, perfect).

**Unified Configuration Structure:**

.. code:: yaml

  co2_budget:
    emissions_scope: All greenhouse gases - (CO2 equivalent)
    values: absolute  # How to interpret per-period values: "absolute" (Gt CO2/year) or "fraction" (% of 1990)
    # Upper emissions constraints
    upper:
      enable: true

      # Option 1: Per-period hard caps (Gt CO2/year)
      2030: 2.04  # 55% reduction by 2030 (Fit for 55)
      2040: 0.45  # 90% by 2040
      2050: 0.00  # climate-neutral by 2050

      # Option 2: Total cumulative budget with automatic distribution
      total: 40.0  # Total budget in Gt CO2 over planning horizon
      distribution: ex0  # Distribution mode (see below)

    # Lower emissions constraints (optional)
    lower:
      enable: true
      2030: 0.5  # Minimum emissions floor for 2030 (Gt CO2/year)

**Configuration Options:**

1. **Per-period caps**: Specify absolute emissions caps (Gt CO₂/year) for specific years.
   When both a per-period cap and a distributed budget exist for the same horizon, the
   minimum (more restrictive) value is used to ensure budget constraints aren't loosened.

2. **Total budget with distribution**: Specify a cumulative total budget and how it should
   be distributed across planning horizons.

   - For **perfect foresight**: If ``distribution`` is null, ``total`` is applied as a single
     cumulative constraint over all periods. If ``distribution`` is specified, the budget is
     distributed across periods using the specified mode (same as myopic/overnight).

   - For **myopic/overnight**: ``total`` is split across periods using the specified
     ``distribution`` mode.

**Distribution Modes:**

- ``ex0``: Exponential decay with initial growth rate r=0 (commonly used)
- ``ex1``: Exponential decay with initial growth rate r=1
- ``be1``: Beta decay with parameter 1
- ``be2``: Beta decay with parameter 2
- ``linear``: Linear decrease from current emissions to zero
- ``equal``: Equal distribution across all periods

The exponential decay follows:

.. math::
  e(t) = e_0 (1+ (r+m)t) e^{-mt}

where r is the initial linear growth rate, e_0 is the initial emission level (typically
from 2019), and the decay parameter m is determined by imposing that the integral of
the path equals the specified budget.

**Example Configurations:**

*Example 1: Per-period caps with fraction values (default)*

.. code:: yaml

  co2_budget:
    values: fraction  # Values as % of 1990 baseline (default)
    upper:
      enable: true
      2030: 0.450  # 45% of 1990 emissions (55% reduction)
      2040: 0.100  # 10% of 1990 emissions (90% reduction)
      2050: 0.000  # Climate neutral

*Example 2: Per-period caps with absolute values*

.. code:: yaml

  co2_budget:
    values: absolute  # Values in Gt CO2/year
    upper:
      enable: true
      2030: 2.04  # Gt CO2/year
      2040: 0.45
      2050: 0.00

*Example 3: Total budget for perfect foresight*

.. code:: yaml

  foresight: perfect
  planning_horizons: [2030, 2040, 2050]

  co2_budget:
    upper:
      enable: true
      total: 40.0  # Cumulative constraint over all periods
      distribution: null  # Ignored for perfect foresight

    lower:
      enable: true
      2030: 0.5  # Prevent too-fast decarbonization

*Example 4: Total budget for myopic/overnight with distribution*

.. code:: yaml

  foresight: myopic
  planning_horizons: [2030, 2040, 2050]

  co2_budget:
    upper:
      enable: true
      total: 40.0
      distribution: ex0  # Required for myopic/overnight

*Example 5: Mixed approach (caps + total budget)*

.. code:: yaml

  co2_budget:
    upper:
      enable: true
      2030: 2.04  # Per-period cap for 2030
      total: 40.0  # Applied to other periods
      distribution: ex0

*Example 6: Selective year constraints*

.. code:: yaml

  planning_horizons: [2030, 2040, 2050]

  co2_budget:
    values: fraction
    upper:
      enable: true
      2030: 0.450  # Only constrain 2030
      2050: 0.000  # And 2050
      # 2040 not specified - no upper constraint for 2040
    lower:
      enable: true
      2030: 0.200  # Prevent too-fast decarbonization in 2030 only
      # 2040 and 2050 not specified - no lower constraints

**Constraint Application Logic:**

The model applies CO₂ constraints for each planning horizon following this priority:

1. **Per-period caps take precedence**: If a per-period cap is specified for a horizon, it's
   always applied. When both a per-period cap and a distributed budget value exist for the
   same year, the model uses the **minimum** (more restrictive) value.

2. **Distributed budget as fallback**: If no per-period cap exists for a horizon, the model
   uses the distributed budget value (if total budget with distribution is configured).

3. **Perfect foresight cumulative constraint**: For perfect foresight only, if ``total``
   is specified without ``distribution`` (null), an additional cumulative constraint is
   applied over all periods at the final horizon.

4. **Lower bounds**: When ``lower.enable: true``, lower bounds are applied only for planning
   horizons where values are explicitly specified. Unspecified horizons have no lower constraint.
   Note: Lower-only constraints (without upper) are not supported; at least one upper constraint
   must be specified.

**Relationship to Published Research:**

The paper `Speed of technological transformations required in Europe to achieve
different climate goals (2022) <https://doi.org/10.1016/j.joule.2022.04.016>`__
defines CO₂ budgets corresponding to global temperature increases (1.5°C – 2°C).
Global carbon budgets are converted to European budgets assuming equal per-capita
distribution (6.43% share for Europe). These can be specified using the ``total``
and ``distribution`` parameters. For example, a 1.5°C target corresponds to
approximately 25.7 Gt CO₂ for Europe, which can be configured as:

.. code:: yaml

  co2_budget:
    upper:
      enable: true
      total: 25.7  # Gt CO2 for 1.5°C target
      distribution: ex0

See Supplemental Note S1 of the paper for detailed derivations of the distribution
functions.

**Emissions Scope:**

The ``emissions_scope`` parameter determines which greenhouse gas(es) are accounted for
when calculating historical baseline emissions (e.g., 1990 emissions) and applying budget
constraints. This parameter corresponds to the ``Pollutant_name`` field in the EEA UNFCCC
emissions database.

Available options:

- ``All greenhouse gases - (CO2 equivalent)`` - All greenhouse gases in CO₂-equivalent,
  including CO₂, CH₄, N₂O, and fluorinated gases (HFCs, PFCs, SF₆, NF₃) (**default**)
- ``CO2`` - Carbon dioxide only
- ``CH4`` - Methane emissions only
- ``N2O`` - Nitrous oxide emissions only

The choice of emissions scope affects:

1. **Baseline calculations**: Historical emissions used to calculate budget fractions
   (e.g., "55% reduction from 1990" requires knowing 1990 emissions)
2. **Budget conversion**: When using ``cbXX`` format, the total budget is scaled based
   on the selected emissions scope
3. **Constraint application**: The optimization constraints in the network model apply
   to the selected emissions scope

.. note::
   The default uses all greenhouse gases in CO₂-equivalent to align with EU climate policy
   targets (e.g., Fit for 55: 55% GHG reduction by 2030, climate neutrality by 2050). This
   ensures budget calculations use the same baseline as official climate targets. While the
   energy system model primarily tracks and optimizes CO₂ emissions from energy use, using
   all GHGs for baseline calculations provides proper context for the model's emission
   constraints relative to economy-wide climate goals. Non-energy GHG emissions (e.g., from
   agriculture, industrial processes) are typically handled as exogenous assumptions or
   boundary conditions.


General myopic code structure
---------------------------------

The myopic workflow iterates through ``planning_horizons`` inside a single
:mod:`scripts/compose_network` entry point. For each horizon the rule performs:

1. Load ``networks/clustered.nc`` together with all derived assets (busmaps,
   demand profiles, sectoral inputs).
2. If ``existing_capacities.enabled`` is ``true``, merge the brownfield assets
   built before the base year via ``add_existing_capacities`` inside
   ``compose_network.py``.
3. When ``w.horizon`` is not the first element of ``planning_horizons``, import
   ``RESULTS/networks/solved_{prev}.nc`` (myopic) or
   ``networks/composed_{prev}.nc`` (perfect) as the warm-start baseline.
4. Persist the assembled network as ``resources/{run}/networks/composed_{horizon}.nc``.

After composition, :mod:`scripts/solve_network` optimises every horizon to
produce ``results/{run}/networks/solved_{horizon}.nc``. Downstream plotting and
reporting stages such as :mod:`scripts/make_summary` consume these solved files
directly, so all foresight modes share the same results directory structure.

Rule overview
--------------

- :mod:`compose_network`

  Performs the per-horizon assembly of simplified assets, existing capacities,
  and warm-start seeding within a single entry point. Custom extensions should
  hook into the clearly marked sections of ``scripts/compose_network.py``.

  .. note::

     Earlier workflows used dedicated ``add_*`` and ``prepare_*`` rules. Keep
     bespoke injections inside ``compose_network.py`` to stay aligned with the
     streamlined pipeline.

- :mod:`solve_network`

  Applies the configured solver to each composed network and writes
  ``results/{run}/networks/solved_{horizon}.nc`` alongside solver logs stored
  in ``results/{run}/logs/solve_network``.

- :mod:`make_summary`

  Converts the solved networks into CSV/plot outputs for downstream dashboards.

  .. note::

     Legacy helper scripts such as ``make_summary_perfect.py`` distinguished
     between foresight modes. The refactored workflow routes every summary
     through ``scripts/make_summary.py``.
