.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _foresight:

#####################
Foresight Options
#####################

Planning horizons
=================================

All foresight modes share the same ``planning_horizons`` list at the top of the
configuration file. The workflow iterates over this list and composes/solves a
network for each year:

- ``resources/{run}/networks/composed_{horizon}.nc`` — output of
  :mod:`compose_network`
- ``results/{run}/networks/solved_{horizon}.nc`` — output of
  :mod:`scripts/solve_network`

Iteration behaviour depends on the foresight mode:

- ``overnight`` expects a single value and does **not** reuse previous years.
- ``myopic`` requires at least two horizons and feeds the solved network from
  the previous year ``RESULTS/networks/solved_{prev}.nc`` into the next
  ``compose_network`` call.
- ``perfect`` also iterates sequentially but consumes
  ``networks/composed_{prev}.nc`` to build the full multi-period optimisation, after which the complete network is solved.



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

The ``{horizon}`` wildcard indicates the year for which the network is optimized.
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
``co2_budget`` section in ``config/config.yaml`` for all three foresight modes (overnight, myopic, perfect).

  co2_budget:
    emissions_scope: CO2
    relative: true  # true = fraction of 1990 baseline, false = absolute (Gt CO2/year)
    upper: null     # null | scalar | {year: value}
    lower: null     # null | scalar | {year: value}

**Constraint semantics:**

- ``upper``/``lower`` as a **dict**: apply constraints only for explicitly listed years
  (omitted years are unconstrained). In perfect foresight, these constraints are attached
  to the corresponding investment period. Use ``<year>: null`` to explicitly clear an
  inherited default bound for a specific year.
- ``upper``/``lower`` as a **scalar**: add a single constraint as a budget across all investment periods.

  - **overnight**: one CO₂ constraint for the solved model.
  - **myopic**: one CO₂ constraint for each solved horizon.
  - **perfect**: one total CO₂ constraint across all horizons (budget), added when composing the final horizon.  

**Example 1: Per-period caps with relative values**

.. code:: yaml

  co2_budget:
    relative: true
    upper:
      2030: 0.450  # 45% of 1990 emissions
      2040: 0.100
      2050: 0.000

**Example 2: Scalar cap (absolute values)**

.. code:: yaml

  co2_budget:
    relative: false
    upper: 2.04  # Gt CO2/year, same semantics across foresight modes (see above)

**Example 3: Selective year constraints**

.. code:: yaml

  planning_horizons: [2030, 2040, 2050]

  co2_budget:
    relative: false
    upper:
      2030: 2.0
      2050: 0.0
      # 2040 not specified -> no upper constraint for 2040
    lower:
      2030: 0.5

**Emissions scope:**

The ``emissions_scope`` parameter determines which greenhouse gas(es) are accounted for
when calculating the 1990 baseline for ``relative: true`` and applying the corresponding
budget constraints. This parameter corresponds to the ``Pollutant_name`` field in the EEA
UNFCCC emissions database.

Only available options currently is `CO2`, other options that could potentially be tracked 

- ``All greenhouse gases - (CO2 equivalent)`` - All greenhouse gases in CO₂-equivalent,
  including CO₂, CH₄, N₂O, and fluorinated gases (HFCs, PFCs, SF₆, NF₃)
- ``CH4`` - Methane emissions only
- ``N2O`` - Nitrous oxide emissions only

The choice of emissions scope affects:

1. **Baseline calculations**: Historical emissions used to calculate budget fractions
   (e.g., "55% reduction from 1990" requires knowing 1990 emissions)
2. **Constraint application**: The optimization constraints in the network model apply
   to the selected emissions scope

.. note::
   The default uses ``CO2`` emissions scope, tracking only direct energy-related CO₂ emissions.
   This aligns with the model's primary focus on optimizing energy system decarbonization.
   While EU climate policy targets (e.g., Fit for 55) include all greenhouse gases in CO₂-equivalent,
   the energy system model constraints apply specifically to energy-related CO₂. For analyses
   requiring broader GHG scope alignment with EU targets, the ``emissions_scope`` parameter
   can be set to ``All greenhouse gases - (CO2 equivalent)``. Non-energy GHG emissions
   (e.g., from agriculture, industrial processes) are typically handled as exogenous
   assumptions or boundary conditions.


General myopic code structure
---------------------------------

The myopic workflow iterates through ``planning_horizons`` inside a single
:mod:`scripts/compose_network` entry point. For each horizon the rule performs:

1. Load ``networks/clustered.nc`` together with all derived assets (busmaps,
   demand profiles, sectoral inputs).
2. If ``existing_capacities.enabled`` is ``true``, insert the brownfield assets  
   built before the base year via ``add_existing_capacities`` inside
   :mod:`compose_network`.  
3. When ``w.horizon`` is not the first element of ``planning_horizons``, import
   ``RESULTS/networks/solved_{prev}.nc`` (myopic) or
   ``networks/composed_{prev}.nc`` (perfect) as the solution from the previous planning horizon.  
4. Store the composed network as ``resources/{run}/networks/composed_{horizon}.nc``.  

After composition, :mod:`solve_network` optimises every horizon to  
produce ``results/{run}/networks/solved_{horizon}.nc``. Downstream plotting and  
reporting stages such as :mod:`make_summary` consume these solved files  
directly.  

Rule overview
--------------

- :mod:`compose_network`

  Performs the per-horizon assembly of simplified assets, existing capacities,
  and investments from solved previous planning horizons.  
  hook into the clearly marked sections of ``scripts/compose_network.py``.

  .. note::

    Earlier workflows used several add_* and prepare_* rules, which are now unified in a single :mod:`compose_network` rule."

- :mod:`solve_network`

  Solves each composed network and writes and writes 
  ``results/{run}/networks/solved_{horizon}.nc`` alongside solver logs stored
  in ``results/{run}/logs/solve_network``.

- :mod:`make_summary`

  Converts the solved networks into CSVs.

  .. note::

    Legacy helper scripts such as ``make_summary_perfect.py`` distinguished
    between foresight modes. The refactored workflow routes every summary
    through ``scripts/make_summary.py``.
