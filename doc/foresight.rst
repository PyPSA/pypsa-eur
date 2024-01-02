..
  SPDX-FileCopyrightText: 2021-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _foresight:

#####################
Foresight Options
#####################

.. _overnight:

Overnight (greenfield) scenarios
================================

The default is to calculate a rebuilding of the energy system to meet demand, a so-called overnight or greenfield approach.

In this case, the ``planning_horizons`` parameter specifies the reference year for exogenously given transition paths (e.g. the level of steel recycling).
It does not affect the year for cost and technology assumptions, which is set separately in the config.

.. code:: yaml

  scenario:
    planning_horizons:
    - 2050

  costs:
    year: 2030

For running overnight scenarios, use in the ``config/config.yaml``:

.. code:: yaml

  foresight: overnight

.. _perfect:

Perfect foresight scenarios
===========================

.. warning::

  Perfect foresight is currently implemented as a first test version.

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

The ``{planning_horizons}`` wildcard indicates the year in which the network is
optimized. For a myopic optimization, this is equivalent to the investment year.
To set the investment years which are sequentially simulated for the myopic
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

The total carbon budget for the entire transition path can be indicated in the
`sector_opts
<https://github.com/PyPSA/pypsa-eur-sec/blob/f13902510010b734c510c38c4cae99356f683058/config.default.yaml#L25>`_
in ``config/config.yaml``. The carbon budget can be split among the
``planning_horizons`` following an exponential or beta decay. E.g. ``'cb40ex0'``
splits a carbon budget equal to 40 Gt :math:`_{CO_2}` following an exponential
decay whose initial linear growth rate r is zero. They can also follow some
user-specified path, if defined `here
<https://github.com/PyPSA/pypsa-eur-sec/blob/413254e241fb37f55b41caba7264644805ad8e97/config.default.yaml#L56>`_.
The paper `Speed of technological transformations required in Europe to achieve
different climate goals (2022) <https://doi.org/10.1016/j.joule.2022.04.016>`__
defines CO_2 budgets corresponding to global temperature increases (1.5C – 2C)
as response to the emissions. Here, global carbon budgets are converted to
European budgets assuming equal-per capita distribution which translates into a
6.43% share for Europe. The carbon budgets are in this paper distributed
throughout the transition paths assuming an exponential decay. Emissions e(t) in
every year t are limited by

.. math::
  e(t) = e_0 (1+ (r+m)t) e^{-mt}

where r is the initial linear growth rate, which here is assumed to be r=0, and
the decay parameter m is determined by imposing the integral of the path to be
equal to the budget for Europe. Following this approach, the CO_2 budget is
defined. Following the same approach as in this paper, add the following to the
``scenario.sector_opts`` E.g.  ``-cb25.7ex0`` (1.5C increase) Or ``cb73.9ex0``
(2C increase). See details in Supplemental Note S1 `Speed of technological
transformations required in Europe to achieve different climate goals (2022)
<https://doi.org/10.1016/j.joule.2022.04.016>`__.


General myopic code structure
---------------------------------

The myopic code solves the network for the time steps included in
``planning_horizons`` in a recursive loop, so that:

1. The existing capacities (those installed before the base year are added as
   fixed capacities with p_nom=value, p_nom_extendable=False). E.g. for
   baseyear=2020, capacities installed before 2020 are added. In addition, the
   network comprises additional generator, storage, and link capacities with
   p_nom_extendable=True. The non-solved network is saved in
   ``results/run_name/networks/prenetworks-brownfield``.

The base year is the first element in ``planning_horizons``. Step 1 is
implemented with the rule add_baseyear for the base year and with the rule
add_brownfield for the remaining planning_horizons.

2. The 2020 network is optimized. The solved network is saved in
   ``results/run_name/networks/postnetworks``

3. For the next planning horizon, e.g. 2030, the capacities from a previous time
   step are added if they are still in operation (i.e., if they fulfil planning
   horizon <= commissioned year + lifetime). In addition, the network comprises
   additional generator, storage, and link capacities with
   p_nom_extendable=True. The non-solved network is saved in
   ``results/run_name/networks/prenetworks-brownfield``.

Steps 2 and 3 are solved recursively for all the planning_horizons included in
``config/config.yaml``.

Rule overview
--------------

- rule add_existing baseyear

  The rule add_existing_baseyear loads the network in
  ‘results/run_name/networks/prenetworks’ and performs the following operations:

  1. Add the conventional, wind and solar power generators that were installed
     before the base year.

  2. Add the heating capacities that were installed before the base year.

  The existing conventional generators are retrieved from the `powerplants.csv
  file
  <https://pypsa-eur.readthedocs.io/en/latest/preparation/build_powerplants.html?highlight=powerplants>`__
  generated by pypsa-eur which, in turn, is based on the `powerplantmatching
  <https://github.com/FRESNA/powerplantmatching>`__ database.

  Existing wind and solar capacities are retrieved from `IRENA annual statistics
  <https://www.irena.org/Statistics/Download-Data>`__ and distributed among the
  nodes in a country proportional to capacity factor. (This will be updated to
  include capacity distributions closer to reality.)

  Existing heating capacities are retrieved from the report `Mapping and
  analyses of the current and future (2020 - 2030) heating/cooling fuel
  deployment (fossil/renewables)
  <https://ec.europa.eu/energy/studies/mapping-and-analyses-current-and-future-2020-2030-heatingcooling-fuel-deployment_en?redir=1>`__.

  The heating capacities are assumed to have a lifetime indicated by the
  parameter lifetime in the configuration file, e.g 25 years. They are assumed
  to be decommissioned linearly starting on the base year, e.g., from 2020 to
  2045.

  Then, the resulting network is saved in
  ``results/run_name/networks/prenetworks-brownfield``.

- rule add_brownfield

  The rule add_brownfield loads the network in
  ``results/run_name/networks/prenetworks`` and performs the following
  operation:

  1. Read the capacities optimized in the previous time step and add them to the
     network if they are still in operation (i.e., if they fulfill planning
     horizon < commissioned year + lifetime)

  Then, the resulting network is saved in
  ``results/run_name/networks/prenetworks_brownfield``.
