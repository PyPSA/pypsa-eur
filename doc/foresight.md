<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Foresight Options {#foresight}

## Overnight (greenfield) scenarios {#overnight}

The default is to calculate a rebuilding of the energy system to meet demand, a so-called overnight or greenfield approach.

In this case, the `planning_horizons` parameter specifies the reference year for exogenously given transition paths (e.g. the level of steel recycling).
It does not affect the year for cost and technology assumptions, which is set separately in the config.

```yaml
scenario:
  planning_horizons:
  - 2050

costs:
  year: 2030
```

For running overnight scenarios, use in the `config/config.yaml`:

```yaml
foresight: overnight
```

## Perfect foresight scenarios {#perfect}

!!! warning

    Perfect foresight is currently implemented as an experimental test version.

For running perfect foresight scenarios, you can adjust the
 `config/config.perfect.yaml`:

```yaml
foresight: perfect
```

## Myopic foresight scenarios {#myopic}

The myopic code can be used to investigate progressive changes in a network, for
instance, those taking place throughout a transition path. The capacities
installed in a certain time step are maintained in the network until their
operational lifetime expires.

The myopic approach was initially developed and used in the paper [Early
decarbonisation of the European Energy system pays off (2020)](https://www.nature.com/articles/s41467-020-20015-4) and later further
extended in [Speed of technological transformations required in Europe to
achieve different climate goals (2022)](https://doi.org/10.1016/j.joule.2022.04.016). The current implementation
complies with the PyPSA-Eur-Sec standard working flow and is compatible with
using the higher resolution electricity transmission model [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) rather than a one-node-per-country
model.

The current code applies the myopic approach to generators, storage technologies
and links in the power sector. It furthermore applies it to the space and water
heating sector (e.g., the share of district heating and reduced space heat
demand), industry processes (e.g., steel, direct reduced iron, and aluminum
production via primary route), the share of fuel cell and battery electric
vehicles in land transport, and the hydrogen share in shipping (see
[Supply and Demand](supply_demand.md) for further information).

The following subjects within the land transport and biomass currently do not
evolve with the myopic approach:

- The percentage of electric vehicles that allow demand-side management and
  vehicle-to-grid services.

- The annual biomass potential (default year and scenario for which potential is
  taken is 2030, as defined in config)

```yaml
{{ yaml_section("biomass.year") }}
```

### Configuration

For running myopic foresight transition scenarios, set in `config/config.yaml`:

```yaml
foresight: myopic
```

The following options included in the `config/config.yaml` file  are relevant for the
myopic code.

The `{planning_horizons}` wildcard indicates the year in which the network is
optimized. For a myopic optimization, this is equivalent to the investment year.
To set the investment years which are sequentially simulated for the myopic
investment planning, select for example:

```yaml
{{ yaml_section("scenario.planning_horizons", source="test/config.myopic.yaml") }}
```

**existing capacities**

Grouping years indicates the bins limits for grouping the existing capacities of
different technologies. Note that separate bins are defined for the power and
heating plants due to different data sources.

`grouping_years_power: [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020,
2025, 2030]`

`grouping_years_heat: [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2019]`

**threshold capacity**

If for a technology, node, and grouping bin, the capacity is lower than
threshold_capacity, it is ignored.

`threshold_capacity: 10`

**conventional carriers**

Conventional carriers indicate carriers used in the existing conventional
technologies.

    conventional_carriers:

    \- lignite

    \- coal

    \- oil

    \- uranium

### Options

The total carbon budget for the entire transition path can be indicated in the
[sector_opts](https://github.com/PyPSA/pypsa-eur-sec/blob/f13902510010b734c510c38c4cae99356f683058/config.default.yaml#L25)
in `config/config.yaml`. The carbon budget can be split among the
`planning_horizons` following an exponential or beta decay. E.g. `'cb40ex0'`
splits a carbon budget equal to 40 Gt $_{CO_2}$ following an exponential
decay whose initial linear growth rate r is zero. They can also follow some
user-specified path, if defined [here](https://github.com/PyPSA/pypsa-eur-sec/blob/413254e241fb37f55b41caba7264644805ad8e97/config.default.yaml#L56).
The paper [Speed of technological transformations required in Europe to achieve
different climate goals (2022)](https://doi.org/10.1016/j.joule.2022.04.016)
defines CO_2 budgets corresponding to global temperature increases (1.5C -- 2C)
as response to the emissions. Here, global carbon budgets are converted to
European budgets assuming equal-per capita distribution which translates into a
6.43% share for Europe. The carbon budgets are in this paper distributed
throughout the transition paths assuming an exponential decay. Emissions e(t) in
every year t are limited by

$$
e(t) = e_0 (1+ (r+m)t) e^{-mt}
$$

where r is the initial linear growth rate, which here is assumed to be r=0, and
the decay parameter m is determined by imposing the integral of the path to be
equal to the budget for Europe. Following this approach, the CO_2 budget is
defined. Following the same approach as in this paper, add the following to the
`scenario.sector_opts` E.g.  `-cb25.7ex0` (1.5C increase) Or `cb73.9ex0`
(2C increase). See details in Supplemental Note S1 [Speed of technological
transformations required in Europe to achieve different climate goals (2022)](https://doi.org/10.1016/j.joule.2022.04.016).

### National CO2 budgets

In addition to the system-wide carbon budget, individual countries can be given
their own CO2 budget to model differentiated, country-specific climate targets
(for example, Austria's path to climate neutrality by 2040 following its
*Klimaschutzgesetz* targets). National budgets are available in myopic foresight
and are configured under ``solving: constraints: co2_budget_national`` as a
fraction of each country's 1990 emissions, keyed by country code and planning
horizon:

.. code:: yaml

  solving:
    constraints:
      co2_budget_national:
        AT:
          2020: 0.67   # 67% of 1990 emissions
          2030: 0.34
          2040: 0.00   # net-zero
          2050: -0.05  # net-negative

Leave the mapping empty (the default) to disable the feature. The absolute
1990 reference is taken from ``co2_totals.csv`` for the emission sectors enabled
in the ``sector`` configuration, scaled to the number of years represented by
the planning horizon.

For every listed country and planning horizon a constraint is added that
balances all emissions at the ``co2 atmosphere`` bus:

1. all links of the country connecting to the ``co2 atmosphere`` bus are
   identified, across all of their bus ports;
2. emissions (positive and negative) are accounted for in the region where the
   technology causes them, so emissions in one region of a country can be
   compensated by negative emissions (CCS, DAC, CCU) in another region of the
   *same* country;
3. aviation emissions are scaled by the domestic-to-total aviation ratio
   (from ``energy_totals.csv``) so that international aviation is excluded;
4. the resulting sum is constrained to be below the national budget:
   :math:`\sum \text{emissions} \le \text{budget}_{ct}`.

Unlike the PyPSA-DE implementation this feature is derived from, no distinction
is made between synthetic and fossil fuels: all emissions are booked where
combustion occurs, and a country cannot use synthetic-fuel production to offset
emissions in another country. The constraint can be combined with the global
``co2limit``/carbon budget, which continues to cap system-wide emissions.

### General myopic code structure

The myopic code solves the network for the time steps included in
`planning_horizons` in a recursive loop, so that:

1. The existing capacities (those installed before the base year are added as
   fixed capacities with p_nom=value, p_nom_extendable=False). E.g. for
   baseyear=2020, capacities installed before 2020 are added. In addition, the
   network comprises additional generator, storage, and link capacities with
   p_nom_extendable=True. The non-solved network is saved in
   `resources/run_name/networks`.

The base year is the first element in `planning_horizons`. Step 1 is
implemented with the rule add_baseyear for the base year and with the rule
add_brownfield for the remaining planning_horizons.

2. The 2020 network is optimized. The solved network is saved in
   `results/run_name/networks`

3. For the next planning horizon, e.g. 2030, the capacities from a previous time
   step are added if they are still in operation (i.e., if they fulfil planning
   horizon <= commissioned year + lifetime). In addition, the network comprises
   additional generator, storage, and link capacities with
   p_nom_extendable=True. The non-solved network is saved in
   `results/run_name/networks`.

Steps 2 and 3 are solved recursively for all the planning_horizons included in
`config/config.yaml`.

### Rule overview

- rule add_existing baseyear

  The rule add_existing_baseyear loads the network in
  `resources/run_name/networks` and performs the following operations:

  1. Add the conventional, wind and solar power generators that were installed
     before the base year.

  2. Add the heating capacities that were installed before the base year.

  The existing conventional generators are retrieved from the [powerplants.csv
  file](https://pypsa-eur.readthedocs.io/en/latest/preparation/build_powerplants.html?highlight=powerplants)
  generated by pypsa-eur which, in turn, is based on the [powerplantmatching](https://github.com/PyPSA/powerplantmatching) database.

  Existing wind and solar capacities are retrieved from [IRENA annual statistics](https://www.irena.org/Statistics/Download-Data) and distributed among the
  nodes in a country proportional to capacity factor. (This will be updated to
  include capacity distributions closer to reality.)

  Existing heating capacities are retrieved from the report [Mapping and
  analyses of the current and future (2020 - 2030) heating/cooling fuel
  deployment (fossil/renewables)](https://ec.europa.eu/energy/studies/mapping-and-analyses-current-and-future-2020-2030-heatingcooling-fuel-deployment_en?redir=1).

  The heating capacities are assumed to have a lifetime indicated by the
  parameter lifetime in the configuration file, e.g 25 years. They are assumed
  to be decommissioned linearly starting on the base year, e.g., from 2020 to
  2045.

  Then, the resulting network is saved in
  `resources/run_name/networks`.

- rule add_brownfield

  The rule add_brownfield loads the network and reads the capacities optimized
  in the previous time step and add them to the
     network if they are still in operation (i.e., if they fulfill planning
     horizon < commissioned year + lifetime)
