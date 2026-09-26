<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Electricity {#design-electricity}

{{ scope() }}

Electricity is the backbone of the model [@PyPSAEur]. Every region has an
electricity bus at transmission level, and in the sector-coupled model most
other carriers are produced from or converted into electricity.

## Demand

Hourly electricity demand per country is taken from historical load statistics
of the transmission system operators and the Open Power System Data platform
[@muehlenpfordtOpenPower2020], filled and cleaned where the records have gaps,
and optionally supplemented with synthetic time series for years or countries
without records.

Within a country, demand is distributed to regions with a spatial key. For
households and services the key combines population and economic activity;
where available, gridded estimates of annual electricity consumption refine
the picture. Industrial electricity demand is placed at the sites of
energy-intensive industry instead, see [Industry](industry.md).

=== "Sector-coupled"

    Two parts are subtracted from the historical load and modelled
    separately: the electricity used today for heating, so that power-to-heat
    can be optimised on its own, and the electricity used in industry, so
    that further electrification of industry can be represented. Cooling
    demand stays in the electricity demand and is not changed.

=== "Electricity-only"

    The historical load is used as it is, distributed with the same spatial
    key. Nothing is subtracted for heating or industry.

## Supply

### Renewables

Wind and solar capacities are decided by the optimisation within limits set by
land availability and the weather.

**Potentials.** For each technology and region, the eligible land is
determined from land cover [@europeanenvironmentagencyeeaCorineLand], nature
protection areas [@europeanenvironmentagencyeeaNatura2000], bathymetry
[@gebcoGEBCO2014] and shipping lanes, together with distance criteria such
as a minimum distance from shore for offshore wind. The eligible area is
multiplied with an allowed deployment density to obtain the installable
capacity [@mckennaHighresolutionLargescale2022; @boschTemporallyExplicit2018].
The cost of connecting new capacity to the grid follows from the average
distance of the eligible sites to the substation, and for offshore wind from
the distance to the landfall. Offshore wind is split into near-shore sites
with AC connection, far-shore sites with DC connection and deep-water sites
for floating turbines. Rooftop solar potential is derived from population
density [@bodisHighresolutionGeospatial2019]. Each region can further be
split into resource classes so that sites of different quality within one
region are represented separately.

**Time series.** Hourly capacity factors are computed from gridded weather
data with [atlite](https://atlite.readthedocs.io)
[@hofmannAtliteLightweight2021]: wind speeds through turbine power curves,
solar irradiation through panel models with fixed or single-axis tracking
mounts, and hydro inflow from surface runoff scaled to historical generation
statistics. The grid-cell time series are aggregated to regions in proportion
to the mean capacity factor of each cell, which assumes that new capacity is
built where the resource is best.

**Hydropower.** Run-of-river, reservoir and pumped-hydro plants are taken from
the power plant database with their existing capacities and are not
expandable. Reservoir inflow follows the runoff-based time series.

**Geothermal.** Enhanced geothermal systems can generate electricity, and heat
for district heating, where the geological potential allows. Their potential
and cost per region come from a gridded assessment.

### Conventional plants

Existing thermal power plants come from the open
[powerplantmatching](https://github.com/PyPSA/powerplantmatching) database
[@gotzensPerformingEnergy2019], which merges several public sources into one
list with location, fuel, technology, capacity and commissioning year.
Existing wind and solar capacities are taken from the plant-level records in
the same database, or optionally from national statistics, and placed within
a country in proportion to resource quality and potential.

New gas turbines in open and combined cycle can be built. By default only gas
turbines may be built new among the thermal plants; nuclear, coal, lignite and
oil plants enter transition pathways as existing capacities.

=== "Sector-coupled"

    Thermal plants are links from a fuel bus, so that their fuel can be
    fossil, biogenic or synthetic, see [Methane](methane.md). Further
    options to generate electricity from other carriers are hydrogen
    turbines and fuel cells, combined heat and power plants, methanol-fired
    plants and oxy-fuel gas plants with inherent carbon capture.

=== "Electricity-only"

    There are no fuel buses. Conventional plants are generators with a fuel
    price and a carbon intensity per carrier, and the set of carriers that may
    be built new is chosen in the configuration. Optional operational detail
    for these plants, such as an operating reserve margin and linearised unit
    commitment, is only available here.

## Conversion

Electricity is the input to most other carriers: heat pumps and resistive
heaters make heat, electrolysers make hydrogen, and electric vehicles,
electrified industrial processes and direct air capture draw power. These
links are described on the pages of the carriers they produce. The flexibility
they add, by shifting their consumption to hours of abundant wind and solar,
is one of the main effects of sector coupling on the electricity system
[@brownSynergiesSector2018a].

## Storage

Electricity can be stored in stationary batteries and in the existing
pumped-hydro plants. Batteries are represented with separate energy and power
components, so that duration and capacity are sized independently, with
distinct costs for utility-scale and home batteries. Several battery
chemistries and other electricity storage technologies with different typical
durations, from hours to days, can be offered to the optimisation.

=== "Sector-coupled"

    Electric vehicle batteries add a large flexible storage, see
    [Transport](transport.md), and hydrogen with its underground storage acts
    as seasonal storage when it is re-electrified, see
    [Hydrogen](hydrogen.md).

=== "Electricity-only"

    Batteries and hydrogen storage are storage units with a fixed
    energy-to-power ratio per technology, or optionally stores with separate
    charging and discharging links. There is no hydrogen carrier beyond this
    storage.

## Transport

### Transmission grid {#transmission}

The transmission grid is derived from OpenStreetMap [@xiongModellingHighvoltage2025].
Substations, lines, cables and HVDC links at and above the highest voltage
levels are extracted, their tags cleaned and standardised, and substations
close to each other merged into one bus. Missing circuit counts are filled by
heuristics, and each line is assigned a standard line type of its voltage
level [@oedingElektrischeKraftwerke2011], from which impedance and thermal
rating follow. HVDC links keep their reported ratings and are connected to
the nearest AC bus through converter links. The dataset has been validated
against the transmission system operators' inventory statistics. An
alternative base network can be built from the ENTSO-E reference grid
[@tyndp2018].

Power flow on the AC network follows linearised power flow equations, so that
Kirchhoff's voltage law is respected [@horschLinearOptimal2018]; HVDC links
are transport links. A security margin reduces the usable line capacity to
approximate N-1 safe operation, and dynamic line rating from weather data can
be enabled.

Grid expansion is endogenous. Existing corridors can be reinforced, and the
total expansion can be capped relative to today's grid volume or cost.
Resistive losses on AC lines are approximated piecewise linearly, and losses
on HVDC links per distance.

### Distribution grid {#distribution}

=== "Sector-coupled"

    The grid below transmission level is not represented topologically.
    Instead, each region has a low-voltage bus connected to the transmission
    bus by a link whose capacity is optimised. Electricity demand, rooftop
    solar, heat pumps, resistive heaters, home batteries and electric vehicle
    chargers connect to the low-voltage side; large generators and storage
    connect to the transmission side. Distribution capacity is then expanded
    only where local generation and demand do not balance, and distribution
    losses are represented as a fixed share.

=== "Electricity-only"

    There is no distribution level; all technologies connect to the
    transmission bus.

## Further reading

- Rules: [build_electricity_demand][], [build_renewable_profiles][],
  [build_hydro_profile][], [build_powerplants][], [base_network][],
  [cluster_network][], [build_line_rating][]
- Configuration: [electricity](../configuration.md#electricity_cf),
  [renewable](../configuration.md#renewable_cf), [lines](../configuration.md#lines_cf),
  [load](../configuration.md#load_cf)
