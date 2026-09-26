<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Industry {#design-industry}

{{ scope(electricity="off") }}

Industry differs from the other sectors in two ways. Its energy demand is tied
to the production of specific materials, and part of its emissions does not
come from burning fuels but from the chemistry of the process itself. The
model therefore builds industrial demand bottom-up from material production
and subsector-specific energy intensities
[@neumannPotentialRole2023; @materialeconomicsIndustrialTransformation2019].

## Demand

The method has three steps.

1. **Energy intensities.** For each industrial subsector, the energy demand
   and process emissions per unit of material output are derived from
   statistics of today's plants
   [@jrcIDEESIntegrated2017]. A second
   set of intensities describes the subsector after the transformations
   listed below, for example a fully electrified furnace. Between the two,
   the intensities are interpolated along the planning years.
2. **Production.** Material production per country is taken from statistics
   for a reference year and held constant, except where an exogenous route
   change is assumed, for example a switch from primary to secondary steel.
3. **Location.** Country totals are distributed to regions using a database
   of georeferenced industrial sites [@manzGeoreferencedIndustrial2018],
   supplemented by plant trackers for steel and cement and a list of ammonia
   plants, so that a cement region gets cement demand and a steel region gets
   steel demand. Where a subsector has no sites in a country, population is
   used instead.

The result is a demand per region for each energy carrier consumed by
industry: electricity, methane, hydrogen, solid biomass, coal and coke,
naphtha, methanol, ammonia, low-temperature heat and, as a separate stream,
process emissions. The choice between fossil and synthetic supply of these
carriers is then made by the optimisation.

### Subsectors

**Iron and steel.** Today's primary route reduces ore with coke in blast
furnaces and carries large process emissions; the secondary route melts scrap
in electric arc furnaces. A third route reduces ore with hydrogen and finishes
in an electric arc furnace, which avoids the process emissions
[@voglAssessmentHydrogen2018; @hybritSummaryFindings2021]. The shares of the
three routes are exogenous per planning year; the model then sees a demand
for hydrogen, electricity and methane.

**Chemicals.** Ammonia, methanol, chlorine and high-value chemicals are
separated from the aggregate statistics because their production routes
differ [@bazzanellaLowCarbon2017; @leviMappingGlobal2018]. Ammonia demand is
met from an explicit ammonia carrier, see [Ammonia](ammonia.md). Methanol and
chlorine follow fixed production statistics. High-value chemicals need a
hydrocarbon feedstock, naphtha, whose origin, fossil, synthetic or biogenic,
is optimised. Recycling and reuse of plastics reduce the primary feedstock
and the associated process emissions [@circular_economy; @meysAchievingNetzero2021; @kullmannValueRecycling2022]; their shares are
exogenous. The carbon in plastics that are neither recycled nor reused is
either released, permanently landfilled, or burned in waste-to-energy plants
with or without capture, in prescribed shares.

**Non-metallic minerals.** Cement carries the largest process emissions in
European industry from the calcination of limestone
[@fennellDecarbonizingCement2021]. These emissions can be captured, which the
optimisation decides on cost. High-temperature heat for kilns comes from
methane. Ceramics and glass are assumed to be electrified
[@furszyferdelrioDecarbonizingGlass2022].

**Non-ferrous metals.** Aluminium dominates the subsector. The primary route
from bauxite is electricity intensive and carries process emissions; the
secondary route remelts scrap. Their shares are exogenous. Other metals are
electrified but keep their process emissions.

**Other subsectors.** Pulp and paper, food, textiles, machinery, transport
equipment, wood products and the remaining industries have no process
emissions [@sovacoolDecarbonizingFood2021]. Their low-temperature heat comes
from biomass and the rest is electrified.

## Supply

Where a process has both a fossil and an electric variant, the electric one is
assumed [@lechtenbohmerDecarbonisingEnergy2016]. Low- and medium-temperature
process heat is supplied by solid biomass; high-temperature heat by methane
or electricity where electric processes exist
[@rehfeldtBottomupEstimation2018; @naeglerQuantificationEuropean2015].
Hydrogen for high-temperature heat is not represented. Which of the demanded
carriers is produced how, from fossil, biogenic or synthetic sources, is
decided on the carrier pages: [Hydrogen](hydrogen.md),
[Methane](methane.md), [Oil](oil.md), [Methanol](methanol.md) and
[Biomass](biomass.md).

Process emissions are tracked separately from energy-related emissions. They
can be captured at a given capture rate, which turns them into a captured
CO~2~ stream for usage or sequestration, see
[Carbon management](carbon.md). Emissions from refining fossil oil can be
added in proportion to fossil oil consumption.

## Agriculture, forestry and fishing

Energy demand of agriculture, forestry and fishing per country is split into
electricity, low-temperature heat and machinery fuel. Electricity and heat are
constant loads distributed by population, with heat assigned to rural
services heating. Machinery is supplied by oil and, with an exogenous share
per planning year, by electricity.

## Further reading

- Rules: [build_industry_sector_ratios][],
  [build_industrial_production_per_country][],
  [build_industrial_production_per_country_tomorrow][],
  [build_industrial_distribution_key][], [build_industrial_energy_demand_per_node][]
- Configuration: [industry](../configuration.md#industry_cf)
