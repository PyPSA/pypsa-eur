<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Foresight {#foresight}

{{ scope() }}

Foresight describes how the model treats time beyond the operational year:
whether it plans a single target year from scratch, walks through a sequence
of investment years, or optimises all years at once. All three modes share the
same workflow and the same list of planning horizons; they differ in what each
horizon knows about the others.

## Planning horizons

Every run has a list of planning horizons, the years for which a network is
composed and solved. A horizon fixes the year of the cost and technology
assumptions and the year of all exogenous transition paths, such as the share
of electric vehicles, the share of district heating or the route mix in steel
making. These paths are given as values per year and interpolated for the
horizons in between, so that a pathway can be run at any spacing of years.
With a single horizon the model is a snapshot of one year; with several it is
a pathway.

!!! warning "One horizon for overnight"
    Overnight planning accepts a single planning horizon. Several horizons
    require myopic or perfect foresight.

## Overnight

Overnight, or greenfield, planning builds the system for one horizon from
scratch. Assets that cannot be rebuilt are carried in as given: the existing
transmission grid, hydropower and nuclear plants, and today's wind and solar
capacities as a lower bound. Everything else is optimised as if it could be
constructed overnight. In the electricity-only model, all existing thermal
plants are carried in at their capacity. This is the default and
the simplest mode: no assumptions about the existing stock and its ageing are
needed, and the result shows the cost-optimal system for the target year
regardless of the path to it.

## Myopic

Myopic planning walks through the horizons in sequence
[@victoriaEarlyDecarbonisation2020]. The first horizon
starts from the existing stock of power plants, heating systems and renewable
capacities. Power plants carry their commissioning year and a technology
lifetime and retire when the sum is reached; heating systems, for which no
plant-level ages exist, are retired in equal shares over an assumed lifetime.
In each horizon the model adds capacity to meet that year's demand and
emission limit, knowing only the present. Assets built in one horizon stay in
the network until their lifetime expires, so the next horizon inherits them.

This mode captures path dependency and stranded assets, but it can be
short-sighted: it cannot anticipate a tighter emission limit in a later
horizon. The exogenous transition paths evolve along the horizons; a few
inputs, such as the biomass potential scenario, are held fixed.

!!! note "Phase-outs follow lifetimes"
    Configured phase-out years for nuclear or coal apply in perfect foresight
    only. In myopic planning, plants retire when their lifetime ends.

## Perfect foresight

Perfect foresight optimises all horizons together in one multi-period problem
[@SpeedTechnological2022].
Investments in one year are made with full knowledge of the demands, costs and
emission limits of all later years, so the result is the cost-optimal pathway
under the assumption of an omniscient planner. Emission limits can be set per
year or as one budget over the whole pathway, which lets the model decide when
to decarbonise. Announced phase-outs of nuclear or coal in specific countries
are enforced by the year they take effect.

## Discount rates

Two discount rates play different roles. The financial discount rate turns
overnight investment costs into annuities and reflects the cost of capital; a
high rate penalises capital-intensive technologies such as wind and solar
relative to fuel-intensive ones. The social discount rate weighs the costs of
different years against each other in a pathway and expresses how much less a
cost in the future counts than the same cost today. It matters in perfect
foresight, where all years are optimised together, and when pathway costs are
compared after the fact.

## Emission limits across modes

Emission limits are set in the same way for all modes, either as absolute
amounts or relative to a historical baseline year, and either per horizon or
as one number. The difference is in how they bind: in overnight and myopic
planning a single number is the cap of every solved horizon; in perfect
foresight it is a budget on the sum over all horizons. A cumulative budget
derived from a temperature target can therefore be handed to perfect
foresight as one sum, or turned into per-horizon limits along a chosen
trajectory, for example a steady decline, for the other modes. Which
greenhouse gases are counted in the baseline is a choice, with energy-related
CO~2~ as the default.

## Further reading

- Rules: [compose_network][], [solve_network][], [add_existing_baseyear][],
  [add_brownfield][], `prepare_perfect_foresight`
- Configuration: [foresight](../configuration.md#foresight_cf),
  [planning_horizons](../configuration.md#planning_horizons_cf),
  [existing_capacities](../configuration.md#existing_capacities_cf),
  [co2_budget](../configuration.md#co2_budget_cf)
