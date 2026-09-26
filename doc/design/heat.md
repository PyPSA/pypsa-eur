<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Heating {#design-heat}

{{ scope(electricity="off") }}

Building heat is the largest energy demand in Europe after transport and has a
strong seasonal profile [@zeyenMitigatingHeat2021]. The model resolves space
and water heating of the residential and services sectors per region and
distinguishes buildings that are served by district heating from those with
individual heating systems.

## Demand

**Annual totals.** Annual heat demand per country is taken from energy
balances [@jrcIDEESIntegrated2017; @eurostatEnergyBalances2021] as useful energy, that is the heat delivered
rather than the fuel burned, and split into space heating and hot water, and
into residential and services buildings.

**Time series.** Space heat demand is spread over the year in proportion to
daily heating degree days computed from ambient temperature, so that cold days
carry more demand. Within the day, standard load profiles [@bdewBDEWHeat2021]
distribute the demand over the hours, with different shapes for weekdays and
weekends and for residential and services buildings. Hot water demand is
constant over the year.

**Heat systems.** Within a country, heat demand is split between urban and
rural areas using population density, and each part is distributed to regions
in proportion to its population. Dense urban areas can be served by district
heating. This results in five heat systems per region, each with its own
supply options:

| Heat system | Buildings served |
|---|---|
| Urban central | District heating in dense urban areas, residential and services |
| Residential urban decentral | Residential buildings in urban areas without district heating |
| Services urban decentral | Services buildings in urban areas without district heating |
| Residential rural | Residential buildings in areas with low population density |
| Services rural | Services buildings in areas with low population density, including agriculture |

!!! note "Three heat buses by default"
    The residential and services buses of the same type are merged by default
    to keep the problem small, leaving one district heating, one urban
    decentral and one rural bus per region.

The share of urban demand on district heating is
exogenous and can grow along the planning years up to a potential. District
heating networks incur lump-sum distribution losses. Spatial data on district
heating areas is used to estimate where heat sources and large storage can be
connected.

**Renovation.** Renovation of the thermal envelope reduces space heat demand.
It can be prescribed as an exogenous reduction, or it can be optimised per
region and heat system [@zeyenMitigatingHeat2021]. In the endogenous case,
costs per unit of energy saved are estimated from the building stock and
renovation costs, and several insulation depths are offered to the
optimisation as a step-wise supply curve, so that it can trade insulation
against supply capacity.

## Supply

### Decentral heating

Individual buildings can be heated with air-sourced heat pumps, ground-sourced
heat pumps in rural areas, resistive heaters, gas, oil and biomass boilers,
micro combined heat and power units, and solar thermal collectors. The
coefficient of performance of heat pumps follows the temperature difference
between the source and the sink [@staffellReviewDomestic2012]; it is lower on
cold days and more stable for ground-sourced units, whose source temperature
varies less. Solar thermal yield follows irradiation and collector
orientation.

### District heating

District heating networks have more options because they can host large
plants. Combined heat and power plants, fired with solid biomass or methane by
default and optionally with other fuels or with non-recyclable waste, run with
or without carbon capture in back-pressure mode with a fixed ratio of
electricity to heat. Large heat pumps, resistive heaters, boilers, fuel cells
and enhanced geothermal systems complete the picture.

Large heat pumps can draw on several sources. Besides ambient air, which is
the default, the model can exploit geothermal heat, rivers and the sea, each
limited by a regional potential derived from geological and hydrological data
within reach of the district heating areas. Their coefficient of performance
is approximated from the source temperature and the network supply
temperature with a thermodynamic model. The supply temperature itself is
approximated from ambient temperature following measured reference curves
[@pieperAssessmentCombination2019]: it is highest on cold days and lowest on mild days, and it
can decline over the planning years as networks modernise. Sources that are
hot enough can feed the network directly without a heat pump.

Waste heat from conversion processes is a further source: electrolysis, fuel
cells, Fischer-Tropsch synthesis, methanolisation and ammonia synthesis reject
heat that can be fed into district heating where these plants are located.

In transition pathways the existing heating stock matters. Installed
capacities of boilers, heat pumps and resistive heaters are taken from a
European survey, assigned a lifetime and decommissioned gradually, so that the
model replaces them over time, see [Foresight](foresight.md).

## Storage

**Water tanks** provide short-term storage for individual buildings, and
large **pit storages** provide seasonal storage for district heating. Both are
modelled as stores with standing losses; pit storage capacity can follow the
network temperatures [@sorknaesSimulationMethod2018], and supplemental heating covers the
times when the stored water is not hot enough. **Aquifer storage** is
available where suitable aquifers underlie district heating areas
[@jacksonAquiferThermal2024].

**Building thermal mass** can optionally be used for demand-side management
of residential space heating: supply may run ahead of demand within a day, as
long as the buffer is empty again at fixed hours, so that buildings cannot act
as seasonal storage. Surplus heat, for example from waste heat in summer,
can be vented so that heat supply never has to be curtailed.

## Further reading

- Rules: [build_daily_heat_demand][], [build_hourly_heat_demand][],
  [build_district_heat_share][], [build_cop_profiles][build_cop_profiles.run],
  [build_central_heating_temperature_profiles][build_central_heating_temperature_profiles.run],
  [build_ptes_operations][build_ptes_operations.run],
  [build_ates_potentials][], [build_retro_cost][], [build_solar_thermal_profiles][]
- Configuration: [sector](../configuration.md#sector_cf), in particular
  `district_heating`, `heat_pump_sources`, `tes` and `retrofitting`
