<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Transport {#design-transport}

{{ scope(electricity="off") }}

Transport comprises land transport, aviation and shipping. The fuel mix of
transport is exogenous, because it follows vehicle fleets and policy rather
than energy system cost, and it is given per planning year so that it can
trace a transition path [@brownSynergiesSector2018a].

## Demand

Annual energy demands per country and mode come from energy balances
[@jrcIDEESIntegrated2017; @eurostatEnergyBalances2021]. Road and rail transport are combined into one
land transport demand, except for electrified rail, which is already part of
the electricity demand. Aviation covers domestic and international flights,
and a demand factor allows scenarios with growing or shrinking air traffic.
Shipping is split into domestic and international traffic.

**Location.** Land transport and domestic shipping are distributed to
regions by population. International shipping demand is distributed by the
trade volume of the ports in each region [@worldbankWorldBank].

**Time series.** Hourly land transport demand follows weekly traffic count
profiles [@bundesanstaltfurstrassenwesenAutomatischeZahlstellen2021], shifted
to the local time of each country. A temperature dependence adds the energy
for heating and cooling the vehicle on cold and hot days. Aviation and
shipping demands are constant over the year.

## Supply

**Land transport** is split into three drivetrains whose shares are
prescribed: battery electric vehicles, fuel cell vehicles and internal
combustion engines. Each drivetrain converts the same transport service with
its own efficiency, so that electrification reduces final energy demand.
Electric vehicles draw from the low-voltage electricity bus, fuel cell
vehicles add a hydrogen demand at the regional hydrogen bus, and combustion
engines add an oil demand, see [Oil](oil.md).

**Aviation** consumes kerosene. It is an oil demand that can be met by fossil
kerosene or synthetic fuels, and optionally by kerosene made from methanol.

**Shipping** can run on oil, methanol or hydrogen, with prescribed shares per
planning year. The energy for liquefying hydrogen can be included where
hydrogen is used.

## Storage

Electric vehicles are more than a load. Their batteries are represented as an
aggregated storage per region with a charging connection whose availability
follows the driving profile: cars that drive cannot charge. Smart charging
lets a share of the fleet shift charging within the day, subject to a minimum
state of charge in the morning so that the batteries cannot act as seasonal
storage. Vehicle-to-grid lets the same share feed electricity back. The cost
of the vehicle batteries is not part of the system cost, since they are
bought for mobility.

## Further reading

- Rules: [build_transport_demand][], [build_shipping_demand][],
  [build_energy_totals][]
- Configuration: [sector](../configuration.md#sector_cf), in particular the
  `land_transport_*`, `shipping_*` and `aviation_*` settings
