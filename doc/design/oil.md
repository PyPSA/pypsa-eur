<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Oil {#design-oil}

{{ scope(electricity="partial") }}

Oil products stand for naphtha, kerosene, diesel and petrol alike. They are
one carrier whose fossil, biogenic and synthetic origins compete at the same
bus, so that the model can replace fossil oil without deciding in advance by
what.

=== "Sector-coupled"

    Oil is a carrier with the demands, supply routes and regional demand
    resolution described on this page.

=== "Electricity-only"

    Oil appears only as existing power plants with a fuel price.

## Demand

Naphtha is the feedstock of the chemicals industry [@leviMappingGlobal2018].
Kerosene fuels aviation. Oil products also fuel combustion engine vehicles,
ships and agricultural machinery in the shares prescribed for those sectors,
see [Transport](transport.md) and [Industry](industry.md). Oil boilers in
buildings and existing oil power plants complete the demand.

## Supply

Oil products come from fossil crude, optionally with the emissions of
refining, or from synthetic routes: Fischer-Tropsch synthesis from hydrogen
and captured CO~2~, liquefaction of solid biomass, and electrobiofuels that
combine biomass with hydrogen to raise the fuel yield. Biomass liquefaction
can be built with carbon capture. Today's biofuels from unsustainable
feedstocks also feed the oil supply and are phased out, see
[Biomass](biomass.md). Imports at a fixed price can be enabled. The share of
fossil and synthetic supply is a result of the optimisation.

## Transport

Liquids are cheap to transport, so supply is copperplated at a single
European node. Oil demand is placed in the regions by default, which keeps
track of where the fuel is burned without constraining its transport.

## Further reading

- Configuration: [sector](../configuration.md#sector_cf), in particular
  `regional_oil_demand`, `oil_refining_emissions`, `biomass_to_liquid` and
  `electrobiofuels`
