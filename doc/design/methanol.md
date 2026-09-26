<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Methanol {#design-methanol}

{{ scope(electricity="off") }}

Methanol is a liquid hydrogen derivative that is easy to store and ship. It
serves as a shipping fuel, a chemical feedstock and, optionally, as an
intermediate towards hydrogen, electricity and kerosene.

## Demand

Methanol is demanded by shipping in the share prescribed per planning year,
see [Transport](transport.md), and by the chemicals industry, where methanol
production follows statistics, see [Industry](industry.md).

## Supply

Methanol is synthesised from hydrogen and captured CO~2~, or produced from
solid biomass with or without carbon capture. The waste heat of synthesis can
supply district heating. A biomass route boosted with hydrogen, as
electrobiofuels are for oil, is not represented for methanol. Methanol can also be imported from outside Europe at
a fixed price, see [Model overview](overview.md#design-overview).

## Conversion

Optionally, methanol can be reformed to hydrogen, burned in power plants for
electricity, or converted to kerosene for aviation. These routes matter when
methanol is imported: they turn an easily shipped liquid into the carriers
that are needed inland.

## Transport

Methanol is copperplated at a single European node. Its demand can
optionally be placed in the regions, like oil demand.

## Further reading

- Configuration: the `methanol` block, `regional_methanol_demand` and
  `imports` under [sector](../configuration.md#sector_cf)
