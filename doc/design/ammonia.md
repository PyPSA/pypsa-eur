<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Ammonia {#design-ammonia}

{{ scope(electricity="off") }}

Ammonia is both a chemical, mostly for fertiliser, and a potential energy
carrier, since it is easier to ship and store than hydrogen
[@wangGreeningAmmonia2018]. The model can represent it as a carrier of its
own or fold it into the hydrogen and electricity it is made from.

## Demand

Ammonia demand comes from the chemicals industry, where production per
country follows statistics
[@unitedstatesgeologicalsurveyAmmoniaProduction2021] and is located at the
existing plants, see [Industry](industry.md). With the ammonia carrier
disabled, this demand is translated into its hydrogen and electricity
equivalents instead.

## Supply

Ammonia is synthesised from hydrogen and electricity in Haber-Bosch plants,
whose capacity is optimised. The waste heat of synthesis can supply district
heating. Ammonia is one of the carriers that can be imported from outside
Europe at a fixed price, see [Model overview](overview.md#design-overview).

## Conversion

Ammonia can be cracked back to hydrogen, which is one way to turn imported
ammonia into hydrogen for other uses, see [Hydrogen](hydrogen.md).

## Storage

Ammonia is stored in tanks, which are far cheaper per unit of energy than
hydrogen tanks, so that stored ammonia can buffer hydrogen in bound form.

## Transport

Ammonia is copperplated by default: one European bus serves all regions.
Optionally it can be resolved regionally, in which case each region has its
own bus.

## Further reading

- Rules: [build_industrial_distribution_key][],
  [build_industrial_energy_demand_per_node][]
- Configuration: the `ammonia` and `imports` settings under
  [sector](../configuration.md#sector_cf)
