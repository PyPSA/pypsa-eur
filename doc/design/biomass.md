<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Biomass {#design-biomass}

{{ scope(electricity="partial") }}

Biomass is the one renewable resource that can supply high-temperature heat,
dispatchable power and carbon at once, and it is scarce. How much of it the
model may use, and for what, is therefore a defining assumption
[@millingerDiversityBiomass2025; @MILLINGER2022120016].

=== "Sector-coupled"

    Solid biomass, biogas and municipal waste are carriers with regional
    potentials and the conversion routes described on this page.

=== "Electricity-only"

    Biomass appears only as existing power plants with a fuel price.

## Demand

Solid biomass supplies low- and medium-temperature process heat in industry
and fuels combined heat and power plants and boilers in heating. Biogas is
demanded indirectly through the methane it is upgraded to. These demands are
described on the [Industry](industry.md), [Heating](heat.md) and
[Methane](methane.md) pages.

## Supply

Regional potentials are taken from a European assessment of bioenergy
potentials [@ruizENSPRESOOpen2019] that lists many feedstocks per
statistical region for several years and availability scenarios. Each
feedstock is assigned to one of four classes:

| Class | Typical feedstocks | Use in the model |
|---|---|---|
| Solid biomass | Agricultural and forestry residues, landscape care wood | Combustion and conversion |
| Biogas | Manure, sewage sludge | Upgraded to methane |
| Municipal solid waste | Biodegradable waste | Waste-to-energy, optional |
| Not included | Energy crops, roundwood | Excluded |

The default excludes crops that compete with food and primary wood whose
sustainability is contested [@bentsenCarbonDebt2017]. Two transitions
overlap along the planning years. Feedstocks that are used today but excluded
in the future, such as energy crops and roundwood, are represented explicitly
as unsustainable solid biomass, biogas and bioliquids with a minimum use that
is phased out, so that early horizons reflect today's consumption. At the
same time the sustainable potential is phased in, so that it is not fully
available in the first horizons. Imports of solid biomass at a fixed price
can be allowed with a cap. No further biomass is assumed to arrive from
outside Europe.

Potentials are mapped from statistical regions to model regions by area
overlap.

## Conversion

Solid biomass can be burned in combined heat and power plants with or without
carbon capture, converted to liquid fuels, to methanol, to hydrogen or to
synthetic natural gas, and combined with hydrogen in electrobiofuels.
Unsustainable bioliquids feed the oil supply directly. Biogas is upgraded to
methane before it enters the gas system. Which conversion routes are built is
a result of the optimisation; the carbon captured from biogenic sources
counts as negative emissions when it is sequestered, see
[Carbon management](carbon.md).

## Transport

Solid biomass can be one European pool or regionally resolved. When resolved,
transport between regions can be free or priced per distance with
country-specific costs, which reflects road transport of bulky feedstocks.
Biogas follows the gas network if that is resolved.

## Further reading

- Rules: [build_biomass_potentials][], [build_biomass_transport_costs][]
- Configuration: [biomass](../configuration.md#biomass_cf) and the
  `biomass_*` settings under [sector](../configuration.md#sector_cf)
