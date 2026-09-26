<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Carbon management {#design-carbon}

{{ scope(electricity="partial") }}

In deep decarbonisation scenarios, carbon dioxide becomes a carrier that is
captured, transported, used and stored [@hofmannH2CO2Network2025]. The model tracks it
explicitly, which is what makes net-zero and net-negative scenarios
consistent. Capture, usage and sequestration are represented as separate
technologies rather than under one umbrella, so that each can be optimised on
its own and the interplay between them becomes visible: captured carbon can
feed fuel synthesis or go underground, and the model chooses.

=== "Sector-coupled"

    CO~2~ is a carrier with an atmosphere bus, capture technologies, usage,
    pipelines and sequestration as described on this page.

=== "Electricity-only"

    There is no CO~2~ carrier. Emissions follow from the carbon intensity of
    each generator's fuel, and the emission limit applies to their sum.

## Accounting

Every combustion and process step that releases CO~2~ is linked to an
atmosphere bus, and every capture step withdraws from it. The emission limit
is a constraint on the net flow into that bus. It can be set as an absolute
amount or relative to a historical baseline year, per planning year or, in
perfect foresight, as a budget over the whole pathway, see
[Foresight](foresight.md). Fossil, biogenic and synthetic fuels are
distinguished so that burning synthetic fuel made from captured CO~2~ is
neutral and sequestering biogenic CO~2~ is negative.

## Capture

CO~2~ can be captured at point sources and from the air:

- industrial process emissions, for example from cement calcination
  [@kuramochiComparativeAssessment2012];
- fuel combustion in industry, in combined heat and power plants and in gas
  power plants with oxy-fuel combustion;
- hydrogen production by steam methane reforming;
- biomass conversion routes and biogas upgrading;
- **direct air capture** [@breyerCarbonDioxide2020], which consumes
  electricity and low-temperature heat from district heating.

For each point source the model can choose between the variant with and
without capture; capture is never forced. Each technology has a capture rate,
and the remainder still reaches the atmosphere. Capture consumes electricity,
and direct air capture also heat. Captured CO~2~ arrives at a regional
CO~2~ bus by default, or at one European bus if CO~2~ is not resolved
spatially.

Because most capture options are attached to plants whose main product is
something else, heat, electricity, hydrogen or fuel, the cost of capturing one
more tonne differs widely between them and depends on how the plant is
operated. The optimisation therefore establishes an implicit merit order of
capture: cheap sources such as process emissions are used first, and expensive
ones such as direct air capture only when the emission limit or the demand
for carbon requires it. With a CO~2~ network, cheap capture in one region can
serve usage or sequestration in another.

## Usage

Captured CO~2~ is the carbon in synthetic methane, methanol and Fischer-Tropsch
fuels, see [Methane](methane.md), [Oil](oil.md) and [Methanol](methanol.md). These
drop-in fuels are not limited in quantity; their volume follows from the
demands they serve and from the cost of the hydrogen and carbon they need.
Using captured carbon closes the loop only if the carbon came from the air or
from biomass; the accounting above makes sure that fossil carbon released
later still counts. The clearest case is naphtha: the carbon locked into
plastics is released at their end of life unless they are landfilled or
burned with capture, see [Industry](industry.md).

## Storage

Two kinds of storage are distinguished. **Tanks** hold CO~2~ for a short time
to decouple capture from usage or transport. **Sequestration** removes CO~2~
permanently in geological formations such as saline aquifers and depleted gas
and oil fields [@martin-robertsCarbonCapture2021]. The annual amount that
may be sequestered is limited, and the limit can rise over the planning
years. Storage sites can be resolved regionally from a geological database
of potentials, with a choice between offshore formations only, which are
larger and far from settlements, and onshore formations as well.

!!! warning "The sequestration limit shapes the result"
    The limit is a deliberate choice: it keeps sequestration for
    hard-to-abate emissions and prevents it from becoming a backstop that
    offsets fossil fuel use which could be avoided. Results depend strongly
    on it.

## Transport

CO~2~ can be transported in a pipeline network between regions, onshore and
submarine, whose capacity is optimised, with candidate routes following the
existing electricity network. An optional compression stage with its
electricity demand can be added for transport in dense phase. Alternatively
CO~2~ is copperplated, with a transport and storage cost per tonne. Venting
captured CO~2~ back to the atmosphere can be allowed as an escape valve.

## Further reading

- Rules: [build_co2_sequestration_potentials][],
  [build_clustered_co2_sequestration_potentials][], [build_co2_totals][]
- Configuration: [co2_budget](../configuration.md#co2_budget_cf) and the
  `co2_*`, `dac` and `cc_fraction` settings under
  [sector](../configuration.md#sector_cf)
