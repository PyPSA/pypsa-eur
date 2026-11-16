<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Rules Overview

The pages in this section document the rules of the [workflow](../introduction.md#workflow),
grouped by topic and roughly in the order in which they run.

| Page | What it covers | Used by |
|------|----------------|---------|
| [Retrieving Data](retrieve.md) | Downloads of all external datasets, and the rules that need credentials | both |
| [Regions and Weather](regions.md) | Country, administrative and offshore shapes, weather cutouts | both |
| [Transmission Grid](grid.md) | Base network, transmission projects, line rating, simplification and clustering | both |
| [Electricity Demand and Supply](electricity.md) | Electricity demand, power plants, renewable potentials and profiles | both |
| [Population and Energy Balances](balances.md) | Population layouts, energy and emission totals per country and sector | both, mostly sector-coupled |
| [Heating](heating.md) | Heat demand, district heating, heat pumps, thermal storage, renovation | sector-coupled |
| [Transport](transport.md) | Land transport, aviation and shipping demand, EV profiles | sector-coupled |
| [Industry](industry.md) | Industrial production and energy demand per region | sector-coupled |
| [Costs and Prices](costs.md) | Technology costs, fuel and CO~2~ prices | both |
| [Biomass, Gas and Carbon](carriers.md) | Biomass potentials, gas network, salt caverns, CO~2~ storage | sector-coupled |
| [Composing and Solving](solving.md) | Assembling the network per planning horizon and solving it | both |
| [Plotting and Summaries](plotting.md) | Maps, time series plots and summary tables of the solved networks | both |

## Reading a rule entry

Each entry is generated from the rule definition in `rules/*.smk`
(e.g. [build_renewable_profiles][]):

- The summary below the rule name is the rule's docstring. Snakemake also
  prints it as the job message.
- Parts in braces are [wildcards](../wildcards.md), filled in from the
  requested file name.
- An input marked *depends on configuration* is chosen when the workflow runs.
- Changing one of the listed settings makes Snakemake rerun the rule.

`snakemake --list-rules` prints all rules with their summaries. Whole workflow
stages are bundled in [collection targets](../introduction.md#collection-targets),
and the full rule graph is shown on the [start page](../index.md#workflow).
