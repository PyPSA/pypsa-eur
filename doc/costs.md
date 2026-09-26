<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Techno-economic assumptions {#costs}

{{ scope() }}

Every technology in the model is described by a small set of parameters:
investment cost, fixed and variable operation and maintenance cost, fuel
cost, efficiency, lifetime, discount rate and carbon dioxide intensity. They
come from one shared database so that all carriers and sectors are built on
consistent assumptions.

## Source

The database is maintained in the separate repository
[technology-data](https://github.com/pypsa/technology-data), which compiles
values from public sources, above all the technology catalogues of the
Danish Energy Agency [@DEA], and projects them to future years. PyPSA-Eur
retrieves one table per year. A pinned version of the database can be
selected so that results stay reproducible, see
[Managing data versions](data_sources.md#managing_data_versions).

## Cost years

By default each planning horizon uses the cost assumptions of its own year,
so that a pathway sees falling costs for maturing technologies. A single
reference year can be set instead, which is useful for sensitivity studies
where only the demand or the emission limit should change between runs.

## Annualisation

The optimisation compares investments with operating costs over one year.
Overnight investment costs are therefore turned into an annuity with the
discount rate $r$ over the economic lifetime $n$ using the annuity factor

$$
a = \frac{1-(1+r)^{-n}}{r},
$$

and the fixed operation and maintenance cost is added to obtain the annual
capital cost of one unit of capacity. The discount rate expresses the cost of
capital, see [Foresight](design/foresight.md#discount-rates) for its
interplay with the social discount rate of pathways. Marginal costs follow
from fuel cost, variable operation and maintenance cost and efficiency.

## Overriding assumptions

Selected values can be overridden in the configuration, for example to test
a cheaper electrolyser or a higher gas price, without editing the database.
Missing values are filled with defaults.

## Further reading

- Configuration: [costs](configuration.md#costs_cf)
