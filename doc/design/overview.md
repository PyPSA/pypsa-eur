<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Model overview {#design-overview}

PyPSA-Eur is an open model of the European energy system at the resolution of
the electricity transmission network [@PyPSAEur]. In its electricity-only form
it covers generation, storage and transmission of electricity. In its
sector-coupled form it adds the heating, transport and industry sectors
together with the carriers that connect them: hydrogen, ammonia, methane,
oil, methanol, biomass and carbon dioxide
[@brownSynergiesSector2018a; @neumannPotentialRole2023].

The model builds a linear optimisation problem that plans infrastructure from
open data. Given demands, weather, potentials and techno-economic assumptions,
it decides which capacities to build and how to operate them so that the total
system cost is minimised while all demands are met and emission limits are
respected. Every step from raw data to solved model is a rule in a
[Snakemake](https://snakemake.github.io/) workflow [@snakemake], so the whole
chain is reproducible and can be re-run when data or assumptions change.

## Energy and carbon flows

The figure below shows how energy and carbon circulate in the sector-coupled
model. Primary energy enters on the left as weather-dependent renewables,
fossil fuels, biomass and imports. Conversion technologies link the carrier
buses in the middle. Final demands leave on the right by sector. Carbon
dioxide is tracked as a carrier of its own, from emission and capture to usage,
transport and sequestration.

![Energy and carbon flows in the sector-coupled model](../img/multisector_figure.png)

## Model variants

The workflow can stop after the electricity system or continue to the full
sector-coupled model, and it can plan a single year or a pathway of several
planning years, see [Foresight](foresight.md). Every design page starts with a
row of badges that says in which variants its content applies:

{{ scope(electricity="on", sector="partial", overnight="off") }}

Blue badges are the model variants, teal badges the foresight modes. A filled
badge means the page applies to that variant, an outlined badge means it
applies in a reduced form that the page explains in a pair of tabs, and a
struck badge means it does not apply. The matrix below gives the same
information for the whole model.

| Feature | Electricity-only | Sector-coupled | Depends on foresight |
|---|:-:|:-:|---|
| Electricity demand from historical load | :material-check: | :material-check: as residual after heat and industry are removed | |
| Wind, solar and hydro potentials and time series | :material-check: | :material-check: | |
| Conventional plants as explicit fuel carriers | :material-close: fuel price per generator | :material-check: | |
| Batteries and pumped hydro | :material-check: | :material-check: | |
| Transmission grid expansion | :material-check: | :material-check: | |
| Distribution grid level | :material-close: | :material-check: | |
| Heating, transport, industry, agriculture | :material-close: | :material-check: sectors can be switched off individually | |
| Hydrogen, ammonia, methane, oil, methanol, biomass | :material-close: hydrogen only as storage | :material-check: | |
| Carbon capture, usage, transport and sequestration | :material-close: emission limit only | :material-check: | |
| Existing stock with ages and retirement | :material-check: | :material-check: | myopic and perfect foresight |
| Phase-out of nuclear and coal by year | :material-check: | :material-check: | perfect foresight |
| Emission budget over the pathway | :material-check: | :material-check: | perfect foresight |

## Modelling principles

**Buses, links and stores.** The model is built from the components of
[PyPSA](https://docs.pypsa.org). Each carrier at each location is a bus that
must balance in every hour. Generators inject energy from outside the model,
loads withdraw it, stores and storage units shift it in time, and links convert
between carriers with an efficiency, for example a heat pump from electricity
to heat or an electrolyser from electricity to hydrogen. Transmission lines and
pipelines are links or lines between buses of the same carrier at different
locations. Most of the design decisions described in this section are choices
about which buses exist and which links connect them.

**Endogenous and exogenous.** The optimisation decides capacities and
dispatch of all technologies that are declared extendable, and with them the
fuel mix of most conversion steps. Some transformations are fixed exogenously
because they are driven by policy or consumer choice rather than by system
cost, for example the share of electric vehicles in land transport, the
production route of steel, or the share of buildings on district heating. These
exogenous shares are given per planning year, so that they can follow a
transition path. Each page states which quantities are optimised and which are
prescribed.

**Energy imports.** Europe is not closed. Besides fossil fuels, which enter
at a price, the model can import green hydrogen and its derivatives, ammonia,
methanol, synthetic methane and oil, as well as solid biomass, from outside
the modelled countries at fixed prices and with an optional cap on the total.
Whether imports are cheaper than domestic production is then part of the
optimisation. The carrier pages note where imports connect.

**Data-driven.** All inputs come from open datasets, from the grid topology and
weather reanalysis to energy balances and industrial site locations. The
[data sources](../data_sources.md) page lists them with their licences; the
rule documentation describes how each dataset is processed, and the
[configuration](../configuration.md) reference lists the settings. This
section explains what is represented and why, not how it is computed.
