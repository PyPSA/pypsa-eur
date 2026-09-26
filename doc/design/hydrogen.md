<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Hydrogen {#design-hydrogen}

{{ scope(electricity="partial") }}

Hydrogen links the electricity system to industry, transport and the synthetic
fuels [@neumannPotentialRole2023; @staffellRoleHydrogen2019a]. Every region
has a hydrogen bus.

=== "Sector-coupled"

    Hydrogen is a full carrier with demands, production routes, pipelines
    and storage as described on this page.

=== "Electricity-only"

    Hydrogen exists only as an electricity storage option: a storage unit
    with a fixed energy-to-power ratio, or a store with electrolyser and fuel
    cell links.

## Demand

Hydrogen is consumed as a feedstock in industry, for direct reduction of iron
and for ammonia synthesis, and as an input to synthetic methane, methanol and
Fischer-Tropsch fuels. In transport it fuels fuel cell vehicles and,
optionally, ships. These demands are either fixed by the industry and
transport pathways or follow from the optimisation of the synthetic fuel
routes. Hydrogen can also be re-electrified in fuel cells or hydrogen
turbines, which lets it act as long-duration storage for the power system.

## Supply

Several production routes compete:

- **Electrolysis** splits water with electricity. Its waste heat can be fed
  into district heating.
- **Steam methane reforming** converts methane, which may itself be fossil,
  biogenic or synthetic. A variant with carbon capture reduces its emissions.
- **Biomass gasification** produces hydrogen from solid biomass, with or
  without carbon capture, see [Biomass](biomass.md).
- **Ammonia cracking** and **methanol reforming** recover hydrogen from
  carriers that are easier to ship, which matters when imports are enabled.

Hydrogen and its derivatives can also be imported from outside Europe at a
fixed price, with an optional cap on the total, see
[Model overview](overview.md#design-overview). The split between routes and
the installed capacities are results of the optimisation.

## Storage

Hydrogen can be stored in underground salt caverns where the geology allows
[@caglayanTechnicalPotential2020] and, far more expensively, in steel tanks
where it does not. Cavern potentials are regional and can be restricted to
onshore sites or sites near the shore to avoid brine disposal at sea. The
constraint is usually where caverns exist, not how much they can hold.

## Transport

Hydrogen moves through pipelines. New pipelines can be built along corridors
where an electricity or gas connection exists today. Optionally, existing gas
pipelines can be retrofitted [@gasforclimateEuropeanHydrogen2020], which is
cheaper but tied to the gas network: for every unit of gas capacity taken out
of service, a fixed share becomes available for hydrogen on the same route.
If the gas network is not resolved regionally, the retrofit option still
exists as a potential per corridor. The electricity demand of compression per
distance can be represented. Whether a European hydrogen backbone emerges,
and from which pipelines, is a result [@neumannPotentialRole2023].

## Further reading

- Rules: [build_salt_cavern_potentials][], [build_gas_network][],
  [cluster_gas_network][]
- Configuration: [sector](../configuration.md#sector_cf), in particular
  `H2_network`, `H2_retrofit`, `hydrogen_underground_storage`, `SMR` and
  `imports`
