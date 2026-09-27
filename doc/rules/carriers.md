<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Biomass, Gas and Carbon

The rules on this page prepare the regional potentials and infrastructure of
the carriers that are not electricity: biomass potentials and transport costs,
the gas network with its entry points, salt caverns for hydrogen storage and
geological CO~2~ sequestration potentials. See [Biomass](../design/biomass.md),
[Methane](../design/methane.md), [Hydrogen](../design/hydrogen.md) and [Carbon
management](../design/carbon.md) in the design section.

## Biomass

{{ rules("build_biomass_potentials", "build_biomass_transport_costs") }}

## Gas network and hydrogen storage

{{ rules("build_gas_network", "build_gas_input_locations", "cluster_gas_network", "build_salt_cavern_potentials") }}

## Carbon sequestration

{{ rules("build_co2_sequestration_potentials", "build_clustered_co2_sequestration_potentials") }}
