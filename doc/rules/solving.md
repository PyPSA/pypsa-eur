<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Composing and Solving Networks

One rule assembles the network for each planning horizon from the clustered
electricity network and the prepared sector data, and one rule solves it. Both
run the same way for electricity-only and sector-coupled studies and for all
[foresight](../design/foresight.md) modes; the foresight mode decides which
inputs the composition step needs. Dispatch-only analyses on an already solved
network fix the expanded capacities and re-solve operation, optionally in a
rolling horizon.

## Temporal aggregation

{{ rules("time_aggregation") }}

## Composing

{{ rules("compose_network") }}

## Solving

{{ rules("solve_network", "solve_operations_network") }}

## Library modules

These modules are not rules of their own. Their functions are imported and
called by [compose_network][] in the order the foresight mode requires.

### `add_electricity`

::: add_electricity

### `prepare_network`

::: prepare_network

### `prepare_sector_network`

::: prepare_sector_network

### `add_existing_baseyear`

::: add_existing_baseyear

### `add_brownfield`

::: add_brownfield

### `prepare_perfect_foresight`

::: prepare_perfect_foresight
