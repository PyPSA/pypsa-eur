<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Solving Networks

After generating and clustering the networks, [compose_network][] produces
`networks/composed_{horizon}.nc` for each planning horizon. These files are
then passed to the single [solve_network][] rule, which runs
`scripts/solve_network.py` regardless of whether the study is electricity-only
or sector-coupled. Dispatch-only analyses on an already solved network are
available through [solve_operations_network][], which fixes the expanded
capacities and re-solves operation, optionally in a rolling horizon manner via
`solving.operations`.

## Rule `solve_network` {#solve}

::: solve_network

## Rule `solve_operations_network` {#solve_operations}

::: solve_operations_network
