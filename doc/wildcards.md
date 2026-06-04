<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Wildcards {#wildcards}

It is easy to run PyPSA-Eur for multiple scenarios using the wildcards feature of `snakemake`.
Wildcards allow to generalise a rule to produce all files that follow a regular expression pattern
which e.g. defines one particular scenario. One can think of a wildcard as a parameter that shows
up in the input/output file names of the `Snakefile` and thereby determines which rules to run,
what data to retrieve and what files to produce.

!!! note
    Detailed explanations of how wildcards work in `snakemake` can be found in the
    [relevant section of the documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards).

## The `{cutout}` wildcard {#cutout_wc}

The `{cutout}` wildcard facilitates running the rule [build_cutout][]
for all cutout configurations specified under `atlite: cutouts:`.
These cutouts will be stored in a folder specified by `{cutout}`.

## The `{technology}` wildcard {#technology}

The `{technology}` wildcard specifies for which renewable energy technology to produce availability time
series and potentials using the rule [build_renewable_profiles][].
It can take the values `onwind`, `offwind-ac`, `offwind-dc`, `offwind-float`, and `solar` but **not** `hydro`
(since hydroelectric plant profiles are created by a different rule)

## The `{clusters}` wildcard {#clusters}

The `{clusters}` wildcard specifies the number of buses a detailed
network model should be reduced to in the rule [cluster_network][].
The number of clusters must be lower than the total number of nodes
and higher than the number of countries. However, a country counts twice if
it has two asynchronous subnetworks (e.g. Denmark or Italy).

## The `{opts}` wildcard {#opts}

The `{opts}` wildcard is used for electricity-only studies. It triggers
optional constraints, which are activated in either [prepare_network][] or
the [solve_network][] step. It may hold multiple triggers separated by `-`,
i.e. `Co2L-3h` contains the `Co2L` trigger and the `3h` switch. There are
currently:

| Trigger | Description |
|---------|-------------|
| `Co2L` | Add an pointwise or global CO2 limit constraint |
| `Co2L{p}` | Set CO2 limit to `p` times the 1990 CO2 emissions |
| `Ep` | Add an pointwise or global emission pricing constraint |
| `Ep{p}` | Set the emission price to `p` EUR/tCO2 |
| `{n}h` | Temporal averaging/down-sampling by `n` hours |
| `{n}seg` | Apply `n` time segments using the tsam package |

!!! note
    The wildcard options are processed in ``scripts/prepare_network.py`` and ``scripts/solve_network.py``.

## The `{sector_opts}` wildcard {#sector_opts}

!!! warning
    More comprehensive documentation for this wildcard will be added soon.
    To really understand the options here, look in scripts/prepare_sector_network.py

The `{sector_opts}` wildcard is only used for sector-coupling studies.

| Trigger | Description |
|---------|-------------|
| `cb{x}ex{y}` | Set a carbon budget of `x` Gt CO2 with exponential decay starting at linear growth rate `y` |
| `Co2L{p}` | Set CO2 limit to `p` times the 1990 CO2 emissions |
| `{n}h` | Temporal averaging/down-sampling by `n` hours |
| `{n}seg` | Apply `n` time segments using the tsam package |

!!! note
    The wildcard options for sector-coupled studies are processed in ``scripts/prepare_sector_network.py``.
    For more details, refer to the source code.

## The `{planning_horizons}` wildcard {#planning_horizons}

!!! warning
    More comprehensive documentation for this wildcard will be added soon.

The `{planning_horizons}` wildcard is only used for sector-coupling studies.
It takes years as values, e.g. 2020, 2030, 2040, 2050.
