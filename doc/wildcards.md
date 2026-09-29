<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Wildcards {#wildcards}

!!! note
    If you are migrating from an earlier release, see [migration](migration.md) for
    guidance on translating legacy wildcards and targets.

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

## The `{horizon}` wildcard {#planning_horizons}

The `{horizon}` wildcard is used for multi-period optimization studies. It
takes years as values, e.g. 2030, 2040, 2050.
