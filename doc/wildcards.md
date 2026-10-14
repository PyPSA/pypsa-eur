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
It can take the values `onwind`, `offwind-ac`, `offwind-dc`, `offwind-float`, `solar` and `solar-hsat` but **not** `hydro`
(since hydroelectric plant profiles are created by [build_hydro_profile][]).

## The `{horizon}` wildcard {#planning_horizons}

The `{horizon}` wildcard is the planning year of a network, e.g. 2030, 2040 or 2050.
It takes the values listed under `planning_horizons` and appears in all
horizon-specific files such as `composed_{horizon}.nc` and `solved_{horizon}.nc`.

## The `{run}` wildcard {#run}

The `{run}` wildcard appears in `resources/` and `results/` paths when
`run: scenarios: enable:` is `true`. It takes the scenario names defined in the
scenario file (see [run](configuration.md#run_cf)).

## The `{country}` wildcard {#country}

The `{country}` wildcard takes two-letter country codes from `countries:`. It is
used by rules that retrieve or process data per country, such as
`retrieve_osm_data_raw`, [clean_osm_data][] and [build_osm_boundaries][].

## The `{year}` wildcard {#year}

The `{year}` wildcard takes the weather years covered by `snapshots:`. It is
used by rules that retrieve yearly data files, such as
`retrieve_seawater_temperature` and `retrieve_hera_data`.

## The `{carrier}` wildcard {#carrier}

The `{carrier}` wildcard selects the bus carrier for which a map is plotted, for
example `AC`, `H2` or `heat`. It is used by [plot_balance_map][],
`plot_balance_map_interactive` and [plot_heat_source_map][]. The carriers to
plot are set under `plotting:`.

## The `{kind}` wildcard {#kind}

The `{kind}` wildcard takes the values `energy` and `heat` and selects which
totals are distributed to regions by [build_population_weighted_energy_totals][].
