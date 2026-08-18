<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Migration Guide {#migration}

This guide helps you move an existing PyPSA-Eur setup to the streamlined
workflow ([#1838](https://github.com/PyPSA/pypsa-eur/pull/1838)). Other changes
that were released at the same time are listed in the
[release notes](release_notes.md).

## What changed

Scenario information used to live in filenames, such as
`elec_s_37_lv1.25_3H_2030.nc`. It now lives in the configuration file. The
workflow runs through five stages with fixed, readable names:

```
base → simplified → clustered → composed_{horizon} → solved_{horizon}
```

Only one wildcard is left: `{horizon}`, the planning year.

All three foresight modes use these same stages. They differ only in how they
loop over the horizons:

- `overnight` composes and solves a single horizon.
- `myopic` composes and solves one horizon after the other. Each `compose_network`
  reads the previous horizon's *solved* network.
- `perfect` composes every horizon in sequence, each one reading the previous
  horizon's *composed* network, and then solves them all together.

See [foresight](foresight.md) for the details.

## How to use this guide

Everyone upgrading should read [Before you start](#before-you-start) and
[Rewrite your config](#rewrite-your-config). If you maintain a fork or soft-fork,
also read [Migrating a fork](#fork-migration): the scripts and rules were
restructured, so you will get merge conflicts.

!!! danger "An outdated config does not raise an error"
    This is the one thing to remember. Renamed keys are not rejected, they are
    ignored. Your old config still validates, your run still succeeds — but it
    models something you did not ask for, because every ignored key silently falls
    back to its default. For example:

    - a leftover `scenario:` block leaves you with one horizon (2050) and 50 clusters,
    - a leftover flat `co2_budget:` mapping leaves you on the default CO₂ trajectory,
    - a leftover `resolution_elec:` leaves you at hourly resolution, which is slow,
    - a config without a `sector:` block builds a full sector-coupled model.

    Work through [Rewrite your config](#rewrite-your-config) key by key
    instead of waiting for an error message. There will not be one. If your results
    look odd afterwards, check [Symptoms and causes](#symptoms).

## Before you start {#before-you-start}

**Your retrieved data and cutouts stay valid.** Nothing about the input data or the
weather cutouts changed, so there is no need to download or cut anything again.
There are two exceptions, both files that you supply yourself and that were keyed by
the old `{clusters}` wildcard. If you use `clustering.mode: custom_busmap` or
`custom_busshapes`, rename your files to
`data/busmaps/simplified_{n_clusters}_{base_network}.csv` and
`data/busshapes/simplified_{n_clusters}_{base_network}.geojson`.

**Your `resources/` and `results/` directories are outdated.** Almost every filename
changed, so Snakemake does not recognise the existing files and rebuilds everything
from `base.nc` onwards. The old files are not deleted, they are just left behind.
Keep them for comparison or delete them once you are happy with the new run, but
plan for one full rebuild. If you only want the new summaries and do not want to
solve again, see [reusing solved networks](#reusing-solved-networks).

**Do not expect identical numbers.** A few changes go beyond renaming and change the
model itself. Read [Behaviour changes](#behaviour-changes) before you compare a new
run to an old one.

**Your analysis scripts need attention.** The summary CSVs changed shape, and the
per-horizon `csvs/individual/*.csv` files are gone. See
[Summary outputs](#summaries).

Then migrate in this order:

1. [Declare your model type](#model-type) — electricity-only or sector-coupled.
2. [The `scenario` block](#scenario-block) — horizons and cluster count.
3. [`opts` and `sector_opts`](#opts-translation) — every option token becomes a config key.
4. [CO₂ constraints](#co2-budget) — restructured, and the units changed.
5. [Temporal resolution](#temporal) — one setting instead of two.
6. [Other renamed keys](#other-keys) — and the new keys you should decide on.
7. [File paths](#files) — in your own scripts.
8. [Snakemake targets](#targets) — in your commands.

Steps 1 to 6 are config edits, and they end with a
[complete before-and-after example](#example). After that, read
[Behaviour changes](#behaviour-changes) and [Summary outputs](#summaries) before
you compare results or run your analysis scripts.

!!! tip "Learn from a working config"
    If you prefer to start from a complete file,
    `config/test/config.electricity.yaml` is the clearest electricity-only example
    in the repository, and `config/test/config.myopic.yaml` shows the multi-horizon
    sector-coupled form.

## Rewrite your config {#rewrite-your-config}

### Declare your model type {#model-type}

Whether you build an electricity-only or a sector-coupled model used to depend on
the target you typed: `solve_elec_networks` or `solve_sector_networks`. It is now a
config key, `sector.enabled`. Both model types share the same rules and the same
default target.

!!! danger "Electricity-only models must opt out explicitly"
    `sector.enabled` defaults to `true`, and an old electricity-only config has no
    `sector:` block at all. Such a config therefore validates and then builds and
    solves a full sector-coupled model, with all the data retrieval and runtime that
    this implies. Say it explicitly:

    ```yaml
    sector:
      enabled: false
    ```

### The `scenario` block is gone {#scenario-block}

The `scenario:` section used to list the values of the `{clusters}`, `{opts}`,
`{sector_opts}` and `{planning_horizons}` wildcards. All four wildcards were
removed, so the block is not read any more. Move its content as follows:

| Old | New |
|--------|-----|
| `scenario.planning_horizons` | `planning_horizons` (top level) |
| `scenario.clusters` | `clustering.cluster_network.n_clusters` |
| `scenario.opts` | explicit config keys, see [below](#opts-translation) |
| `scenario.sector_opts` | explicit config keys, see [below](#opts-translation) |

`planning_horizons` is now a top-level list of years, and it defines the values of
the `{horizon}` wildcard. A single year may be written as a scalar
(`planning_horizons: 2030`); it is turned into a list for you.

`n_clusters` takes one integer, not a list, because the wildcard sweep is gone. To
compare several cluster counts, write one scenario per count under `run.scenarios`
(see [configuration](configuration.md)). For administrative clustering, set
`clustering.mode: administrative` instead of the former `clusters: [adm]`.

!!! warning "Overnight runs need exactly one horizon"
    With `foresight: overnight`, `planning_horizons` must contain exactly one year.
    Otherwise the workflow stops while building the DAG, before any job starts, with
    `Overnight optimization can only be run for a single planning horizon.` An
    overnight setup could previously list several years and solve each one
    independently. Use `run.scenarios` for that now.

!!! warning "What a scenario may not change"
    The collection targets and the default target are built from the *base* config.
    Therefore `foresight` and `planning_horizons` must be identical in every
    scenario. Set them at the top level, and run different foresight modes or
    horizon sets as separate workflows with their own `run.name`. Everything that
    rules read through `config_provider` may vary per scenario, including
    `clustering.cluster_network.n_clusters`, `co2_budget`, `sector` and `solving`.

### Translating `opts` and `sector_opts` {#opts-translation}

Every option token has a config equivalent. The translation layer
(`_helpers.update_config_from_wildcards`) was deleted, so leftover option strings
have no effect at all.

`opts` tokens:

| Old token | New config |
|--------------|-----------|
| `Co2L<x>` | `co2_budget.upper` (see [CO₂](#co2-budget)) |
| `<n>h` | `clustering.temporal.averaging: <n>` |
| `<n>seg` | `clustering.temporal.segmentation: <n>` |
| `CH4L<x>` | `electricity.gaslimit_enable: true` + `electricity.gaslimit` |
| `Ep<x>` | `costs.emission_prices.enable: true` + `costs.emission_prices.co2` |
| `Ept` | `costs.emission_prices.dynamic: true` |
| `ATK` / `ATKc` | `electricity.autarky.enable` / `electricity.autarky.by_country` |
| `lv<x>` / `lc<x>` | `electricity.transmission_limit` (e.g. `vopt`, `v1.25`, `c1.25`) |
| `<carrier>+<comp>+<p\|e\|c\|m><factor>` | `adjustments.electricity.factor` |

`sector_opts` tokens:

| Old token | New config |
|--------------|-----------|
| `T`, `H`, `B`, `I`, `A` | `sector.transport`, `heating`, `biomass`, `industry`, `agriculture` |
| `CCL`, `EQ`, `BAU`, `SAFE` | `solving.constraints.CCL`, `EQ`, `BAU`, `SAFE` |
| `<n>h` / `<n>sn` / `<n>seg` | `clustering.temporal.averaging` / `representative` / `segmentation` |
| `decentral` | `sector.electricity_transmission_grid: false` |
| `noH2network` | `sector.H2_network: false` |
| `nowasteheat` | the `sector.use_*_waste_heat` flags, set to `0` |
| `nodistrict` | `sector.district_heating.progress: {<year>: 0, …}` |
| `dist<f>` | `sector.electricity_distribution_grid: true` + `sector.electricity_distribution_grid_cost_factor` |
| `biomasstransport` | `sector.biomass_transport: true` |
| `linemaxext<n>` | `lines.s_nom_max_extension` and `links.p_nom_max_extension` |
| `sdr<x>` | `costs.social_discountrate` |
| `seq<x>` | `sector.co2_sequestration_potential` |
| `<carrier>+<comp>+<p\|e\|c\|m><factor>` | `adjustments.sector.factor` |
| `Co2L<x>`, `cb<x>` | `co2_budget.upper` |
| `CF+<key>+<subkey>+<value>` | set that key directly in the config |
| `cb<x>ex`, `cb<x>be` | no replacement: `build_carbon_budget` is not called any more, so write the resulting per-year values into `co2_budget.upper` yourself |

Three tokens did a unit conversion for you that you now have to do yourself:

- `linemaxext<n>` was given in GW and multiplied by `1e3`, so `linemaxext20`
  becomes `20000` (MW).
- `CH4L<x>` was multiplied by `1e6`, so `CH4L200` becomes `gaslimit: 200000000`.
- `sdr<x>` was a percentage divided by `100`, so `sdr2` becomes
  `social_discountrate: 0.02`.

!!! note
    `sector.district_heating.progress` has to stay a year-to-value mapping; a bare
    `0` is rejected by the schema. Config files are deep-merged, so list every
    horizon you optimise, for example `progress: {2030: 0, 2040: 0, 2050: 0}`.

### CO₂ constraints {#co2-budget}

CO₂ is now handled entirely through `co2_budget`. Three separate mechanisms feed
into it:

| Old | New |
|--------|-----|
| `co2_budget: {<year>: <fraction>}` (flat mapping) | `co2_budget.upper: {<year>: <fraction>}` |
| `electricity.co2limit_enable` / `co2limit` / `co2base` | `co2_budget.upper` with `relative: false` |
| `energy.emissions` | `co2_budget.emissions_scope` |

The new block looks like this:

```yaml
co2_budget:
  emissions_scope: CO2     # was energy.emissions
  relative: true           # upper/lower are fractions of 1990 emissions
  upper:
    2030: 0.45
    2050: 0.0
  lower:                   # optional floor, unset by default
```

There are three traps here.

**The units of absolute limits changed.** With `relative: false`, `upper` and
`lower` are in **Gt CO₂ per year**, while the old `electricity.co2limit` was in
tonnes. Divide by `1e9`, so `co2limit: 77500000.0` becomes `upper: 0.0775`. If you
copy the old number unchanged, your limit is a billion times too loose.

**Electricity-only runs are now capped by default.** Previously an
electricity-only model had no CO₂ constraint unless you set
`electricity.co2limit_enable: true`. Now the same `co2_budget` block applies to
electricity-only runs, and the shipped default ends at `2050: 0.0`. An
electricity-only run with `planning_horizons: [2050]` is therefore capped at zero
emissions and may be infeasible. To get the old, unconstrained behaviour, set
`co2_budget: {upper: null}`.

**Partial overrides are merged, not replaced.** Config files are deep-merged, so
writing only `co2_budget: {upper: {2030: 0.4}}` keeps the other six default years.
Set the years you do not want to `null`.

A mapping and a scalar mean two different things. A mapping constrains each listed
year on its own. A scalar is a single budget across all periods; under perfect
foresight it becomes one cumulative constraint attached to the final horizon
instead of a per-period cap. See [foresight](foresight.md) for the exact semantics
in each mode.

!!! warning "List every horizon you optimise"
    A horizon that is missing from `co2_budget.upper` gets **no CO₂ constraint at
    all**. The workflow only logs an info line and continues. Earlier releases
    interpolated between the neighbouring years instead. So make sure that every
    entry of `planning_horizons` also appears in `upper` (and in `lower`, if you use
    it).

### Temporal resolution {#temporal}

`clustering.temporal.resolution_elec` and `resolution_sector` are replaced by three
integer options under `clustering.temporal`, of which you may set only one:

- `averaging: <n>` averages over `n` hours,
- `segmentation: <n>` aggregates the year into `n` segments using `tsam`,
- `representative: <n>` keeps every `n`-th snapshot.

The same setting now applies to electricity-only and to sector-coupled runs. If you
used two different values before, you have to pick one. Setting more than one option
raises `clustering.temporal: only one of averaging, segmentation and representative
may be set`.

| Old string | New |
|---------------|-----|
| `resolution_elec: 24h` / `resolution_sector: 24H` | `averaging: 24` |
| `resolution_sector: 8760h` | `averaging: 8760` |
| `resolution_*: 4380seg` | `segmentation: 4380` |
| `resolution_*: 3sn` | `representative: 3` |
| `resolution_*: false` | leave all three `false` |

!!! note
    `averaging` and `representative` reproduce the previous numbers exactly.
    `segmentation` does not: the segments are now derived from the renewable profile,
    electricity demand, heat demand and solar thermal resources rather than from the
    time series of the fully prepared network, so segment boundaries and weights
    differ.

### Other renamed keys {#other-keys}

| Old | New |
|--------|-----|
| `lines.max_extension` | `lines.s_nom_max_extension` |
| `links.max_extension` | `links.p_nom_max_extension` (default 30000 MW, now applies to sector runs too) |
| `solving.options.rolling_horizon` (+ `horizon`, `overlap`) | `solving.operations.rolling_horizon` (+ `horizon`, `overlap`) |

The rolling-horizon settings moved into their own block because they configure
`solve_operations_network`. The key `solving.options.rolling_horizon` still exists,
but it now only affects the capacity-expansion rule `solve_network`. This way you
can run capacity expansion first and rolling-horizon dispatch afterwards.

Two keys are new and worth a decision:

- `existing_capacities.phase_outs` holds the national policy phase-outs that used to
  be hardcoded, as a list of `{carriers, countries, year}` rules. They cap the
  lifetimes of conventional assets under **perfect foresight only**: generators in
  electricity-only runs, links in sector-coupled ones. The defaults reproduce the
  previous behaviour. Careful: lists are replaced as a whole and not merged, so your
  own list replaces all defaults, including the ones you did not want to drop.
- `run.default_target_rule` selects the collect rule that `snakemake` builds when
  you name no target. It defaults to `all`.

Finally, `solving.check_objective.expected_value` now accepts either a single value
or a `{<horizon>: value}` mapping, so every `solved_{horizon}.nc` of a multi-horizon
run is checked against its own benchmark. The default tolerances were tightened
(`atol` 1e6 → 1e4, `rtol` 0.01 → 0.001).

### Complete example {#example}

**Before:**

```yaml
foresight: myopic

scenario:
  clusters: [37]
  opts: [Co2L0.45-24h]
  sector_opts: [""]
  planning_horizons: [2030, 2040]

energy:
  emissions: CO2

clustering:
  temporal:
    resolution_elec: 24h
    resolution_sector: 24h

lines:
  max_extension: 20000
links:
  max_extension: 20000
```

**After:**

```yaml
foresight: myopic
planning_horizons: [2030, 2040]

clustering:
  cluster_network:
    n_clusters: 37
  temporal:
    averaging: 24

co2_budget:
  emissions_scope: CO2
  relative: true
  upper:
    2030: 0.45
    2040: 0.45

lines:
  s_nom_max_extension: 20000
links:
  p_nom_max_extension: 20000

sector:
  enabled: true
```

## Behaviour changes {#behaviour-changes}

The changes below are not renames. Even with an equivalent config the model itself
changes, so do not expect to reproduce old numbers unless you act on them. Each item
says how to get the old behaviour back, where that is possible.

**Overnight runs now use per-horizon costs.** Every layer is costed with
`costs_{horizon}_processed.csv`, and `costs.year` no longer influences composition.
Previously an overnight run used `costs.year` (default 2050) regardless of the
horizon, and in myopic and perfect runs the electricity layer used `costs.year` while
the sector layer used the horizon. Wherever those two years differed, all costs,
efficiencies and lifetimes change.
*Old behaviour:* set `planning_horizons: [<old costs.year>]`.

**Electricity-only runs are CO₂-capped by default**, and absolute caps are given in
Gt instead of tonnes. A default electricity-only run therefore changes from
unconstrained to zero-emission.
*Old behaviour:* `co2_budget: {upper: null, lower: null}`, or divide your old tonne
value by `1e9`.

**Relative budgets are always measured against 1990.**
`energy.base_emissions_year` no longer influences the CO₂ baseline; it still controls
`build_co2_totals`. For electricity-only runs the baseline still covers every sector
that is flagged in `sector:`, which makes a relative budget almost non-binding.
*Old behaviour:* nothing to do for sector runs with a 1990 baseline. Otherwise
rescale `upper` by hand or use `relative: false`. For electricity-only runs, switch
the `sector.transport`, `heating`, `industry` and `agriculture` flags off.

**A missing horizon drops the CO₂ constraint** instead of interpolating. That
horizon is then unconstrained, and in myopic mode so is every later one.
*Old behaviour:* list every horizon explicitly in `co2_budget.upper`.

**Sector-coupled HVDC extension headroom rises from 20 to 30 GW.** In addition, a
finite `lines.s_nom_max` or `links.p_nom_max` now binds in sector runs, where it used
to be ignored.
*Old behaviour:* `links.p_nom_max_extension: 20000`, and `lines.s_nom_max: .inf` /
`links.p_nom_max: .inf`.

**Both `adjustments` blocks are applied**, always with the current horizon, and after
transmission costs are set. Sector runs now also pick up `adjustments.electricity`,
electricity-only runs also pick up `adjustments.sector`, and Line/Link
`capital_cost` adjustments now take effect.
*Old behaviour:* keep entries only in the block that your model type used before, and
drop Line/Link `capital_cost` adjustments.

**`segmentation` is computed from the profile and demand resources** rather than from
the prepared network's time series, so segment boundaries and weights differ, and with
them capacity factors and the objective.
*Old behaviour:* use `averaging` or `representative`; both are numerically unchanged.

**Electricity-only myopic and perfect foresight now run.** An electricity-only config
used to ignore `foresight` and always solved a single period. Now `foresight: myopic`
or `perfect` gives you a genuinely different model, with an existing fleet and with
capacity accumulating across horizons.
*Old behaviour:* use `foresight: overnight` with a single horizon.

### If you use perfect foresight {#perfect}

Sector-coupled perfect foresight works again; earlier releases stopped with `PyPSA
versions >=1.0 are not supported for perfect foresight`. The assembly was rewritten
rather than translated line by line, so treat perfect-foresight results as new
results rather than comparable ones. Three points in particular:

- **AC lines stay lines.** Previously every AC line was replaced by a DC link with a
  flat 3%/1000 km efficiency and one vintage per investment period. Now the lines
  survive as a single extendable set that is shared across all periods, with PyPSA's
  resistance-based loss model (`solving.options.transmission_losses`) and the
  linearised power flow. Transmission is thus one decision for the whole horizon
  instead of a staged build-out, which changes the objective, the line capacities and
  the cross-border flows. You cannot recover the old behaviour from the config; the
  conversion code no longer exists.
- **Per-horizon maps are not produced.** In this mode only the summary CSVs and the
  `graphs/*.svg` summary plots are created.
- **The summary CSVs changed**, both in content and in which files exist. See
  [Summary outputs](#summaries).

## Summary outputs and your analysis scripts {#summaries}

This is the change that is most likely to break your own post-processing.

The two-stage summary is gone. `make_summary` no longer writes per-horizon
`csvs/individual/<metric>_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv`
files, and `make_global_summary` and `make_cumulative_costs` were removed as separate
rules. A single `make_summary` writes `results/csvs/<metric>.csv` directly, built from
all solved networks at once via `pypsa.NetworkCollection`.

This has two consequences for your scripts:

- **`csvs/individual/` does not exist any more.** The planning horizons are now
  **columns** of the aggregated CSV. To read a single horizon, slice that column:

  ```python
  pd.read_csv("results/csvs/costs.csv", index_col=[0, 1, 2])["2030"]
  ```

  The column label is the horizon year and comes back from CSV as a string. How deep
  the index is depends on the metric: three levels for `costs`, four for
  `nodal_costs`, one for `metrics`.
- **The aggregated CSVs kept their paths but changed their header.** They used to
  carry a three-level `(cluster, opt, planning_horizon)` column MultiIndex. A reader
  that still passes `header=[0, 1, 2]` does not fail, it silently mis-parses the new
  single-row header.

Under perfect foresight the set of files changes as well. Summaries are now computed
from PyPSA `statistics` instead of the deleted `make_summary_perfect.py`, so the
reported costs are undiscounted and not split by period activity. Concretely:

- `csvs/supply.csv` and `csvs/supply_energy.csv` are replaced by
  `csvs/energy_balance.csv` and `csvs/nodal_energy_balance.csv`, which are not
  column-for-column equivalents.
- `csvs/price_statistics.csv` and `csvs/co2_emissions.csv` have no replacement. You
  can recover the CO₂ figures from the `co2` rows of `energy_balance.csv`.

If you run perfect foresight, also read
[If you use perfect foresight](#perfect).

One per-solve output is no longer written: the config dumps
`results/configs/config.*.yaml`. The effective config is still stored inside each
solved network as `n.meta`. `solving.options.store_model` still writes the linopy
model, now to `results/models/solved_{horizon}.nc`.

### Reusing solved networks from an old run {#reusing-solved-networks}

Since the filenames changed, Snakemake finds no solved network and would solve
everything again. If you only want the new summaries, rename your old solved networks
to `results/<run>/networks/solved_{horizon}.nc` and ask for the summary:

```console
$ snakemake results/<run>/csvs/costs.csv --configfile <configfile>
```

The rename is all it takes. Snakemake does not rebuild the upstream chain of an
output that already exists, so the missing `composed_*.nc` and `resources/` files are
not created. Note that `make_summary` reads every horizon of `planning_horizons` in
one job — only the last one under perfect foresight — so all of them have to be
there.

## Files and directories {#files}

The paths below leave out the run directory. It is inserted as `resources/<run>/…`,
`results/<run>/…`, `logs/<run>/…` and `benchmarks/<run>/…` when `run.name` is set:
as a literal directory, or as a `{run}` wildcard if `run.scenarios.enable` is also
true. A non-empty `run.prefix` adds one more level in front. With the shipped
defaults, where `run.prefix` and `run.name` are empty, the paths are exactly as
shown.

### Network stages

| Old | New |
|--------|-----|
| `resources/networks/base_s.nc` | `resources/networks/simplified.nc` |
| `resources/networks/base_s_{clusters}.nc` | `resources/networks/clustered.nc` |
| `resources/networks/base_s_{clusters}_elec.nc` | *removed*, folded into `composed_{horizon}.nc` |
| `resources/networks/base_s_{clusters}_elec_{opts}.nc` | `resources/networks/composed_{horizon}.nc` |
| `resources/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc` | `resources/networks/composed_{horizon}.nc` |
| `results/networks/base_s_{clusters}_elec_{opts}.nc` | `results/networks/solved_{horizon}.nc` (electricity-only) |
| `results/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc` | `results/networks/solved_{horizon}.nc` (sector-coupled) |
| `results/networks/base_s_{clusters}_elec_{opts}_op.nc` | `results/networks/operations_{horizon}.nc` (opt-in dispatch re-solve) |

### Regions, busmaps and other resources

| Old | New |
|--------|-----|
| `regions_onshore.geojson` / `regions_offshore.geojson` | `onshore_regions_base.geojson` / `offshore_regions_base.geojson` |
| `regions_onshore_base_s.geojson` / `regions_offshore_base_s.geojson` | `onshore_regions_simplified.geojson` / `offshore_regions_simplified.geojson` |
| `regions_onshore_base_s_{clusters}.geojson` / `regions_offshore_base_s_{clusters}.geojson` | `onshore_regions.geojson` / `offshore_regions.geojson` |
| `busmap_base_s.csv` | `busmap_simplify_network.csv` |
| `busmap_base_s_{clusters}.csv` | `busmap_cluster_network.csv` |
| `linemap_base_s_{clusters}.csv` | `linemap_cluster_network.csv` |
| — | `busmap.csv` (new: base buses → clustered buses, from `chain_busmaps`) |
| `powerplants_s_{clusters}.csv` | `powerplants.csv` |
| `pop_layout_base_s.csv` / `pop_layout_base_s_{clusters}.csv` | `pop_layout_simplified.csv` / `pop_layout.csv` |
| `electricity_demand_base_s.nc` | `electricity_demand_simplified.nc` |
| — | `electricity_demand.nc` (new: clustered demand, from `cluster_electricity_demand`) |
| `availability_matrix_{clusters}_{technology}.nc` | `availability_matrix_{technology}.nc` |
| `profile_{clusters}_{technology}.nc` | `profile_{technology}.nc` |
| `regions_by_class_{clusters}_{technology}.geojson` | `regions_by_class_{technology}.geojson` |
| `solar_rooftop_potentials_s_{clusters}.csv` | `solar_rooftop_potentials.csv` |
| `gas_input_locations_s_{clusters}.geojson` | `gas_input_locations.geojson` |
| `costs_{planning_horizons}[_processed].csv` | `costs_{horizon}[_processed].csv` |
| `data/busmaps/base_s_{clusters}_{base_network}.csv` | `data/busmaps/simplified_{n_clusters}_{base_network}.csv` |
| `data/busshapes/base_s_{clusters}_{base_network}.geojson` | `data/busshapes/simplified_{n_clusters}_{base_network}.geojson` |

Watch out for the flipped word order in the region files: `regions_onshore_*` became
`onshore_regions*`.

Beyond this table, the rule for sector resources is simple: `_s_{clusters}` is
dropped and `{planning_horizons}` becomes `{horizon}`. For example
`biomass_potentials_s_{clusters}_{planning_horizons}.csv` becomes
`biomass_potentials_{horizon}.csv`, and
`cop_profiles_base_s_{clusters}_{planning_horizons}.nc` becomes
`cop_profiles_{horizon}.nc`.

### Maps and graphs

| Old | New |
|--------|-----|
| `resources/maps/power-network.pdf` | `resources/maps/base_network.pdf` |
| `resources/maps/power-network-s-{clusters}.pdf` | `resources/maps/clustered_network.pdf` |
| `results/maps/static/base_s_…-costs-all_{planning_horizons}.pdf` | `results/maps/static/power_network_{horizon}.pdf` |
| `results/maps/static/base_s_…-h2_network_{planning_horizons}.pdf` | `results/maps/static/h2_network_{horizon}.pdf` |
| `results/maps/interactive/base_s_…_{planning_horizons}-balance_map_{carrier}.html` | `results/maps/interactive/balance_map_{carrier}_{horizon}.html` |
| `results/graphs/costs.pdf` | `results/graphs/costs.svg` |
| `results/graphs/energy.pdf` | `results/graphs/energy.svg` |
| `results/graphs/balances-energy.pdf` | `results/graphs/balances_energy.svg` (note the underscore) |
| `results/graphs/cop_profiles_s_{clusters}_{planning_horizons}.html` | `results/graphs/cop_profiles_{horizon}.html` |
| `results/graphics/balance_timeseries/s_{clusters}_…_{planning_horizons}/` | `results/graphs/balance_timeseries_{horizon}/` |
| `results/graphics/heatmap_timeseries/s_{clusters}_…_{planning_horizons}/` | `results/graphs/heatmap_timeseries_{horizon}/` |
| `results/graphics/interactive_bus_balance/s_{clusters}_…_{planning_horizons}/` | `results/graphics/interactive_bus_balance/{horizon}/` |

Note that `plot_summary` writes SVG now, not PDF.

## Snakemake targets {#targets}

Every foresight mode and both model types support the default target, so you can
always run the whole workflow with:

```console
$ snakemake --configfile <configfile>
```

The collection rules were unified:

| Old rule | New rule |
|-------------|----------|
| `cluster_networks` | `cluster_networks` (now `networks/clustered.nc` and `busmap.csv`) |
| `prepare_elec_networks`, `prepare_sector_networks` | `compose_networks` |
| `solve_elec_networks`, `solve_sector_networks`, `solve_sector_networks_perfect` | `solve_networks` |
| `plot_power_networks_clustered` | `plot_power_networks` |
| — | `solve_operations_networks` (new) |

`solve_networks` only asks for the solved network of the last horizon; the earlier
horizons follow from the dependency chain. `compose_networks` expands over the full
`planning_horizons` list, and so does `process_costs` — except under
`foresight: overnight`, where it still derives a single file from `costs.year`.
Composition always reads `costs_{horizon}_processed.csv`, so on an overnight run this
collect target can build a cost file that the run itself never uses.

If you target files directly:

| Old | New |
|--------|-----|
| `snakemake resources/networks/base_s_37.nc` | `snakemake resources/networks/clustered.nc` |
| `snakemake resources/networks/base_s_37_elec_Co2L-24h.nc` | `snakemake resources/networks/composed_2030.nc` |
| `snakemake results/networks/base_s_37_elec_Co2L-24h.nc` | `snakemake results/networks/solved_2030.nc` |
| `snakemake results/networks/base_s_37_Co2L-24h_Co2L0-1H-T-H-B-I_2030.nc` | `snakemake results/networks/solved_2030.nc` |

The `**config["scenario"]` expansion pattern is gone. Collection rules take the
`{horizon}` wildcard from `config["planning_horizons"]`.

## Symptoms and causes {#symptoms}

| Symptom | Likely cause |
|---------|--------------|
| Run builds hydrogen and industry inputs, or takes far longer than the old electricity-only run | `sector.enabled` defaults to `true`; set it to `false`. Note that the default target still builds heat-pump COP profiles even for electricity-only runs — use `run.default_target_rule: solve_networks` to avoid them |
| Runs only 2050, or only one horizon | `planning_horizons` is still inside `scenario:`; move it to the top level |
| Runs at 50 clusters regardless of the config | `scenario.clusters` was not moved to `clustering.cluster_network.n_clusters` |
| Model solves at hourly resolution and takes forever | `resolution_elec`/`resolution_sector` was not translated to `averaging`/`segmentation`/`representative` |
| Electricity-only run infeasible right after upgrading | the default `co2_budget` now applies; set `co2_budget: {upper: null}` or a realistic cap |
| CO₂ cap seems to do nothing | the absolute value is still in tonnes — divide by `1e9`; or the horizon is missing from `co2_budget.upper` |
| Grid expansion larger than expected | `max_extension` was not renamed, so the default 20000/30000 MW applies |
| `MissingInputException` on a `data/busmaps/…` or `data/busshapes/…` path | rename your custom busmap or bus shapes file, see [before you start](#before-you-start) |
| `ValueError: Overnight optimization can only be run for a single planning horizon.` | `foresight: overnight` with several horizons; use one horizon or `run.scenarios` |
| `ValidationError` on `clustering.cluster_network.n_clusters` | `n_clusters` must be an integer; `all` and `null` are rejected by the schema |
| `ValidationError` on `sector.district_heating.progress` | it has to stay a year-to-value mapping, not a scalar |
| `ValueError: clustering.temporal: only one of averaging, segmentation and representative may be set, got …` | more than one temporal option is set |
| Analysis script reads garbage from `csvs/*.csv` | the header is a single row now, not a three-level MultiIndex; see [summary outputs](#summaries) |
| `FileNotFoundError` or `MissingInputException` on a `csvs/individual/…` path | that stage no longer exists; slice the horizon column instead |

## Migrating a fork {#fork-migration}

This section is for maintainers of forks and soft-forks. The preparation stage was
consolidated from six rules into one, so the conflicts concentrate in the
`add_*`/`prepare_*` scripts and in the deleted `solve_*.smk` files.

Work in an order that gets you to a workflow that parses as early as possible:

1. Migrate your config templates and scenario YAMLs first, including every
   `config/test/*.yaml` that you carry. Nothing parses until they have a top-level
   `planning_horizons`, a decision on `sector.enabled` and an integer
   `clustering.cluster_network.n_clusters`. Work through
   [Rewrite your config](#rewrite-your-config) for each of them.
2. Accept the deletion of `rules/solve_*.smk` and of the three summary scripts. Keep
   your patches to them in a scratch diff instead of trying to merge them
   ([Rule files](#fork-rules), [Scripts](#fork-scripts)).
3. Remove `{clusters}`, `{opts}` and `{sector_opts}` from your own rules, rename
   `{planning_horizons}` to `{horizon}`, and update your `resources()`/`RESULTS`
   paths ([Fixing file references](#fork-paths)).
4. Get `snakemake -n --configfile <your test config>` to build a DAG, fixing helper
   renames as they show up ([Configuration access](#fork-config)).
5. Re-apply your input and param additions in `get_compose_inputs` and in the
   `compose_network` rule's `params:` ([Porting code](#fork-porting)).
6. Re-apply your script patches inside the new `main()` functions, respecting the new
   call order ([Scripts](#fork-scripts)).
7. Run the test matrix ([Verifying the merge](#fork-tests)).

### Rule files {#fork-rules}

**Deleted, together with all their rules:**

| File | Rules |
|------|-------|
| `rules/solve_electricity.smk` | `solve_network`, `solve_operations_network` |
| `rules/solve_overnight.smk` | `solve_sector_network` |
| `rules/solve_myopic.smk` | `add_existing_baseyear`, `add_brownfield`, `solve_sector_network_myopic` |
| `rules/solve_perfect.smk` | `add_existing_baseyear`, `prepare_perfect_foresight`, `solve_sector_network_perfect`, `make_summary_perfect` |

The input helpers that were defined in those files went with them:
`input_profile_tech_brownfield`, `input_network_year` and
`input_networks_make_summary_perfect`. Their equivalents are now assembled in
`get_compose_inputs`. The `ruleorder: add_existing_baseyear > add_brownfield`
disambiguation is not needed any more, because both are branches inside one rule.

**Added:**

- `rules/compose.smk` — `compose_network` for every foresight mode and both model
  types, plus the input-assembly function `get_compose_inputs(w)`.
- `rules/solve.smk` — `solve_network` (all modes) and `solve_operations_network`.

**Removed from surviving files:** `add_electricity`, `prepare_network`
(`build_electricity.smk`); `prepare_sector_network` (`build_sector.smk`);
`make_global_summary`, `make_cumulative_costs` (`postprocess.smk`).

**Renamed:** `plot_power_network_clustered` → `plot_clustered_network`;
`build_clustered_solar_rooftop_potentials` → `build_solar_rooftop_potentials`.

**Split:** `build_shapes` → `build_shapes` + `build_offshore_shapes` +
`build_nuts3_shapes`; `build_energy_totals` → `build_energy_totals` +
`build_co2_totals` + `build_transformation_output_coke`.

**New in surviving files:** `chain_busmaps`, `cluster_electricity_demand`.

!!! note "Where foresight branching lives now"
    The foresight-conditional *includes* are gone: the `Snakefile` includes
    `rules/compose.smk` and `rules/solve.smk` unconditionally. Foresight branching
    itself is not gone. It now takes the form of `if config["foresight"] …` blocks
    around rule definitions in `rules/postprocess.smk` and around the output lists in
    the `Snakefile`. Put your own conditional rules there. Note that these blocks read
    the parse-time `config` dict and not `config_provider`, so a scenario that
    overrides `foresight` or `planning_horizons` does not change which rules exist or
    what the collect rules ask for.

### Scripts {#fork-scripts}

**Deleted:** `make_summary_perfect.py`, `make_global_summary.py` and
`make_cumulative_costs.py`, all folded into `make_summary.py`, where cumulative costs
are computed by `calculate_cumulative_costs()`.
`plot_power_network_perfect.py` is deleted too: `plot_power_network.py` covers every
foresight mode.

**Renamed:** `build_clustered_solar_rooftop_potentials.py` →
`build_solar_rooftop_potentials.py` (content unchanged).

**Added:** `compose_network.py` (the driver), `co2_budget.py`, `chain_busmaps.py`,
`cluster_electricity_demand.py`, `build_co2_totals.py`,
`build_transformation_output_coke.py`, `build_nuts3_shapes.py`,
`build_offshore_shapes.py`.

**Turned into libraries.** `add_electricity.py`, `add_existing_baseyear.py`,
`add_brownfield.py`, `prepare_network.py`, `prepare_sector_network.py` and
`prepare_perfect_foresight.py` still exist, but Snakemake does not execute them any
more. Their `if __name__ == "__main__"` blocks became `main()` functions that
`compose_network.py` imports and calls. The signatures differ, so check the one you
are patching:

```python
# add_electricity.py, add_existing_baseyear.py
def main(n: pypsa.Network, inputs, params, costs: pd.DataFrame) -> None

# prepare_network.py
def main(n: pypsa.Network, inputs, params, costs: pd.DataFrame, nyears: float) -> None

# prepare_sector_network.py
def main(n: pypsa.Network, inputs, params, costs: pd.DataFrame,
         nyears: float, current_horizon: int) -> None

# add_brownfield.py
def main(n: pypsa.Network, n_previous: pypsa.Network, inputs, params,
         current_horizon: int, renewable_carriers: list[str]) -> None

# prepare_perfect_foresight.py — returns the merged multi-period network
def main(n: pypsa.Network, n_previous: pypsa.Network | None, params,
         current_horizon: int) -> pypsa.Network
```

`compose_network.py` imports them under descriptive aliases and calls them in this
order:

```python
add_electricity_components(n, inputs, params, costs)
if sector.enabled:
    add_sector_components(n, inputs, params, costs, nyears, horizon)
apply_temporal_aggregation(n, inputs, params)
if foresight != "overnight" and first horizon:
    add_existing_capacities(n, inputs, params, costs)
    if foresight == "myopic":
        adjust_renewable_capacity_limits(n, str(horizon), renewable_carriers)
if foresight == "myopic" and not first horizon:
    apply_brownfield(n, n_previous, inputs, params, horizon, renewable_carriers)
prepare_network_for_solving(n, inputs, params, costs, nyears)
if foresight == "perfect":
    n = prepare_perfect_foresight(n, n_previous, params, horizon)
    apply_phase_outs(n, params.existing_capacities["phase_outs"], horizons)
apply_co2_budget_constraints(n, inputs, params, nyears, horizon)
maybe_adjust_costs_and_potentials(n, adjustments["electricity"], horizon)
maybe_adjust_costs_and_potentials(n, adjustments["sector"], horizon)
```

!!! warning "The order changed, not just the packaging"
    The chain used to be `add_electricity` → `prepare_network` →
    `prepare_sector_network` → `add_existing_baseyear`/`add_brownfield`. Now the body
    of `prepare_network` runs *after* the sector layer and *after* existing and
    brownfield capacities, temporal aggregation runs on the composed network, and the
    CO₂ limits and both `adjustments` blocks are applied last. So a patch in
    `prepare_network.py` that assumed an electricity-only network will now see a full
    sector network with existing vintages.

### Porting code inside a `main()` {#fork-porting}

Inside these functions the `snakemake` object is not available. Use the `inputs` and
`params` arguments instead:

| Old pattern | New pattern |
|----------------|-------------|
| `snakemake.input["powerplants"]` | `inputs["powerplants"]` |
| `snakemake.input.fuel_price` | `inputs.fuel_price` |
| `snakemake.input.network` | the `n` argument — `compose_network.py` opens the clustered network and passes it in |
| the previous horizon's network | `inputs.network_previous` |
| `snakemake.params.costs` | `params.costs` |
| `snakemake.config["sector"]["enabled"]` | `params.sector["enabled"]` |

Config values reach the scripts only through `params`, so anything you used to read
via `snakemake.config[...]` needs an entry in the `compose_network` rule:

```python
rule compose_network:
    params:
        your_custom_param=config_provider("path", "to", "config"),
```

To add a custom input file, extend `get_compose_inputs(w)` in `rules/compose.smk`. It
returns the whole input dictionary for the rule. Since `compose_network` declares
`input: unpack(get_compose_inputs)` and nothing else, there is no second place to
look.

### Fixing file references {#fork-paths}

Remove the scenario wildcards `{clusters}`, `{opts}` and `{sector_opts}` from your
paths, and rename `{planning_horizons}` to `{horizon}`. All other wildcards are
unchanged: `{run}` still exists when `run.name` is set and scenario management is
enabled, and `{technology}`, `{carrier}`, `{year}`, `{country}` and `{cutout}` are
untouched.

| Old | New |
|--------|-----|
| `profile_{clusters}_{technology}.nc` | `profile_{technology}.nc` |
| `powerplants_s_{clusters}.csv` | `powerplants.csv` |
| `costs_{planning_horizons}_processed.csv` | `costs_{horizon}_processed.csv` |
| `regions_onshore_base_s_{clusters}.geojson` | `onshore_regions.geojson` |
| `regions_offshore_base_s_{clusters}.geojson` | `offshore_regions.geojson` |

See [Files and directories](#files) for the full mapping, and
[wildcards](wildcards.md) for the wildcards that remain.

!!! warning "A leftover wildcard does not raise"
    `wildcard_constraints` no longer declares `clusters`, `opts`, `sector_opts` or
    `planning_horizons`. A wildcard that you forget to remove therefore falls back to
    Snakemake's default `.+`, which even matches across directory separators. The
    symptom is an `AmbiguousRuleException` or a `MissingInputException` that names a
    nonsensical path. So before anything else, run

    ```console
    $ grep -rn '{clusters}\|{opts}\|{sector_opts}\|{planning_horizons}' rules/ Snakefile
    ```

    and expect no hits outside docstrings and schema field descriptions.

### Configuration access {#fork-config}

Rules get their parameters through `config_provider(...)`, and input functions through
`get_config(w)`. Both resolve the fully merged, scenario-aware config and are cached
per wildcard set (`rules/common.smk`).

!!! warning "`get_config` changed meaning"
    The nested-lookup helper `get_config(config, keys, default)` was renamed to
    `navigate_config(config, keys, default)`. The name `get_config` now takes a
    wildcards object and returns the whole merged config: `get_config(w)`. A fork that
    calls the old signature merges without a textual conflict and then fails at DAG
    build with an opaque error, so grep your fork for `get_config(`. And because the
    result is cached and shared between callers, never mutate it in your own input
    function.

The wildcard-to-config translation layer `_helpers.update_config_from_wildcards` was
deleted along with the option wildcards. The generic `_helpers.parse` helper survives,
but it is no longer used for option strings, so a fork that added its own option
tokens has to expose them as config keys instead.

### Verifying the merge {#fork-tests}

```console
$ pixi run unit-tests
$ pixi run integration-tests
$ pixi run integration-tests-streamlined
```

Add your fork's cases to the new `test/test_compose_network.py`,
`test/test_co2_budget.py`, `test/test_make_summary.py` and
`test/test_prepare_perfect_foresight.py` instead of to the deleted test files. If you
added config keys, run `pixi run generate-config` to refresh
`config/config.default.yaml` and `config/schema.default.json`, otherwise the schema
test fails. The cheapest regression signal is
`solving.check_objective.expected_value` in the test configs — expect to re-baseline
it for your fork.
