<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Migration Guide {#migration}

This guide covers the streamlined workflow refactor
([#1838](https://github.com/PyPSA/pypsa-eur/pull/1838)). It restructures the
workflow around the stage progression `base → simplified → clustered → composed →
solved` and moves all scenario information out of filenames and into the
configuration file. Other changes that landed in the same release are listed in
the [release notes](release_notes.md); this guide is about the refactor.

All three foresight modes now run through the same stages and differ only in how
they loop: overnight composes and solves one horizon, myopic composes and solves
each horizon in turn with `compose_network` reading the previous horizon's *solved*
network, and perfect foresight composes every horizon in sequence — each reading the
previous horizon's *composed* network — and then solves them together. See
[foresight](foresight.md) for details.

Everyone upgrading should read [Before you start](#before-you-start) and
[Rewrite your config](#rewrite-your-config). If you maintain a fork or soft-fork
of PyPSA-Eur, also read [Migrating a fork](#fork-migration) — the script and rule
structure changed substantially and you will get merge conflicts.

!!! danger "Nothing tells you that your old config stopped working"
    An old config still validates without a single error. Keys that were renamed
    or restructured are simply ignored, and the setting they used to control falls
    back to its default — so the run succeeds and models something you did not ask
    for. A leftover `scenario:` block leaves you at `planning_horizons: [2050]` and
    the default `n_clusters: 50`; a leftover flat `co2_budget:` mapping leaves you
    on the default CO₂ trajectory; a leftover `resolution_elec:` leaves you at
    native hourly resolution; and a config without a `sector:` block builds a
    fully sector-coupled model. Work through
    [Rewrite your config](#rewrite-your-config) key by key rather than waiting for
    an error message, and check [Symptoms and causes](#symptoms) if results look
    off.

## Before you start {#before-you-start}

**Your retrieved data and cutouts stay valid.** Nothing about the retrieved input
data or the weather cutouts changed, so you do not need to re-download or re-cut
anything. Two exceptions, both user-supplied files keyed by the old `{clusters}`
wildcard: if you use `clustering.mode: custom_busmap` or `custom_busshapes`,
rename your own files to `data/busmaps/simplified_{n_clusters}_{base_network}.csv`
and `data/busshapes/simplified_{n_clusters}_{base_network}.geojson`.

**Your `resources/` and `results/` directories are stale.** Nearly every
intermediate and final filename changed, so Snakemake will not recognise the
existing files and will rebuild the chain from `base.nc` onwards. Old files are not
deleted; they are simply orphaned. Delete them once you are satisfied with the new
run, or keep them for comparison — but budget for a full rebuild. If you want the
new summaries without re-solving, see
[reusing solved networks](#reusing-solved-networks). Before comparing anything to an
old run, read [Behaviour changes](#behaviour-changes): several changes are not
renames and cannot be undone from the config.

**Your analysis scripts need attention.** The summary CSVs changed shape and the
per-horizon `csvs/individual/*.csv` files no longer exist. See
[Summary outputs](#summaries).

The migration itself works through the following, in order:

1. [Declare your model type](#model-type) — electricity-only or sector-coupled.
2. [The `scenario` block](#scenario-block) — horizons and cluster count.
3. [`opts` and `sector_opts`](#opts-translation) — every option token becomes a config key.
4. [CO₂ constraints](#co2-budget) — restructured, and the units changed.
5. [Temporal resolution](#temporal) — one setting instead of two.
6. [Other renamed keys](#other-keys) — and the new keys you have to decide on.
7. [File paths](#files) — in your own scripts.
8. [Snakemake targets](#targets) — in your commands.

Steps 1 to 6 are all config edits and end with a
[before-and-after example](#example). Then read
[Behaviour changes](#behaviour-changes) and [Summary outputs](#summaries) before you
compare results or run your analysis scripts. If you would rather work from a full
config, `config/test/config.electricity.yaml` is the cleanest end-to-end
electricity-only one in the repository, and `config/test/config.myopic.yaml` shows
the multi-horizon sector-coupled form.

## Rewrite your config {#rewrite-your-config}

### Declare your model type {#model-type}

Whether you build an electricity-only or a sector-coupled model used to be decided
by the target you typed — `solve_elec_networks` versus `solve_sector_networks`.
It is now a configuration key, `sector.enabled`, and both model types share one
set of rules and one default target.

!!! danger "Electricity-only models must opt out explicitly"
    `sector.enabled` defaults to `true`. An old electricity-only config has no
    `sector:` block at all, so after migration it validates and then composes and
    solves a full sector-coupled model, with the runtime and data retrieval that
    implies. Set it explicitly:

    ```yaml
    sector:
      enabled: false
    ```

### The `scenario` block is gone {#scenario-block}

The `scenario:` section defined the values of the `{clusters}`, `{opts}`,
`{sector_opts}` and `{planning_horizons}` wildcards. All four wildcards were
removed, so the block is no longer read.

| Legacy | New |
|--------|-----|
| `scenario.planning_horizons` | `planning_horizons` (top level) |
| `scenario.clusters` | `clustering.cluster_network.n_clusters` |
| `scenario.opts` | explicit config keys, see [below](#opts-translation) |
| `scenario.sector_opts` | explicit config keys, see [below](#opts-translation) |

`planning_horizons` is a top-level list of years and directly controls the
`{horizon}` wildcard. A bare scalar (`planning_horizons: 2030`) is accepted and
normalised to a list.

`n_clusters` takes a single integer, not a list — the wildcard sweep is gone. To
compare several cluster counts, use `run.scenarios` with one scenario per count
(see [configuration](configuration.md)). For administrative clustering set
`clustering.mode: administrative` instead of the former `clusters: [adm]`.

!!! warning "Overnight runs must declare exactly one horizon"
    With `foresight: overnight`, `planning_horizons` must contain exactly one
    year. Otherwise the workflow fails while building the DAG, before any job
    starts, with `Overnight optimization can only be run for a single planning
    horizon.` An overnight setup could previously list several years and solve each
    independently; that sweep now belongs in `run.scenarios`.

!!! warning "What a scenario may not change"
    Collection targets and the default target are built from the *base* config, so
    `foresight` and `planning_horizons` must be the same for every scenario. Set
    them at the top level and run different foresight modes or horizon sets as
    separate workflows with distinct `run.name`. Everything the rules read through
    `config_provider` — including `clustering.cluster_network.n_clusters`,
    `co2_budget`, `sector` and `solving` — may vary per scenario.

### Translating `opts` and `sector_opts` {#opts-translation}

Every option token has a configuration equivalent. The translation layer
(`_helpers.update_config_from_wildcards`) was deleted, so leftover option strings
have no effect whatsoever.

`opts` tokens:

| Legacy token | New config |
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

| Legacy token | New config |
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
| `cb<x>ex`, `cb<x>be` | no replacement: `build_carbon_budget` is no longer called, so specify the resulting per-year values under `co2_budget.upper` |

Three of these tokens carried an implicit unit conversion that you now have to do
yourself: `linemaxext<n>` was given in GW and multiplied by `1e3`, so
`linemaxext20` becomes `20000` (MW); `CH4L<x>` was multiplied by `1e6`, so
`CH4L200` becomes `gaslimit: 200000000`; and `sdr<x>` was a percentage divided by
`100`, so `sdr2` becomes `social_discountrate: 0.02`.

!!! note
    `sector.district_heating.progress` must stay a year-to-value mapping — a bare
    `0` is rejected by the schema. Because config files deep-merge, set every
    horizon you optimise, for example `progress: {2030: 0, 2040: 0, 2050: 0}`.

### CO₂ constraints {#co2-budget}

CO₂ handling is now expressed entirely through `co2_budget`. Three separate legacy
mechanisms feed into it:

| Legacy | New |
|--------|-----|
| `co2_budget: {<year>: <fraction>}` (flat mapping) | `co2_budget.upper: {<year>: <fraction>}` |
| `electricity.co2limit_enable` / `co2limit` / `co2base` | `co2_budget.upper` with `relative: false` |
| `energy.emissions` | `co2_budget.emissions_scope` |

```yaml
co2_budget:
  emissions_scope: CO2     # was energy.emissions
  relative: true           # upper/lower are fractions of 1990 emissions
  upper:
    2030: 0.45
    2050: 0.0
  lower:                   # optional floor, see the warning below
```

Three traps are worth calling out.

**Units changed for absolute limits.** With `relative: false`, `upper` and `lower`
are in **Gt CO₂ per year**, whereas the old `electricity.co2limit` was in tonnes.
Divide by `1e9`: `co2limit: 77500000.0` becomes `upper: 0.0775`. Migrating the
number verbatim gives a limit a billion times too loose.

**Electricity-only runs are now capped by default.** An electricity-only model
previously had no CO₂ constraint unless you set
`electricity.co2limit_enable: true`. Now the same `co2_budget` block applies to
electricity-only runs, and the shipped default ends at `2050: 0.0`. An
electricity-only run at `planning_horizons: [2050]` is therefore capped at zero
emissions and may be infeasible. To reproduce the unconstrained default, set
`co2_budget: {upper: null}`.

**Partial overrides merge, they do not replace.** Config files are deep-merged, so
writing only `co2_budget: {upper: {2030: 0.4}}` leaves the other six default years
in place. Set the years you do not want to `null`.

A mapping and a scalar mean different things: a mapping constrains each listed year
individually, while a scalar is a single budget across all periods — under perfect
foresight it becomes one cumulative constraint attached to the final horizon rather
than a per-period cap. See [foresight](foresight.md) for the exact semantics per
mode.

The optional `lower` bound is reliable only for electricity-only runs, where it is
added as a `primary_energy` constraint whose `>=` sense PyPSA honours. In
sector-coupled runs the constraint type is `co2_atmosphere` (or `Co2Budget` for a
scalar budget under perfect foresight), and those solver routines build
`lhs <= rhs` regardless of the stored sense, so a `lower` value acts as a second,
tighter *upper* limit. Leave `co2_budget.lower` unset for sector-coupled runs.

!!! warning "List every horizon you optimise"
    A horizon missing from `co2_budget.upper` gets **no CO₂ constraint at all** —
    the workflow logs an info line and moves on. Earlier releases interpolated
    between the neighbouring years instead. Make sure every entry of
    `planning_horizons` appears in `upper` (and `lower`, if used).

### Temporal resolution {#temporal}

`clustering.temporal.resolution_elec` and `resolution_sector` are replaced by
three mutually exclusive integer options under `clustering.temporal`. The same
setting now applies to electricity-only and sector-coupled runs, so if you
previously set both to different values you have to choose one. Setting more than
one option raises `clustering.temporal: only one of averaging, segmentation and
representative may be set`.

| Legacy string | New |
|---------------|-----|
| `resolution_elec: 24h` / `resolution_sector: 24H` | `averaging: 24` |
| `resolution_sector: 8760h` | `averaging: 8760` |
| `resolution_*: 4380seg` | `segmentation: 4380` |
| `resolution_*: 3sn` | `representative: 3` |
| `resolution_*: false` | leave all three `false` |

`averaging` averages over `n` hours, `segmentation` aggregates the year into `n`
segments using `tsam`, and `representative` keeps every `n`-th snapshot.

!!! note
    `averaging` and `representative` reproduce the previous numbers exactly.
    `segmentation` does not: the segments are now derived from
    `networks/clustered.nc` rather than from the fully prepared network, so
    segment boundaries and weights differ.

### Other renamed keys {#other-keys}

| Legacy | New |
|--------|-----|
| `lines.max_extension` | `lines.s_nom_max_extension` |
| `links.max_extension` | `links.p_nom_max_extension` (default 30000 MW, now applies to sector runs too) |
| `solving.options.rolling_horizon` (+ `horizon`, `overlap`) | `solving.operations.rolling_horizon` (+ `horizon`, `overlap`) |

The rolling-horizon settings moved into their own block because they configure
`solve_operations_network`. The key `solving.options.rolling_horizon` still
exists, but now affects only the capacity-expansion `solve_network` rule, so
capacity expansion can be followed by rolling-horizon dispatch.

Two more keys are new and worth a decision:

`existing_capacities.phase_outs` holds the national policy phase-outs that used to
be hardcoded, as a list of `{carriers, countries, year}` rules. They cap
conventional asset lifetimes under **perfect foresight only** — generators for
electricity-only runs, links for sector-coupled ones — and the defaults reproduce
the previous behaviour. Note that lists are replaced wholesale rather than merged,
so a list you supply replaces all the defaults, including the ones you did not mean
to drop.

`run.default_target_rule` selects which collect rule `snakemake` builds when you
name no target. It defaults to `all`.

Finally, `solving.check_objective.expected_value` now accepts either a single
value or a `{<horizon>: value}` mapping, so each `solved_{horizon}.nc` of a
multi-horizon run is checked against its own benchmark. Its default tolerances
were tightened (`atol` 1e6 → 1e4, `rtol` 0.01 → 0.001).

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

These are not renames. With an equivalent configuration the model itself changes,
so do not expect to reproduce old numbers without acting on them.

**Overnight runs use per-horizon costs.** Every layer is costed with
`costs_{horizon}_processed.csv`, and `costs.year` no longer influences
composition. Previously an overnight run used `costs.year` (default 2050)
regardless of the horizon, and in myopic and perfect runs the electricity layer
used `costs.year` while the sector layer used the horizon. All costs, efficiencies
and lifetimes change whenever those years differed.
**Reproduce the old behaviour:** set `planning_horizons: [<old costs.year>]`.

**Electricity-only runs are CO₂-capped by default**, and absolute caps are in Gt
rather than tonnes, so a default electricity-only run flips from unconstrained to
zero-emission.
**Reproduce the old behaviour:** `co2_budget: {upper: null, lower: null}`, or
divide your old tonne value by `1e9`.

**Relative budgets are always measured against 1990.**
`energy.base_emissions_year` no longer influences the CO₂ budget baseline (it still
controls `build_co2_totals`). For electricity-only runs the baseline still spans
every sector flagged in `sector:`, which makes a relative budget nearly
non-binding.
**Reproduce the old behaviour:** nothing for sector runs with a 1990 baseline;
otherwise rescale `upper` manually or use `relative: false`. For electricity-only,
switch the `sector.transport`, `heating`, `industry` and `agriculture` flags off.

**Missing horizons drop the CO₂ constraint** instead of interpolating, which
leaves that horizon — and in myopic mode every later one — unconstrained.
**Reproduce the old behaviour:** list every horizon explicitly in
`co2_budget.upper`.

**Sector-coupled HVDC extension headroom rises from 20 to 30 GW**, and a finite
`lines.s_nom_max` or `links.p_nom_max` now binds in sector runs where it was
previously ignored.
**Reproduce the old behaviour:** `links.p_nom_max_extension: 20000`, and
`lines.s_nom_max: .inf` / `links.p_nom_max: .inf`.

**Both `adjustments` blocks are applied**, always with the current horizon, and
after transmission costs are set. Sector runs now also pick up
`adjustments.electricity`, electricity-only runs also pick up `adjustments.sector`,
and Line/Link `capital_cost` adjustments now take effect.
**Reproduce the old behaviour:** keep entries only in the block your model type
used before, and drop Line/Link `capital_cost` adjustments.

**`segmentation` is computed from the clustered network**, so segment boundaries
and weights differ and with them capacity factors and the objective.
**Reproduce the old behaviour:** use `averaging` or `representative`, which are
numerically unchanged.

**Electricity-only myopic and perfect foresight now run.** An electricity-only
config previously ignored `foresight` entirely and always solved a single period;
now `foresight: myopic` or `perfect` produces a completely different model, with an
existing fleet and capacity accumulating across horizons.
**Reproduce the old behaviour:** use `foresight: overnight` with one horizon.

### If you use perfect foresight {#perfect}

Sector-coupled perfect foresight works again: earlier releases aborted with `PyPSA
versions >=1.0 are not supported for perfect foresight`. Because the assembly was
rewritten rather than translated line by line, treat perfect-foresight results as
new rather than comparable. Three specifics:

- **AC lines stay lines.** Previously every AC line was replaced by a DC link with
  a flat 3%/1000 km efficiency and one vintage per investment period. Now the
  lines survive as a single extendable set shared across all periods, with PyPSA's
  resistance-based loss model (`solving.options.transmission_losses`) and the
  linearised power flow. Transmission is one decision for the whole horizon
  instead of a staged build-out, so objective, line capacities and cross-border
  flows change. This cannot be recovered from the config; the conversion code no
  longer exists.
- **Per-horizon maps are not produced.** Only the summary CSVs and the
  `graphs/*.svg` summary plots are generated in this mode.
- **The summary CSVs changed**, both in content and in which files exist. See
  [Summary outputs](#summaries).

## Summary outputs and your analysis scripts {#summaries}

This is the change most likely to break your own post-processing.

The two-stage summary is gone. `make_summary` no longer writes per-horizon
`csvs/individual/<metric>_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv`
files, and `make_global_summary` and `make_cumulative_costs` were removed as
separate rules. A single `make_summary` writes `results/csvs/<metric>.csv`
directly, built from all solved networks at once via `pypsa.NetworkCollection`.

Two consequences:

- **`csvs/individual/` no longer exists.** The planning horizons are now
  **columns** of the aggregated CSV. To read a single horizon, slice the column:
  `pd.read_csv("results/csvs/costs.csv", index_col=[0, 1, 2])["2030"]`. The column
  label is the horizon year, read back from CSV as a string. Index depth varies by
  metric (three levels for `costs`, four for `nodal_costs`, one for `metrics`).
- **The aggregated CSVs kept their paths but changed their header.** They used to
  carry a three-level `(cluster, opt, planning_horizon)` column MultiIndex.
  Readers that still pass `header=[0, 1, 2]` will silently mis-parse the new
  single-row header rather than fail.

Under perfect foresight the set of files also changes. Summaries are recomputed from
PyPSA `statistics` rather than the deleted `make_summary_perfect.py`, so reported
costs are undiscounted and not split by period activity. `csvs/supply.csv` and
`csvs/supply_energy.csv` are replaced by `csvs/energy_balance.csv` and
`csvs/nodal_energy_balance.csv`, which are not column-for-column equivalent, and
`csvs/price_statistics.csv` and `csvs/co2_emissions.csv` have no replacement — CO₂
figures can be recovered from the `co2` rows of `energy_balance.csv`. If you run
perfect foresight, also read
[If you use perfect foresight](#perfect).

Two per-solve outputs are no longer written: the config dumps
`results/configs/config.*.yaml`, and `results/models/*.nc` for
`solving.options.store_model`. The effective configuration is still embedded in
each solved network as `n.meta`.

### Reusing solved networks from an old run {#reusing-solved-networks}

Because the filenames changed, Snakemake sees no solved network and would re-solve
everything. To get the new summaries without re-solving, rename your old solved
networks to `results/<run>/networks/solved_{horizon}.nc` and request the summary:

```console
$ snakemake results/<run>/csvs/costs.csv --configfile <configfile>
```

The rename is all it takes. Snakemake does not rebuild the upstream chain of an
output that already exists, so the missing `composed_*.nc` and `resources/` files
are not re-created. `make_summary` reads every horizon of `planning_horizons` in
one job — only the last one under perfect foresight — so all of them must be
present.

## Files and directories {#files}

Paths below omit the run directory. It is inserted as `resources/<run>/…`,
`results/<run>/…`, `logs/<run>/…` and `benchmarks/<run>/…` when `run.name` is set,
as a literal directory, or as a `{run}` wildcard when `run.scenarios.enable` is
also true. A non-empty `run.prefix` adds a further level in front. With the
shipped defaults — empty `run.prefix` and `run.name` — the paths are exactly as
shown.

### Network stages

| Legacy | New |
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

| Legacy | New |
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

Note the word order flip on the region files: `regions_onshore_*` became
`onshore_regions*`. Beyond the table, the general rule for sector resources is that
`_s_{clusters}` is dropped and `{planning_horizons}` becomes `{horizon}`, for
example `biomass_potentials_s_{clusters}_{planning_horizons}.csv` →
`biomass_potentials_{horizon}.csv` and
`cop_profiles_base_s_{clusters}_{planning_horizons}.nc` →
`cop_profiles_{horizon}.nc`.

### Maps and graphs

| Legacy | New |
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

`plot_summary` now writes SVG rather than PDF.

## Snakemake targets {#targets}

Every foresight mode and both model types support the default target, so you can
always run the whole workflow with:

```console
$ snakemake --configfile <configfile>
```

The collection rules were unified:

| Legacy rule | New rule |
|-------------|----------|
| `cluster_networks` | `cluster_networks` (now `networks/clustered.nc` and `busmap.csv`) |
| `prepare_elec_networks`, `prepare_sector_networks` | `compose_networks` |
| `solve_elec_networks`, `solve_sector_networks`, `solve_sector_networks_perfect` | `solve_networks` |
| `plot_power_networks_clustered` | `plot_power_networks` |
| — | `solve_operations_networks` (new) |

`solve_networks` requests only the last horizon's solved network; the earlier
horizons follow from the dependency chain. `compose_networks` expands over the full
`planning_horizons` list, and so does `process_costs` except under
`foresight: overnight`, where it still derives a single file from `costs.year` —
composition always reads `costs_{horizon}_processed.csv`, so on an overnight run
this collect target can build a cost file the run itself never uses.

The `plot_power_networks` collect rule currently asks for a filename the producing
rule does not write, so it fails at DAG build. Until that is fixed, request
`resources/maps/clustered_network.pdf` or the `plot_clustered_network` rule
directly.

When targeting files directly:

| Legacy | New |
|--------|-----|
| `snakemake resources/networks/base_s_37.nc` | `snakemake resources/networks/clustered.nc` |
| `snakemake resources/networks/base_s_37_elec_Co2L-24h.nc` | `snakemake resources/networks/composed_2030.nc` |
| `snakemake results/networks/base_s_37_elec_Co2L-24h.nc` | `snakemake results/networks/solved_2030.nc` |
| `snakemake results/networks/base_s_37_Co2L-24h_Co2L0-1H-T-H-B-I_2030.nc` | `snakemake results/networks/solved_2030.nc` |

The `**config["scenario"]` expansion pattern is gone; collection rules take the
`{horizon}` wildcard from `config["planning_horizons"]`.

## Symptoms and causes {#symptoms}

| Symptom | Likely cause |
|---------|--------------|
| Run builds hydrogen and industry inputs, or takes far longer than the old electricity-only run | `sector.enabled` defaults to `true`; set it to `false`. Note the default target still builds heat-pump COP profiles even for electricity-only runs — use `run.default_target_rule: solve_networks` to avoid them |
| Runs only 2050, or only one horizon | `planning_horizons` still inside `scenario:`; move it to the top level |
| Runs at 50 clusters regardless of config | `scenario.clusters` not moved to `clustering.cluster_network.n_clusters` |
| Model solves at hourly resolution and takes forever | `resolution_elec`/`resolution_sector` not translated to `averaging`/`segmentation`/`representative` |
| Electricity-only run infeasible right after upgrading | default `co2_budget` now applies; set `co2_budget: {upper: null}` or a real cap |
| CO₂ cap appears to do nothing | absolute value still in tonnes — divide by `1e9`; or the horizon is missing from `co2_budget.upper` |
| Grid expansion larger than expected | `max_extension` not renamed, so the default 20000/30000 MW applies |
| `MissingInputException` on a `data/busmaps/…` or `data/busshapes/…` path | rename your custom busmap or bus shapes file, see [before you start](#before-you-start) |
| `ValueError: Overnight optimization can only be run for a single planning horizon.` | `foresight: overnight` with several horizons; use one horizon or `run.scenarios` |
| `ValidationError` on `clustering.cluster_network.n_clusters` | `n_clusters` must be an integer; `all` and `null` are rejected by the schema |
| `ValidationError` on `sector.district_heating.progress` | it must stay a year-to-value mapping, not a scalar |
| `ValueError: clustering.temporal: only one of averaging, segmentation and representative may be set, got …` | more than one temporal option set |
| Analysis script reads garbage from `csvs/*.csv` | the header is one row now, not a three-level MultiIndex; see [summary outputs](#summaries) |
| `FileNotFoundError` or `MissingInputException` on a `csvs/individual/…` path | that stage no longer exists; slice the horizon column instead |

## Migrating a fork {#fork-migration}

This section is for maintainers of forks and soft-forks. The preparation stage was
consolidated from six rules into one, so conflicts concentrate in the
`add_*`/`prepare_*` scripts and the deleted `solve_*.smk` files.

Work in an order that reaches a parseable workflow as early as possible:

1. Migrate your config templates and scenario YAMLs first, including every
   `config/test/*.yaml` you carry — nothing parses until they have top-level
   `planning_horizons`, a `sector.enabled` decision and an integer
   `clustering.cluster_network.n_clusters`. Work through
   [Rewrite your config](#rewrite-your-config) for each of them.
2. Accept the deletion of `rules/solve_*.smk` and the three summary scripts, and
   keep your patches to them in a scratch diff rather than trying to merge them
   ([Rule files](#fork-rules), [Scripts](#fork-scripts)).
3. Strip `{clusters}`, `{opts}` and `{sector_opts}` from your own rules, rename
   `{planning_horizons}` to `{horizon}`, and update your `resources()`/`RESULTS`
   paths ([Fixing file references](#fork-paths)).
4. Get `snakemake -n --configfile <your test config>` to build a DAG, fixing helper
   renames as they surface ([Configuration access](#fork-config)).
5. Re-apply your input and param additions in `get_compose_inputs` and the
   `compose_network` rule's `params:` ([Porting code](#fork-porting)).
6. Re-apply your script patches inside the `main()` functions, respecting the new
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

The input helpers defined in those files went with them:
`input_profile_tech_brownfield`, `input_network_year` and
`input_networks_make_summary_perfect`. Their equivalents are now assembled in
`get_compose_inputs`. The `ruleorder: add_existing_baseyear > add_brownfield`
disambiguation is no longer needed, because both are branches inside one rule.

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
    The foresight-conditional *includes* are gone — the `Snakefile` includes
    `rules/compose.smk` and `rules/solve.smk` unconditionally. Foresight branching
    itself is not gone: it now takes the form of `if config["foresight"] …` blocks
    wrapping rule definitions in `rules/postprocess.smk` and the output lists in
    the `Snakefile`. Put your own conditional rules there. Note that these blocks
    read the parse-time `config` dict rather than `config_provider`, so a scenario
    that overrides `foresight` or `planning_horizons` does not change which rules
    exist or what the collect rules request.

### Scripts {#fork-scripts}

**Deleted:** `make_summary_perfect.py`, `make_global_summary.py` and
`make_cumulative_costs.py` — all folded into `make_summary.py`, where cumulative
costs are computed by `calculate_cumulative_costs()`.

**Renamed:** `build_clustered_solar_rooftop_potentials.py` →
`build_solar_rooftop_potentials.py` (content unchanged).

**Added:** `compose_network.py` (the driver), `co2_budget.py`, `chain_busmaps.py`,
`cluster_electricity_demand.py`, `build_co2_totals.py`,
`build_transformation_output_coke.py`, `build_nuts3_shapes.py`,
`build_offshore_shapes.py`.

**Refactored into libraries.** `add_electricity.py`, `add_existing_baseyear.py`,
`add_brownfield.py`, `prepare_network.py`, `prepare_sector_network.py` and
`prepare_perfect_foresight.py` are retained but no longer executed by Snakemake.
Their `if __name__ == "__main__"` blocks became `main()` functions that
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
    Previously the chain was `add_electricity` → `prepare_network` →
    `prepare_sector_network` → `add_existing_baseyear`/`add_brownfield`. Now
    `prepare_network`'s body runs *after* the sector layer and *after* existing and
    brownfield capacities, temporal aggregation runs on the composed network, and
    the CO₂ limits and both `adjustments` blocks are applied last. A patch in
    `prepare_network.py` that assumed an electricity-only network will now see a
    full sector network with existing vintages.

### Porting code inside a `main()` {#fork-porting}

Inside these functions the `snakemake` object is not available; use the `inputs`
and `params` arguments instead.

| Legacy pattern | New pattern |
|----------------|-------------|
| `snakemake.input["powerplants"]` | `inputs["powerplants"]` |
| `snakemake.input.fuel_price` | `inputs.fuel_price` |
| `snakemake.input.network` | the `n` argument — `compose_network.py` opens the clustered network and passes it in |
| the previous horizon's network | `inputs.network_previous` |
| `snakemake.params.costs` | `params.costs` |
| `snakemake.config["sector"]["enabled"]` | `params.sector["enabled"]` |

Config values reach the scripts only through `params`, so anything you read via
`snakemake.config[...]` needs a corresponding entry in the `compose_network` rule:

```python
rule compose_network:
    params:
        your_custom_param=config_provider("path", "to", "config"),
```

To add a custom input file, extend `get_compose_inputs(w)` in `rules/compose.smk`;
it returns the whole input dictionary for the rule. Note that `compose_network`
declares `input: unpack(get_compose_inputs)` only, so there is no second place to
look.

### Fixing file references {#fork-paths}

Remove the scenario wildcards `{clusters}`, `{opts}` and `{sector_opts}` from your
paths, and rename `{planning_horizons}` to `{horizon}`. All other wildcards are
unchanged — `{run}` still exists when `run.name` is set and scenario management is
enabled, and `{technology}`, `{carrier}`, `{year}`, `{country}` and `{cutout}` are
untouched.

| Legacy | New |
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
    `planning_horizons`, so a wildcard you forget to remove silently falls back to
    Snakemake's default `.+`, which matches across directory separators. The
    symptoms are an `AmbiguousRuleException` or a `MissingInputException` naming a
    nonsensical path. Before anything else, run

    ```console
    $ grep -rn '{clusters}\|{opts}\|{sector_opts}\|{planning_horizons}' rules/ Snakefile
    ```

    and expect no hits outside docstrings and schema field descriptions.

### Configuration access {#fork-config}

Rules obtain parameters through `config_provider(...)`, and input functions through
`get_config(w)`; both resolve the fully merged, scenario-aware config and are
cached per wildcard set (`rules/common.smk`).

!!! warning "`get_config` changed meaning"
    The nested-lookup helper `get_config(config, keys, default)` was renamed to
    `navigate_config(config, keys, default)`. The name `get_config` now takes a
    wildcards object and returns the whole merged config: `get_config(w)`. A fork
    calling the old signature merges without a textual conflict and then fails at
    DAG build with an opaque error. Grep your fork for `get_config(`. Because the
    result is cached and shared between callers, never mutate it in your own input
    function.

The wildcard-to-config translation layer `_helpers.update_config_from_wildcards`
was deleted along with the option wildcards. The generic `_helpers.parse` helper
survives but is no longer used for option strings, so a fork that added its own
option tokens has to expose them as config keys instead.

### Verifying the merge {#fork-tests}

```console
$ pixi run unit-tests
$ pixi run integration-tests
$ pixi run integration-tests-streamlined
```

Add your fork's cases to the new `test/test_compose_network.py`,
`test/test_co2_budget.py`, `test/test_make_summary.py` and
`test/test_prepare_perfect_foresight.py` rather than to the deleted test files. If
you added configuration keys, run `pixi run generate-config` to refresh
`config/config.default.yaml` and `config/schema.default.json`, otherwise the schema
test fails. `solving.check_objective.expected_value` in the test configs is the
cheapest regression signal — expect to re-baseline it for your fork.
