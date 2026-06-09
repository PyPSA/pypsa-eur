<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->


<a id="config"></a>

# 
# Configuration

PyPSA-Eur has several configuration options which are documented in this section.

<a id="defaultconfig"></a>


## Configuration Files

Any PyPSA-Eur configuration can be set in a `.yaml` file. The default configurations
`config/config.default.yaml` and `config/plotting.default.yaml` are maintained in
the repository and cover all the options that are used/ can be set.

To pass your own configuration, you can create a new file, e.g. `my_config.yaml`,
and specify the options you want to change. They will override the default settings and
options which are not set, will be inherited from the defaults above.

Another way is to use the `config/config.yaml` file, which does not exist in the
repository and is also not tracked by git. But snakemake will always use this file if
it exists. This way you can run snakemake with a custom config without having to
specify the config file each time.

Configuration order of precedence is as follows:
1. Command line options specified with `--config` (optional)
2. Custom configuration file specified with `--configfile` (optional)
3. The `config/config.yaml` file (optional)
4. The default configuration files `config/config.default.yaml` and `config/plotting.default.yaml`

To use your custom configuration file, you need to pass it to the `snakemake` command
using the `--configfile` option:

```console
$ snakemake -call --configfile my_config.yaml
```


!!! warning
    In a previous version of PyPSA-Eur (`<=2025.04.0`), a full copy of the created config
    was stored in the `config/config.yaml` file. This is no longer the case. If the
    file exists, snakemake will use it, but no new copy will be created.


## Accessing configuration inside Snakemake

Rules should **not** access the `snakemake.config` object directly because overrides from
`run.scenarios` are only applied through the helpers in `rules/common.smk`:

- `config_provider("electricity", "extendable_carriers")` returns a callable
  that Snakemake evaluates per wildcard combination. This keeps caching fast and
  ensures the right scenario is used.
- `get_config(w)` materialises the fully merged dictionary for a specific set
  of wildcards. Use this sparingly inside Python helper functions that need to
  read several keys at once.

Reusing these helpers guarantees that documentation examples, rule
implementations, and custom extensions all observe the same precedence rules.


## `version` {#version_cf}

Version of PyPSA-Eur. Descriptive only.

- **Type:** string
- **Default:** `v2026.02.0`

**YAML Syntax**

```yaml
{{ yaml_section("version") }}
```


## `tutorial` {#tutorial_cf}

Switch to retrieve the tutorial data set instead of the full data set.

- **Type:** boolean
- **Default:** `false`

**YAML Syntax**

```yaml
{{ yaml_section("tutorial") }}
```


## `logging` {#logging_cf}

Configuration for top level `logging` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `level` | enum (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`) | `INFO` | Restrict console outputs to all infos, warning or errors only |
| `format` | string | `%(levelname)s:%(name)s:%(message)s` | Custom format for log messages. See [LogRecord ](https://docs.python.org/3/library/logging.html#logging.LogRecord) attributes. |

**YAML Syntax**

```yaml
{{ yaml_section("logging") }}
```


## `remote` {#remote_cf}

"Remote" indicates the address of a server used for data exchange, often for clusters and data pushing/pulling.

Configuration for top level `remote` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `ssh` | string | `""` | Optionally specify the SSH of a remote cluster to be synchronized. |
| `path` | string | `""` | Optionally specify the file path within the remote cluster to be synchronized. |

**YAML Syntax**

```yaml
{{ yaml_section("remote") }}
```


## `run` {#run_cf}

It is common conduct to analyse energy system optimisation models for **multiple scenarios** for a variety of reasons,
e.g. assessing their sensitivity towards changing the temporal and/or geographical resolution or investigating how
investment changes as more ambitious greenhouse-gas emission reduction targets are applied.

The `run` section is used for running and storing scenarios with different configurations which are not covered by [wildcards](wildcards.md).
It determines the path at which resources, networks and results are stored.
Therefore the user can run different configurations within the same directory.

Configuration for top level `run` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `prefix` | string | `""` | Prefix for the run name which is used as a top-layer directory name in the results and resources folders. |
| `name` | string \| list of string | `""` | Specify a name for your run. Results will be stored under this name. If `scenario: enable:` is set to `true`, the name must contain a subset of scenario names defined in `scenario: file:`. If the name is 'all', all defined scenarios will be run. |
| `scenarios` | any |  | Configuration for `run.scenarios` level. |
| ↳ `enable` | boolean | `false` | Switch to select whether workflow should generate scenarios based on `file`. |
| ↳ `file` | string | `config/scenarios.yaml` | Path to the scenario yaml file. The scenario file contains config overrides for each scenario. In order to be taken account, `run: scenarios` has to be set to `true` and `run: name` has to be a subset of top level keys given in the scenario file. In order to automatically create a `scenario.yaml` file based on a combination of settings, alter and use the `config/create_scenarios.py` script in the `config` directory. |
| `disable_progressbar` | boolean | `false` | Switch to select whether progressbar should be disabled. |
| `shared_resources` | any |  | Configuration for `run.shared_resources` level. |
| ↳ `policy` | boolean \| string | `false` | Boolean switch to select whether resources should be shared across runs. If a string is passed, this is used as a subdirectory name for shared resources. If set to 'base', only resources before creating the elec.nc file are shared. |
| ↳ `exclude` | list of string |  | For the case shared_resources=base, specify additional files that should not be shared across runs. |
| `use_shadow_directory` | boolean | `false` | Set to `true` (default) if snakemake shadow directories (`shallow`) should be used. Set to `false` if problems occur. |

**YAML Syntax**

```yaml
{{ yaml_section("run") }}
```


## `foresight` {#foresight_cf}

[planning_horizons](#planning-horizons) has to be set.

Configuration for `foresight` settings.

- **Type:** enum (`overnight`, `myopic`, `perfect`)
- **Default:** `overnight`

**YAML Syntax**

```yaml
{{ yaml_section("foresight") }}
```

!!! note
    If you use myopic or perfect foresight, define at least two values in the
    top-level [planning_horizons](#planning-horizons) list.

!!! note
    The `foresight` setting cannot vary across scenarios defined in
    `run.scenarios`. It is evaluated at workflow parsing time
    to determine which outputs to include. If you need to compare different
    foresight modes, run them as separate workflows with distinct `run.name`.


## `planning_horizons` {#planning-horizons}

Configure planning horizons at the top level rather than through wildcards.
Provide either a single year (for overnight studies) or a list of investment
years that should be simulated sequentially:

```yaml
planning_horizons: [2030, 2040, 2050]
```

Configuration for top level `planning_horizons` settings.

- **Type:** list of integer

**YAML Syntax**

```yaml
{{ yaml_section("planning_horizons") }}
```

- Overnight runs require a single value.
- Myopic runs expect strictly ascending values and continue each horizon from
  `RESULTS/networks/solved_{previous}.nc`.
- Perfect foresight also iterates over the list but reuses
  `networks/composed_{previous}.nc` as the brownfield seed.

!!! note
    Earlier releases derived planning horizons from `scenario` wildcard
    entries. That block is ignored now; define `planning_horizons` at the top
    level and keep scenario sweeps inside `run.scenarios`. See
    [migration](migration.md) for detailed conversion steps.


## `countries` {#countries}

Configuration for `countries` settings.

- **Type:** list of string
- **Default:** `[...]`

**YAML Syntax**

```yaml
{{ yaml_section("countries") }}
```


## `snapshots` {#snapshots_cf}

Specifies the temporal range to build an energy system model for as arguments to [pandas.date_range ](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.date_range.html)

Configuration for `snapshots` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `start` | string \| list of string | `2013-01-01` | Left bound of date range. |
| `end` | string \| list of string | `2014-01-01` | Right bound of date range. |
| `inclusive` | enum (`left`, `right`, `both`) \| null | `left` | Make the time interval closed to the `left`, `right`, or both sides `both` or neither side `None`. |

**YAML Syntax**

```yaml
{{ yaml_section("snapshots") }}
```


## `enable` {#enable_cf}

Switches for some rules and optional features.

Configuration for `enable` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `drop_leap_day` | boolean | `true` | Switch to drop February 29 from all time-dependent data in leap years. |

**YAML Syntax**

```yaml
{{ yaml_section("enable") }}
```


## `co2 budget` {#CO2_budget_cf}

Carbon budgets share one schema for all foresight modes. The `relative` flag
selects whether yearly entries inside `upper`/`lower` are interpreted as
fractions of the 1990 baseline (`true`) or absolute GtCO₂/year
(`false`). Enable `upper` and/or `lower` to enforce those caps only for
the explicitly listed years or a total budget across all [planning_horizons](#planning-horizons).

Configuration for `co2_budget` settings.

- **Type:** dict (str -> number)

**YAML Syntax**

```yaml
{{ yaml_section("co2_budget") }}
```


## `electricity` {#electricity_cf}

Configuration for `electricity` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `voltages` | list of number |  | Voltage levels to consider. |
| `base_network` | enum (`entsoegridkit`, `osm`, `tyndp`) | `osm` | Specify the underlying base network, i.e. GridKit (based on ENTSO-E web map extract), OpenStreetMap (OSM), or TYNDP. |
| `gaslimit_enable` | boolean | `false` | Add an overall absolute gas limit configured in `electricity: gaslimit`. |
| `gaslimit` | number \| boolean | `false` | Global gas usage limit. |
| `operational_reserve` | any |  | Configuration for `electricity.operational_reserve` settings. |
| ↳ `activate` | boolean | `false` | Whether to take operational reserve requirements into account during optimisation. |
| ↳ `epsilon_load` | number | `0.02` | share of total load. |
| ↳ `epsilon_vres` | number | `0.02` | share of total renewable supply. |
| ↳ `contingency` | number | `4000` | Fixed reserve capacity (MW). |
| `max_hours` | any |  | Configuration for `electricity.max_hours` settings. |
| ↳ `battery` | number | `6` | Maximum state of charge capacity of the battery in terms of hours at full output capacity `p_nom`. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| ↳ `li-ion` | number | `6` | Maximum state of charge capacity of the lithium-ion storage in terms of hours at full output capacity `p_nom`. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| ↳ `lfp` | number | `6` | Maximum state of charge capacity of the lithium-ion-LFP storage in terms of hours at full output capacity `p_nom`. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| ↳ `vanadium` | number | `10` | Maximum state of charge capacity of the vanadium-redox-flow storage in terms of hours at full output capacity `p_nom`. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| ↳ `lair` | number | `12` | Maximum state of charge capacity of the liquid-air storage in terms of hours at full output capacity `p_nom`. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| ↳ `pair` | number | `24` | Maximum state of charge capacity of the compressed-air-adiabatic storage in terms of hours at full output capacity `p_nom`. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| ↳ `iron-air` | number | `100` | Maximum state of charge capacity of the iron-air storage in terms of hours at full output capacity `p_nom`. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| ↳ `H2` | number | `168` | Maximum state of charge capacity of the hydrogen storage in terms of hours at full output capacity `p_nom`. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| `extendable_carriers` | any |  | Configuration for `electricity.extendable_carriers` settings. |
| ↳ `Generator` | list of string |  | Defines existing or non-existing conventional and renewable power plants to be extendable during the optimization. Conventional generators can only be built/expanded where already existent today. If a listed conventional carrier is not included in the `conventional_carriers` list, the lower limit of the capacity expansion is set to 0. |
| ↳ `StorageUnit` | list of string |  | Adds extendable storage units at every node/bus after clustering without capacity limits and with zero initial capacity. Supported technologies include battery, H2, li-ion, vanadium, lfp, lair, pair, and iron-air. |
| ↳ `Store` | list of string |  | Adds extendable storage units at every node/bus after clustering without capacity limits and with zero initial capacity. Supported technologies include battery, H2, li-ion, vanadium, lfp, lair, pair, and iron-air. |
| ↳ `Link` | list of string |  | Adds extendable links (H2 pipelines only) at every connection where there are lines or HVDC links without capacity limits and with zero initial capacity. Hydrogen pipelines require hydrogen storage to be modelled as `Store`. |
| `powerplants_filter` | string \| boolean | `(DateOut > 2025 or DateOut != DateOut) and (DateIn < 2026 or DateIn != DateIn)` | Filter query for the default powerplant database. |
| `custom_powerplants` | string \| boolean | `false` | Filter query for the custom powerplant database. |
| `everywhere_powerplants` | list of string |  | List of conventional power plants to add to every node in the model with zero initial capacity. To be used in combination with `extendable_carriers` to allow for building conventional powerplants irrespective of existing locations. |
| `conventional_carriers` | list of string |  | List of conventional power plants to include in the model from `resources/powerplants_s_{clusters}.csv`. If an included carrier is also listed in `extendable_carriers`, the capacity is taken as a lower bound. |
| `renewable_carriers` | list of string |  | List of renewable generators to include in the model. |
| `estimate_renewable_capacities` | any |  | Configuration for `electricity.estimate_renewable_capacities` settings. |
| ↳ `enable` | boolean | `true` | Activate routine to estimate renewable capacities in rule `add_electricity`. This option should not be used in combination with pathway planning `foresight: myopic` or `foresight: perfect` as renewable capacities are added differently in `add_existing_baseyear`. |
| ↳ `from_powerplantmatching` | boolean | `true` | Add renewable capacities from powerplantmatching dataset. |
| ↳ `from_irenastat` | boolean | `false` | Supplement powerplantmatching dataset with heuristics based on country-level renewable capacities from IRENA (IRENASTAT). |
| ↳ `year` | integer | `2024` | Renewable capacities are based on existing capacities reported by IRENA (IRENASTAT) for the specified year. |
| ↳ `expansion_limit` | number \| boolean | `false` | Artificially limit maximum IRENA capacities to a factor. For example, an `expansion_limit: 1.1` means 110% of capacities. If false are chosen, the estimated renewable potentials determine by the workflow are used. |
| ↳ `technology_mapping` | any |  | Configuration for `electricity.estimate_renewable_capacities.technology_mapping` settings. |
| ↳↳ `Offshore` | string | `offwind-ac` | PyPSA-Eur carrier that is considered for existing offshore wind technology (IRENA, GEM). |
| ↳↳ `Onshore` | string | `onwind` | PyPSA-Eur carrier that is considered for existing onshore wind capacities (IRENA, GEM). |
| ↳↳ `PV` | string | `solar` | PyPSA-Eur carrier that is considered for existing solar PV capacities (IRENA, GEM). |
| `estimate_battery_capacities` | boolean | `false` | Enable estimation of existing battery storage capacities. |
| `autarky` | any |  | Configuration for `electricity.autarky` settings. |
| ↳ `enable` | boolean | `false` | Require each node to be autarkic by removing all lines and links. |
| ↳ `by_country` | boolean | `false` | Require each country to be autarkic by removing all cross-border lines and links. `electricity: autarky` must be enabled. |
| `transmission_limit` | string | `vopt` | Limit on transmission expansion. The first part can be `v` (for setting a limit on line volume) or `c` (for setting a limit on line cost). The second part can be `opt` or a float bigger than one (e.g. 1.25). If `opt` is chosen line expansion is optimised according to its capital cost (where the choice `v` only considers overhead costs for HVDC transmission lines, while `c` uses more accurate costs distinguishing between overhead and underwater sections and including inverter pairs). The setting `v1.25` will limit the total volume of line expansion to 25% of currently installed capacities weighted by individual line lengths. The setting `c1.25` will allow to build a transmission network that costs no more than 25 % more than the current system. |

**YAML Syntax**

```yaml
{{ yaml_section("electricity") }}
```


## `atlite` {#atlite_cf}

Define and specify the `atlite.Cutout` used for calculating renewable potentials and time-series. All options except for `features` are directly used as [cutout parameters ](https://atlite.readthedocs.io/en/latest/ref_api.html#cutout).

Configuration for `atlite` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `default_cutout` | string \| list of string | `europe-2013-sarah3-era5` | Defines a default cutout. Can refer to a single cutout or a list of cutouts. |
| `nprocesses` | integer | `1` | Number of parallel processes in cutout preparation. |
| `show_progress` | boolean | `false` | Whether progressbar for atlite conversion processes should be shown. False saves time. |
| `plot_availability_matrix` | boolean | `false` | Whether to plot the landuse availability matrix. |
| `cutouts` | dict (str -> any) |  | Named cutout configurations. |

**YAML Syntax**

```yaml
{{ yaml_section("atlite") }}
```


## `renewable` {#renewable_cf}

### `onwind`

Configuration for onshore wind.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `cutout` | string \| list of string | `default` | Specifies the weather data cutout file(s) to use. |
| `resource` | any |  | Configuration for wind resource settings. |
| ↳ `method` | string | `wind` | A superordinate technology type. |
| ↳ `turbine` | string \| dict (str -> string) |  | Specifies the turbine type and its characteristic power curve. Can be a string or a dictionary with years as keys which denote the year another turbine model becomes available. |
| ↳ `smooth` | boolean | `false` | Switch to apply a gaussian kernel density smoothing to the power curve. |
| ↳ `add_cutout_windspeed` | boolean | `true` | Whether to add cutout windspeed data. |
| `resource_classes` | integer | `1` | Number of resource classes per clustered region. |
| `capacity_per_sqkm` | number | `3` | Allowable density of wind turbine placement. |
| `correction_factor` | number | `1.0` | Correction factor for capacity factor time series. |
| `corine` | boolean \| any |  | CORINE land cover configuration. |
| ↳ `grid_codes` | list of integer |  | Specifies areas according to CORINE Land Cover codes which are generally eligible for wind turbine placement. |
| ↳ `distance` | number | `1000` | Distance in meters to keep from areas specified in `distance_grid_codes`. |
| ↳ `distance_grid_codes` | list of integer |  | Specifies areas according to CORINE Land Cover codes to which wind turbines must maintain a distance specified in the setting `distance`. |
| `luisa` | boolean \| object | `false` | LUISA land cover configuration. |
| `natura` | boolean | `true` | Switch to exclude [Natura 2000 ](https://en.wikipedia.org/wiki/Natura_2000) natural protection areas. Area is excluded if `true`. |
| `excluder_resolution` | number | `100` | Resolution in meters on which to perform geographical eligibility analysis. |
| `clip_p_max_pu` | number | `0.01` | To avoid too small values in the renewables` per-unit availability time series values below this threshold are set to zero. |

Configuration for offshore wind.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `cutout` | string \| list of string | `default` | Specifies the weather data cutout file(s) to use. |
| `resource` | any |  | Configuration for wind resource settings. |
| ↳ `method` | string | `wind` | A superordinate technology type. |
| ↳ `turbine` | string \| dict (str -> string) |  | Specifies the turbine type and its characteristic power curve. Can be a string or a dictionary with years as keys which denote the year another turbine model becomes available. |
| ↳ `smooth` | boolean | `false` | Switch to apply a gaussian kernel density smoothing to the power curve. |
| ↳ `add_cutout_windspeed` | boolean | `true` | Whether to add cutout windspeed data. |
| `resource_classes` | integer | `1` | Number of resource classes per clustered region. |
| `capacity_per_sqkm` | number | `2` | Allowable density of wind turbine placement. |
| `correction_factor` | number | `0.8855` | Correction factor for capacity factor time series. |
| `corine` | boolean \| list of integer | `false` | Specifies areas according to CORINE Land Cover codes which are generally eligible for AC-connected offshore wind turbine placement. |
| `luisa` | boolean \| list of integer | `false` | Specifies areas according to the LUISA Base Map codes which are generally eligible for AC-connected offshore wind turbine placement. |
| `natura` | boolean | `true` | Switch to exclude [Natura 2000 ](https://en.wikipedia.org/wiki/Natura_2000) natural protection areas. Area is excluded if `true`. |
| `ship_threshold` | number | `400` | Ship density threshold from which areas are excluded. |
| `max_depth` | number \| null |  | Maximum sea water depth in meters at which wind turbines can be built. Maritime areas with deeper waters are excluded in the process of calculating the AC-connected offshore wind potential. |
| `min_depth` | number \| null |  | Minimum water depth in meters. |
| `max_shore_distance` | number \| null |  | Maximum distance to the shore in meters above which wind turbines cannot be built. Such areas are excluded in the process of calculating the AC-connected offshore wind potential. |
| `min_shore_distance` | number \| null |  | Minimum distance to the shore in meters below which wind turbines cannot be built. Such areas close to the shore are excluded in the process of calculating the AC-connected offshore wind potential. |
| `excluder_resolution` | number | `200` | Resolution in meters on which to perform geographical eligibility analysis. |
| `clip_p_max_pu` | number | `0.01` | To avoid too small values in the renewables` per-unit availability time series values below this threshold are set to zero. |
| `landfall_length` | number \| string | `20` | Fixed length of the cable connection that is onshorelandfall in km. If 'centroid', the length is calculated as the distance to centroid of the onshore bus. |

Configuration for offshore wind.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `cutout` | string \| list of string | `default` | Specifies the weather data cutout file(s) to use. |
| `resource` | any |  | Configuration for wind resource settings. |
| ↳ `method` | string | `wind` | A superordinate technology type. |
| ↳ `turbine` | string \| dict (str -> string) |  | Specifies the turbine type and its characteristic power curve. Can be a string or a dictionary with years as keys which denote the year another turbine model becomes available. |
| ↳ `smooth` | boolean | `false` | Switch to apply a gaussian kernel density smoothing to the power curve. |
| ↳ `add_cutout_windspeed` | boolean | `true` | Whether to add cutout windspeed data. |
| `resource_classes` | integer | `1` | Number of resource classes per clustered region. |
| `capacity_per_sqkm` | number | `2` | Allowable density of wind turbine placement. |
| `correction_factor` | number | `0.8855` | Correction factor for capacity factor time series. |
| `corine` | boolean \| list of integer | `false` | Specifies areas according to CORINE Land Cover codes which are generally eligible for AC-connected offshore wind turbine placement. |
| `luisa` | boolean \| list of integer | `false` | Specifies areas according to the LUISA Base Map codes which are generally eligible for AC-connected offshore wind turbine placement. |
| `natura` | boolean | `true` | Switch to exclude [Natura 2000 ](https://en.wikipedia.org/wiki/Natura_2000) natural protection areas. Area is excluded if `true`. |
| `ship_threshold` | number | `400` | Ship density threshold from which areas are excluded. |
| `max_depth` | number \| null |  | Maximum sea water depth in meters at which wind turbines can be built. Maritime areas with deeper waters are excluded in the process of calculating the AC-connected offshore wind potential. |
| `min_depth` | number \| null |  | Minimum water depth in meters. |
| `max_shore_distance` | number \| null |  | Maximum distance to the shore in meters above which wind turbines cannot be built. Such areas are excluded in the process of calculating the AC-connected offshore wind potential. |
| `min_shore_distance` | number \| null |  | Minimum distance to the shore in meters below which wind turbines cannot be built. Such areas close to the shore are excluded in the process of calculating the AC-connected offshore wind potential. |
| `excluder_resolution` | number | `200` | Resolution in meters on which to perform geographical eligibility analysis. |
| `clip_p_max_pu` | number | `0.01` | To avoid too small values in the renewables` per-unit availability time series values below this threshold are set to zero. |
| `landfall_length` | number \| string | `20` | Fixed length of the cable connection that is onshorelandfall in km. If 'centroid', the length is calculated as the distance to centroid of the onshore bus. |

Configuration for offshore wind.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `cutout` | string \| list of string | `default` | Specifies the weather data cutout file(s) to use. |
| `resource` | any |  | Configuration for wind resource settings. |
| ↳ `method` | string | `wind` | A superordinate technology type. |
| ↳ `turbine` | string \| dict (str -> string) |  | Specifies the turbine type and its characteristic power curve. Can be a string or a dictionary with years as keys which denote the year another turbine model becomes available. |
| ↳ `smooth` | boolean | `false` | Switch to apply a gaussian kernel density smoothing to the power curve. |
| ↳ `add_cutout_windspeed` | boolean | `true` | Whether to add cutout windspeed data. |
| `resource_classes` | integer | `1` | Number of resource classes per clustered region. |
| `capacity_per_sqkm` | number | `2` | Allowable density of wind turbine placement. |
| `correction_factor` | number | `0.8855` | Correction factor for capacity factor time series. |
| `corine` | boolean \| list of integer | `false` | Specifies areas according to CORINE Land Cover codes which are generally eligible for AC-connected offshore wind turbine placement. |
| `luisa` | boolean \| list of integer | `false` | Specifies areas according to the LUISA Base Map codes which are generally eligible for AC-connected offshore wind turbine placement. |
| `natura` | boolean | `true` | Switch to exclude [Natura 2000 ](https://en.wikipedia.org/wiki/Natura_2000) natural protection areas. Area is excluded if `true`. |
| `ship_threshold` | number | `400` | Ship density threshold from which areas are excluded. |
| `max_depth` | number \| null |  | Maximum sea water depth in meters at which wind turbines can be built. Maritime areas with deeper waters are excluded in the process of calculating the AC-connected offshore wind potential. |
| `min_depth` | number \| null |  | Minimum water depth in meters. |
| `max_shore_distance` | number \| null |  | Maximum distance to the shore in meters above which wind turbines cannot be built. Such areas are excluded in the process of calculating the AC-connected offshore wind potential. |
| `min_shore_distance` | number \| null |  | Minimum distance to the shore in meters below which wind turbines cannot be built. Such areas close to the shore are excluded in the process of calculating the AC-connected offshore wind potential. |
| `excluder_resolution` | number | `200` | Resolution in meters on which to perform geographical eligibility analysis. |
| `clip_p_max_pu` | number | `0.01` | To avoid too small values in the renewables` per-unit availability time series values below this threshold are set to zero. |
| `landfall_length` | number \| string | `20` | Fixed length of the cable connection that is onshorelandfall in km. If 'centroid', the length is calculated as the distance to centroid of the onshore bus. |

Configuration for solar PV.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `cutout` | string \| list of string | `default` | Specifies the weather data cutout file(s) to use. |
| `resource` | any |  | Configuration for solar resource settings. |
| ↳ `method` | string | `pv` | A superordinate technology type. |
| ↳ `panel` | string \| dict (str -> string) | `CSi` | Specifies the solar panel technology and its characteristic attributes. Can be a string or a dictionary with years as keys which denote the year another panel model becomes available. |
| ↳ `orientation` | dict (str -> number) |  | Panel orientation with slope and azimuth. |
| ↳ `tracking` | string \| null |  | Tracking type (e.g., 'horizontal'). |
| `resource_classes` | integer | `1` | Number of resource classes per clustered region. |
| `capacity_per_sqkm` | number | `5.1` | Allowable density of solar panel placement. |
| `correction_factor` | number | `1.0` | A correction factor for the capacity factor (availability) time series. |
| `corine` | boolean \| list of integer |  | Specifies areas according to CORINE Land Cover codes which are generally eligible for solar panel placement. |
| `luisa` | boolean \| list of integer | `false` | Specifies areas according to the LUISA Base Map codes which are generally eligible for solar panel placement. |
| `natura` | boolean | `true` | Switch to exclude [Natura 2000 ](https://en.wikipedia.org/wiki/Natura_2000) natural protection areas. Area is excluded if `true`. |
| `excluder_resolution` | number | `100` | Resolution in meters on which to perform geographical eligibility analysis. |
| `clip_p_max_pu` | number | `0.01` | To avoid too small values in the renewables` per-unit availability time series values below this threshold are set to zero. |

Configuration for hydropower.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `cutout` | string \| list of string | `default` | Specifies the weather data cutout file(s) to use. |
| `carriers` | list of string |  | Specifies the types of hydro power plants to build per-unit availability time series for. 'ror' stands for run-of-river plants, 'PHS' represents pumped-hydro storage, and 'hydro' stands for hydroelectric dams. |
| `PHS_max_hours` | number | `6` | Maximum state of charge capacity of the pumped-hydro storage (PHS) in terms of hours at full output capacity `p_nom`. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| `hydro_max_hours` | string \| number | `energy_capacity_totals_by_country` | Maximum state of charge capacity of the pumped-hydro storage (PHS) in terms of hours at full output capacity `p_nom` or heuristically determined. Cf. [PyPSA documentation ](https://pypsa.readthedocs.io/en/latest/components.html#storage-unit). |
| `flatten_dispatch` | boolean | `false` | Consider an upper limit for the hydro dispatch. The limit is given by the average capacity factor plus the buffer given in `flatten_dispatch_buffer`. |
| `flatten_dispatch_buffer` | number | `0.2` | If `flatten_dispatch` is true, specify the value added above the average capacity factor. |
| `clip_min_inflow` | number | `1.0` | To avoid too small values in the inflow time series, values below this threshold (MW) are set to zero. |
| `eia_norm_year` | boolean \| integer | `false` | To specify a specific year by which hydro inflow is normed that deviates from the snapshots' year. |
| `eia_correct_by_capacity` | boolean | `false` | Correct EIA annual hydro generation data by installed capacity. |
| `eia_approximate_missing` | boolean | `false` | Approximate hydro generation data for years not included in EIA dataset through a regression based on annual runoff. |

**YAML Syntax**

```yaml
{{ yaml_section("renewable") }}
```

**YAML Syntax**

```yaml
{{ yaml_section("renewable.offwind-ac") }}
```

**YAML Syntax**

```yaml
{{ yaml_section("renewable.solar") }}
```

**YAML Syntax**

```yaml
{{ yaml_section("renewable.hydro") }}
```

!!! note
    Notes on `capacity_per_sqkm`. ScholzPhd Tab 4.3.1: 10MW/km^2 and assuming 30% fraction of the already restricted
       area is available for installation of wind generators due to competing land use and likely public
       acceptance issues.

!!! note
    The default choice for corine `grid_codes` was based on Scholz, Y. (2012). Renewable energy based electricity supply at low costs
    development of the REMix model and application for Europe. ( p.42 / p.28)


!!! note
    Notes on `capacity_per_sqkm`. ScholzPhd Tab 4.3.1: 10MW/km^2 and assuming 20% fraction of the already restricted
       area is available for installation of wind generators due to competing land use and likely public
       acceptance issues.

!!! note
    Notes on `correction_factor`. Correction due to proxy for wake losses
    from 10.1016/j.energy.2018.08.153
    until done more rigorously in #153


!!! note
    Notes on `capacity_per_sqkm`. ScholzPhd Tab 4.3.1: 170 MW/km^2 and assuming 1% of the area can be used for solar PV panels.
       Correction factor determined by comparing uncorrected area-weighted full-load hours to those
       published in Supplementary Data to Pietzcker, Robert Carl, et al. "Using the sun to decarbonize the power
       sector -- The economic potential of photovoltaics and concentrating solar
       power." Applied Energy 135 (2014): 704-720.
       This correction factor of 0.854337 may be in order if using reanalysis data.
       for discussion refer to this https://github.com/PyPSA/pypsa-eur/issues/285


## `conventional` {#conventional_cf}

Define additional generator attribute for conventional carrier types. If a
scalar value is given it is applied to all generators. However if a string
starting with "data/" is given, the value is interpreted as a path to a csv file
with country specific values. Then, the values are read in and applied to all
generators of the given carrier in the given country. Note that the value(s)
overwrite the existing values.

Configuration for `conventional` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `unit_commitment` | boolean | `false` | Allow the overwrite of ramp_limit_up, ramp_limit_start_up, ramp_limit_shut_down, p_min_pu, min_up_time, min_down_time, and start_up_cost of conventional generators. Refer to the CSV file 'unit_commitment.csv'. |
| `dynamic_fuel_price` | boolean | `false` | Consider the monthly fluctuating fuel prices for each conventional generator. Refer to the CSV file 'data/validation/monthly_fuel_price.csv'. |
| `fuel_price_rolling_window` | integer | `6` | Monthly rolling mean window for fossil fuel prices smoothing. |
| `nuclear` | dict (str -> string \| number) |  | For any carrier/technology overwrite attributes as listed below. |

**YAML Syntax**

```yaml
{{ yaml_section("conventional") }}
```


## `lines` {#lines_cf}

Configuration for `lines` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `types` | dict (str -> string) |  | Specifies line types to assume for the different voltage levels of the ENTSO-E grid extraction. Should normally handle voltage levels 220, 300, and 380 kV. |
| `s_max_pu` | number | `0.7` | Correction factor for line capacities (`s_nom`) to approximate N-1 security and reserve capacity for reactive power flows. |
| `s_nom_max` | number |  | Global upper limit for the maximum capacity of each extendable line (MW). |
| `max_extension` | number | `20000` | Upper limit for the extended capacity of each extendable line (MW). |
| `length_factor` | number | `1.25` | Correction factor to account for the fact that buses are *not* connected by lines through air-line distance. |
| `reconnect_crimea` | boolean | `true` | Whether to reconnect Crimea to the Ukrainian grid. |
| `under_construction` | enum (`zero`, `remove`, `keep`) | `keep` | Specifies how to handle lines which are currently under construction. |
| `dynamic_line_rating` | any |  | Configuration for `lines.dynamic_line_rating` settings. |
| ↳ `activate` | boolean | `false` | Whether to take dynamic line rating into account. |
| ↳ `cutout` | string \| list of string | `default` | Specifies the weather data cutout file(s) to use. |
| ↳ `correction_factor` | number | `0.95` | Factor to compensate for overestimation of wind speeds in hourly averaged wind data. |
| ↳ `max_voltage_difference` | number \| `False` | `false` | Maximum voltage angle difference in degrees or 'false' to disable. |
| ↳ `max_line_rating` | number \| `False` | `false` | Maximum line rating relative to nominal capacity without DLR, e.g. 1.3 or 'false' to disable. |

**YAML Syntax**

```yaml
{{ yaml_section("lines") }}
```


## `links` {#links_cf}

Configuration for `links` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `p_max_pu` | number | `1.0` | Correction factor for link capacities `p_nom`. |
| `p_min_pu` | number | `-1.0` | Correction factor for link capacities `p_nom`. |
| `p_nom_max` | number |  | Global upper limit for the maximum capacity of each extendable DC link (MW). |
| `max_extension` | number | `30000` | Upper limit for the extended capacity of each extendable DC link (MW). |
| `length_factor` | number | `1.25` | Correction factor to account for the fact that buses are *not* connected by links through air-line distance. |
| `under_construction` | enum (`zero`, `remove`, `keep`) | `keep` | Specifies how to handle lines which are currently under construction. |

**YAML Syntax**

```yaml
{{ yaml_section("links") }}
```


## `transmission_projects` {#transmission_projects_cf}

Allows to define additional transmission projects that will be added to the base network, e.g., from the TYNDP 2020 dataset. The projects are read in from the CSV files in the subfolder of `data/transmission_projects/`. New transmission projects, e.g. from TYNDP 2024, can be added in a new subfolder of transmission projects, e.g. `data/transmission_projects/tyndp2024` while extending the list of `transmission_projects` in the `config.yaml` by `tyndp2024`. The CSV files in the project folder should have the same columns as the CSV files in the template folder `data/transmission_projects/template`.

Configuration for `transmission_projects` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `enable` | boolean | `true` | Whether to integrate this transmission projects or not. |
| `include` | any |  | Configuration for `transmission_projects.include` settings. |
| ↳ `tyndp2020` | boolean | `true` | Whether to integrate the TYNDP 2020 dataset. |
| ↳ `nep` | boolean | `true` | Whether to integrate the German network development plan dataset. |
| ↳ `manual` | boolean | `true` | Whether to integrate the manually added transmission projects. They are taken from the previously existing links_tyndp.csv file. |
| `skip` | list of string |  | Type of lines to skip from all transmission projects. Possible values are: `upgraded_lines`, `upgraded_links`, `new_lines`, `new_links`. |
| `status` | list of string \| dict (str -> list of string) |  | Status to include into the model as list or as dict with name of project and status to include. Possible values for status are `under_construction`, `in_permitting`, `confirmed`, `planned_not_yet_permitted`, `under_consideration`. |
| `new_link_capacity` | enum (`zero`, `keep`) | `zero` | Whether to set the new link capacity to the provided capacity or set it to zero. |

**YAML Syntax**

```yaml
{{ yaml_section("transmission_projects") }}
```


## `transformers` {#transformers_cf}

Configuration for `transformers` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `x` | number | `0.1` | Series reactance in per unit (p.u.), using `s_nom` as base power of the transformer. Overwritten if `type` is specified. |
| `s_nom` | number | `2000.0` | Limit of apparent power which can pass through branch (MVA). Overwritten if `type` is specified. |
| `type` | string | `""` | Specifies transformer types to assume for the transformers of the ENTSO-E grid extraction. |

**YAML Syntax**

```yaml
{{ yaml_section("transformers") }}
```


## `load` {#load_cf}

Configuration for `load` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `fill_gaps` | any |  | Configuration for `load.fill_gaps` settings. |
| ↳ `enable` | boolean | `true` | Whether to fill gaps using interpolation for small gaps and time shift for large gaps. |
| ↳ `interpolate_limit` | integer | `6` | Maximum gap size (consecutive nans) which interpolated linearly. |
| ↳ `time_shift_for_large_gaps` | string | `1w` | Periods which are used for copying time-slices in order to fill large gaps of nans. Have to be valid `pandas` period strings. |
| `manual_adjustments` | boolean | `true` | Whether to adjust the load data manually according to the function in `manual_adjustment`. |
| `scaling_factor` | number | `1.0` | Global correction factor for the load time series. |
| `fixed_year` | integer \| boolean | `false` | To specify a fixed year for the load time series that deviates from the snapshots' year. |
| `supplement_synthetic` | boolean | `true` | Whether to supplement missing data for selected time period should be supplemented by synthetic data from [Zenodo ](https://zenodo.org/records/10820928). |
| `substation_only` | boolean | `true` | Whether to only consider substations for the spatial disaggregation of the per-country electricity demand data. |
| `distribution_key` | any |  | Configuration for `load.distribution_key` settings. |
| ↳ `gdp` | number | `0.6` | Weighting factor for the GDP data in the distribution key. |
| ↳ `population` | number | `0.4` | Weighting factor for the population data in the distribution key. |

**YAML Syntax**

```yaml
{{ yaml_section("load") }}
```


## `energy` {#energy_cf}

Only used for sector-coupling studies.

Configuration for `energy` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `energy_totals_year` | integer | `2023` | The year for the sector energy use. The year must be available in the Eurostat report. |
| `base_emissions_year` | integer | `1990` | The base year for the sector emissions. See [European Environment Agency (EEA) ](https://www.eea.europa.eu/data-and-maps/data/national-emissions-reported-to-the-unfccc-and-to-the-eu-greenhouse-gas-monitoring-mechanism-16). |
| `emissions` | string | `CO2` | Specify which sectoral emissions are taken into account. Data derived from EEA. Currently only CO2 is implemented. |

**YAML Syntax**

```yaml
{{ yaml_section("energy") }}
```

!!! note
    Only used for sector-coupling studies.


## `biomass` {#biomass_cf}

- Manure solid, liquid
- Residues from landscape care
- Bioethanol barley, wheat, grain maize, oats, other cereals and rye
- Sugar from sugar beet
- Miscanthus, switchgrass, RCG
- Willow
- Poplar
- Sunflower, soya seed
- Rape seed
- Fuelwood residues
- FuelwoodRW
- C&P_RW
- Secondary Forestry residues - woodchips
- Sawdust
- Municipal waste
- Sludge

Configuration for `biomass` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `year` | integer | `2030` | Year for which to retrieve biomass potential according to the assumptions of the [JRC ENSPRESO ](https://data.jrc.ec.europa.eu/dataset/74ed5a04-7d74-4807-9eab-b94774309d9f). |
| `scenario` | enum (`ENS_Low`, `ENS_Med`, `ENS_High`) | `ENS_Med` | Scenario for which to retrieve biomass potential. The scenario definition can be seen in [ENSPRESO_BIOMASS ](https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/ENSPRESO/ENSPRESO_BIOMASS.xlsx). |
| `classes` | any |  | Configuration for `biomass.classes` settings. |
| ↳ `solid biomass` | list of string |  | The comodity that are included as solid biomass. |
| ↳ `not included` | list of string |  | The comodity that are not included as a biomass potential. |
| ↳ `biogas` | list of string |  | The comodity that are included as biogas. |
| ↳ `municipal solid waste` | list of string |  | The commodities that are included as municipal solid waste. |
| `share_unsustainable_use_retained` | dict (str -> number) |  | Share of unsustainable biomass use retained using primary production of Eurostat data as reference. |
| `share_sustainable_potential_available` | dict (str -> number) |  | Share determines phase-in of ENSPRESO biomass potentials. |

**YAML Syntax**

```yaml
{{ yaml_section("biomass") }}
```

!!! note
    Only used for sector-coupling studies.

    The list of available biomass is given by the category in [ENSPRESO_BIOMASS](https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/ENSPRESO/ENSPRESO_BIOMASS.xlsx), namely:

    - Agricultural waste


## `solar_thermal` {#solar_thermal_cf}

Only used for sector-coupling studies.

Configuration for `solar_thermal` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `clearsky_model` | enum (`simple`, `enhanced`) | `simple` | Type of clearsky model for diffuse irradiation. |
| `orientation` | any |  | Configuration for `solar_thermal.orientation` settings. |
| ↳ `slope` | number | `45.0` | The angle between the ground and the panels. |
| ↳ `azimuth` | number | `180.0` | The angle between the North and the sun with panels on the local horizon. |
| `cutout` | string | `default` | Name of the cutout to use for solar thermal calculations. |

**YAML Syntax**

```yaml
{{ yaml_section("solar_thermal") }}
```

!!! note
    Only used for sector-coupling studies.


## `existing_capacities` {#existing_capacities_cf}

Only used for sector-coupling studies. The value for grouping years are only used in myopic or perfect foresight scenarios.

Activating `enabled` instructs [compose_network][] to merge the historical
assets stored in `resources/powerplants.csv` into `networks/composed_{horizon}.nc`.

Configuration for `existing_capacities` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `grouping_years_power` | list of integer |  | Intervals to group existing capacities for power. |
| `grouping_years_heat` | list of integer |  | Intervals to group existing capacities for heat. |
| `threshold_capacity` | number | `10` | Capacities (MW) of generators and links below threshold are removed during add_existing_capacities. |
| `default_heating_lifetime` | integer | `20` | Default lifetime for heating technologies (years). |
| `solar_rooftop_ratio` | number | `0.5` | Ratio of existing solar capacity to assign to rooftop vs utility-scale (between 0 and 1). |
| `conventional_carriers` | list of string |  | List of conventional power plants to include in the sectoral network. |

**YAML Syntax**

```yaml
{{ yaml_section("existing_capacities") }}
```

!!! note
    Only used for sector-coupling studies. The value for grouping years are only used in myopic or perfect foresight scenarios.


## `sector` {#sector_cf}

Only used for sector-coupling studies.

??? note "Details"

    Configuration for `sector` settings.
    
    | Property | Type | Default | Description |
    |----------|------|---------|-------------|
    | `transport` | boolean | `true` | Flag to include transport sector. |
    | `heating` | boolean | `true` | Flag to include heating sector. |
    | `biomass` | boolean | `true` | Flag to include biomass sector. |
    | `industry` | boolean | `true` | Flag to include industry sector. |
    | `shipping` | boolean | `true` | Flag to include shipping sector. |
    | `aviation` | boolean | `true` | Flag to include aviation sector. |
    | `agriculture` | boolean | `true` | Flag to include agriculture sector. |
    | `fossil_fuels` | boolean | `true` | Flag to include imports of fossil fuels. |
    | `district_heating` | any |  | Configuration for `sector.district_heating` settings. |
    |   `potential` | number \| dict (str -> number) | `0.6` | Maximum fraction of urban demand which can be supplied by district heating. If given as dictionary, specify one value per country modeled or provide a default value with key `default` to fill values for all unspecified countries. |
    |   `progress` | dict (str -> number) |  | Increase of today's district heating demand to potential maximum district heating share. Progress = 0 means today's district heating share. Progress = 1 means maximum fraction of urban demand is supplied by district heating. |
    |   `district_heating_loss` | number | `0.15` | Share increase in district heat demand in urban central due to heat losses. |
    |   `supply_temperature_approximation` | object |  | Supply temperature approximation settings. |
    |   `ptes` | object |  | Pit thermal energy storage settings. |
    |   `ates` | object |  | Aquifer thermal energy storage settings. |
    |   `heat_source_cooling` | number | `6` | Cooling of heat source for heat pumps. |
    |   `heat_pump_cop_approximation` | object |  | Heat pump COP approximation settings. |
    |   `limited_heat_sources` | object |  | Dictionary with names of limited heat sources (not air). Must be `river_water` / `geothermal` or another heat source in [Manz et al. 2024 ](https://www.sciencedirect.com/science/article/pii/S0960148124001769). |
    |   `direct_utilisation_heat_sources` | list of string |  | List of heat sources for direct heat utilisation in district heating. Must be in the keys of `heat_utilisation_potentials` (e.g. `geothermal`). |
    |   `temperature_limited_stores` | list of string |  | List of names for stores used as limited heat sources. |
    |   `dh_areas` | object |  | District heating areas settings. |
    | `heat_pump_sources` | dict (str -> list of string) |  | Heat pump sources by area. |
    | `residential_heat` | any |  | Configuration for `sector.residential_heat` settings. |
    |   `dsm` | any |  | Configuration for `sector.residential_heat.dsm` settings. |
    |     `enable` | boolean | `false` | Enable residential heat demand-side management that allows heating systems to provide flexibility by shifting demand within configurable time periods. Models building thermal mass as energy storage. |
    |     `direction` | list of string |  | 'overheat-undercool' means both pre-heating and delayed heating are allowed. 'overheat' allows only pre-heating where buildings are heated up above target temperature and then allowed to cool down, while 'undercool' allows only delayed heating where buildings can cool below target temperature and then be heated up again. |
    |     `restriction_value` | dict (str -> number) |  | Maximum state of charge (as fraction) for heat flexibility storage representing available thermal buffer capacity in buildings. Set to 0 for no flexibility or to 1.0 to assume that the entire heating demand can contribute to flexibility. |
    |     `restriction_time` | list of integer |  | Checkpoint hours (0-23) at which heat flexibility storage must return to baseline state of charge, i.e. the residence surplus or missing heat be balanced. Time is the local time for each country and bus. Default: [10, 22] creates 12-hour periods with checkpoints at 10am and 10pm. |
    | `cluster_heat_buses` | boolean | `true` | Cluster residential and service heat buses in [prepare_sector_network.py ](https://github.com/PyPSA/pypsa-eur-sec/blob/master/scripts/prepare_sector_network.py) to one to save memory. |
    | `heat_demand_cutout` | string | `default` | Heat demand cutout. |
    | `bev_dsm_restriction_value` | number | `0.8` | Adds a lower state of charge (SOC) limit for battery electric vehicles (BEV) to manage its own energy demand (DSM). Located in [build_transport_demand.py ](https://github.com/PyPSA/pypsa-eur-sec/blob/master/scripts/build_transport_demand.py). Set to 0 for no restriction on BEV DSM. |
    | `bev_dsm_restriction_time` | number | `7` | Time at which SOC of BEV has to be dsm_restriction_value. |
    | `transport_heating_deadband_upper` | number | `20.0` | The maximum temperature in the vehicle. At higher temperatures, the energy required for cooling in the vehicle increases. |
    | `transport_heating_deadband_lower` | number | `15.0` | The minimum temperature in the vehicle. At lower temperatures, the energy required for heating in the vehicle increases. |
    | `ICE_lower_degree_factor` | number | `0.375` | Share increase in energy demand in internal combustion engine (ICE) for each degree difference between the cold environment and the minimum temperature. |
    | `ICE_upper_degree_factor` | number | `1.6` | Share increase in energy demand in internal combustion engine (ICE) for each degree difference between the hot environment and the maximum temperature. |
    | `EV_lower_degree_factor` | number | `0.98` | Share increase in energy demand in electric vehicles (EV) for each degree difference between the cold environment and the minimum temperature. |
    | `EV_upper_degree_factor` | number | `0.63` | Share increase in energy demand in electric vehicles (EV) for each degree difference between the hot environment and the maximum temperature. |
    | `bev_dsm` | boolean | `true` | Add the option for battery electric vehicles (BEV) to participate in demand-side management (DSM). |
    | `bev_dsm_availability` | number | `0.5` | The share for battery electric vehicles (BEV) that are able to do demand side management (DSM). |
    | `bev_energy` | number | `0.05` | The average size of battery electric vehicles (BEV) in MWh. |
    | `bev_charge_efficiency` | number | `0.9` | Battery electric vehicles (BEV) charge and discharge efficiency. |
    | `bev_charge_rate` | number | `0.011` | The power consumption for one electric vehicle (EV) in MWh. Value derived from 3-phase charger with 11 kW. |
    | `bev_avail_max` | number | `0.95` | The maximum share plugged-in availability for passenger electric vehicles. |
    | `bev_avail_mean` | number | `0.8` | The average share plugged-in availability for passenger electric vehicles. |
    | `v2g` | boolean | `true` | Allows feed-in to grid from EV battery. This is only enabled if BEV demand-side management is enabled, and the share of vehicles participating is V2G is given by `bev_dsm_availability`. |
    | `land_transport_fuel_cell_share` | dict (str -> number) |  | The share of vehicles that uses fuel cells in a given year. |
    | `land_transport_electric_share` | dict (str -> number) |  | The share of vehicles that uses electric vehicles (EV) in a given year. |
    | `land_transport_ice_share` | dict (str -> number) |  | The share of vehicles that uses internal combustion engines (ICE) in a given year. What is not EV or FCEV is oil-fuelled ICE. |
    | `transport_electric_efficiency` | number | `53.19` | The conversion efficiencies of electric vehicles in transport. |
    | `transport_fuel_cell_efficiency` | number | `30.003` | The H2 conversion efficiencies of fuel cells in transport. |
    | `transport_ice_efficiency` | number | `16.0712` | The oil conversion efficiencies of internal combustion engine (ICE) in transport. |
    | `agriculture_machinery_electric_share` | number | `0.5` | The share for agricultural machinery that uses electricity. |
    | `agriculture_machinery_oil_share` | number | `0.5` | The share for agricultural machinery that uses oil. |
    | `agriculture_machinery_fuel_efficiency` | number | `0.7` | The efficiency of electric-powered machinery in the conversion of electricity to meet agricultural needs. |
    | `agriculture_machinery_electric_efficiency` | number | `0.3` | The efficiency of oil-powered machinery in the conversion of oil to meet agricultural needs. |
    | `shipping_hydrogen_liquefaction` | boolean | `false` | Whether to include liquefaction costs for hydrogen demand in shipping. |
    | `shipping_hydrogen_share` | dict (str -> number) |  | The share of ships powered by hydrogen in a given year. |
    | `shipping_methanol_share` | dict (str -> number) |  | The share of ships powered by methanol in a given year. |
    | `shipping_oil_share` | dict (str -> number) |  | The share of ships powered by oil in a given year. |
    | `shipping_methanol_efficiency` | number | `0.46` | The efficiency of methanol-powered ships in the conversion of methanol to meet shipping needs (propulsion). The efficiency increase from oil can be 10-15% higher according to the [IEA ](https://www.iea-amf.org/app/webroot/files/file/Annex%20Reports/AMF_Annex_56.pdf). |
    | `shipping_oil_efficiency` | number | `0.4` | The efficiency of oil-powered ships in the conversion of oil to meet shipping needs (propulsion). Base value derived from 2011. |
    | `aviation_demand_factor` | number | `1.0` | The proportion of demand for aviation compared to today's consumption. |
    | `HVC_demand_factor` | number | `1.0` | The proportion of demand for high-value chemicals compared to today's consumption. |
    | `time_dep_hp_cop` | boolean | `true` | Consider the time dependent coefficient of performance (COP) of the heat pump. |
    | `heat_pump_sink_T_individual_heating` | number | `55.0` | The temperature heat sink used in heat pumps based on DTU / large area radiators. The value is conservatively high to cover hot water and space heating in poorly-insulated buildings. |
    | `reduce_space_heat_exogenously` | boolean | `true` | Influence on space heating demand by a certain factor (applied before losses in district heating). |
    | `reduce_space_heat_exogenously_factor` | dict (str -> number) |  | A positive factor can mean renovation or demolition of a building. If the factor is negative, it can mean an increase in floor area, increased thermal comfort, population growth. The default factors are determined by the [Eurocalc Homes and buildings decarbonization scenario ](http://tool.european-calculator.eu/app/buildings/building-types-area/?levers=1ddd4444421213bdbbbddd44444ffffff11f411111221111211l212221). |
    | `retrofitting` | any |  | Configuration for `sector.retrofitting` settings. |
    |   `retro_endogen` | boolean | `false` | Add retrofitting as an endogenous system which co-optimise space heat savings. |
    |   `cost_factor` | number | `1.0` | Weight costs for building renovation. |
    |   `interest_rate` | number | `0.04` | The interest rate for investment in building components. |
    |   `annualise_cost` | boolean | `true` | Annualise the investment costs of retrofitting. |
    |   `tax_weighting` | boolean | `false` | Weight the costs of retrofitting depending on taxes in countries. |
    |   `construction_index` | boolean | `true` | Weight the costs of retrofitting depending on labour/material costs per country. |
    | `tes` | boolean | `true` | Add option for storing thermal energy in large water pits associated with district heating systems and individual thermal energy storage (TES). |
    | `boilers` | boolean | `true` | Add option for transforming gas into heat using gas boilers. |
    | `resistive_heaters` | boolean | `true` | Add option for transforming electricity into heat using resistive heaters (independently from gas boilers). |
    | `oil_boilers` | boolean | `false` | Add option for transforming oil into heat using boilers. |
    | `biomass_boiler` | boolean | `true` | Add option for transforming biomass into heat using boilers. |
    | `overdimension_heat_generators` | dict (str -> number) |  | Add option for overdimensioning heating systems by a certain factor. This allows them to cover heat demand peaks e.g. 10% higher than those in the data with a setting of 1.1. |
    | `chp` | any |  | Configuration for `sector.chp` settings. |
    |   `enable` | boolean | `true` | Add option for using Combined Heat and Power (CHP). |
    |   `fuel` | list of string |  | Possible options are all fuels which have an existing bus and their CO2 intensity is given in the technology data. Currently possible are "gas", "oil", "methanol", "lignite", "coal" as well as "solid biomass". For all fuels except solid biomass, the techno-economic data from gas CHP is used. For the special case of solid biomass fuel, both CHP plants with and without carbon capture are added. |
    |   `micro_chp` | boolean | `false` | Add option for using gas-fired Combined Heat and Power (CHP) for decentral areas. |
    | `solar_thermal` | boolean | `true` | Add option for using solar thermal to generate heat. |
    | `solar_cf_correction` | number | `0.788457` | The correction factor for the value provided by the solar thermal profile calculations. |
    | `methanation` | boolean | `true` | Add option for transforming hydrogen and CO2 into methane using methanation. |
    | `coal_cc` | boolean | `false` | Add option for coal CHPs with carbon capture. |
    | `dac` | boolean | `true` | Add option for Direct Air Capture (DAC). |
    | `co2_vent` | boolean | `false` | Add option for vent out CO2 from storages to the atmosphere. |
    | `heat_vent` | dict (str -> boolean) |  | Heat venting by area. |
    | `marginal_cost_heat_vent` | number | `0.02` | The marginal cost of heat-venting in all heating systems. |
    | `allam_cycle_gas` | boolean | `false` | Add option to include [Allam cycle gas power plants ](https://en.wikipedia.org/wiki/Allam_power_cycle). |
    | `hydrogen_fuel_cell` | boolean | `true` | Add option to include hydrogen fuel cell for re-electrification. Assuming OCGT technology costs. |
    | `hydrogen_turbine` | boolean | `true` | Add option to include hydrogen turbine for re-electrification. Assuming OCGT technology costs. |
    | `SMR` | boolean | `true` | Add option for transforming natural gas into hydrogen and CO2 using Steam Methane Reforming (SMR). |
    | `SMR_cc` | boolean | `true` | Add option for transforming natural gas into hydrogen and CO2 using Steam Methane Reforming (SMR) and Carbon Capture (CC). |
    | `regional_oil_demand` | boolean | `true` | Spatially resolve oil demand. Set to true if regional CO2 constraints needed. |
    | `regional_coal_demand` | boolean | `false` | Regional coal demand. |
    | `regional_co2_sequestration_potential` | object |  | Add option for regionally-resolved geological carbon dioxide sequestration potentials based on [CO2StoP ](https://setis.ec.europa.eu/european-co2-storage-database_en). |
    | `co2_sequestration_potential` | dict (str -> number) |  | The potential of sequestering CO2 in Europe per year and investment period. |
    | `co2_sequestration_cost` | number | `30` | The cost of sequestering a ton of CO2 (currency/tCO2). |
    | `co2_sequestration_lifetime` | integer | `50` | The lifetime of a CO2 sequestration site (years). |
    | `co2_spatial` | boolean | `true` | Add option to spatially resolve carrier representing stored carbon dioxide. This allows for more detailed modelling of CCUTS, e.g. regarding the capturing of industrial process emissions, usage as feedstock for electrofuels, transport of carbon dioxide, and geological sequestration sites. |
    | `co2_network` | boolean | `true` | Add option for planning a new carbon dioxide transmission network. |
    | `co2_network_cost_factor` | number | `1` | The cost factor for the capital cost of the carbon dioxide transmission network. |
    | `cc_fraction` | number | `0.9` | The default fraction of CO2 captured with post-combustion capture. |
    | `hydrogen_underground_storage` | boolean | `true` | Add options for storing hydrogen underground. Storage potential depends regionally. |
    | `hydrogen_underground_storage_locations` | list of string |  | The location where hydrogen underground storage can be located. Onshore, nearshore, offshore means it must be located more than 50 km away from the sea, within 50 km of the sea, or within the sea itself respectively. |
    | `methanol` | any |  | Configuration for `sector.methanol` settings. |
    |   `regional_methanol_demand` | boolean | `false` | Spatially resolve methanol demand. Set to true if regional CO2 constraints needed. |
    |   `methanol_reforming` | boolean | `false` | Add methanol reforming. |
    |   `methanol_reforming_cc` | boolean | `false` | Add methanol reforming with carbon capture. |
    |   `methanol_to_kerosene` | boolean | `false` | Add methanol to kerosene. |
    |   `methanol_to_power` | dict (str -> boolean) |  | Add different methanol to power technologies. |
    |   `biomass_to_methanol` | boolean | `true` | Add biomass to methanol. |
    |   `biomass_to_methanol_cc` | boolean | `false` | Add biomass to methanol with carbon capture. |
    | `ammonia` | boolean \| string | `true` | Add ammonia as a carrier. It can be either true (copperplated NH3), false (no NH3 carrier) or "regional" (regionalised NH3 without network). |
    | `min_part_load_electrolysis` | number | `0` | The minimum unit dispatch (`p_min_pu`) for electrolysis. |
    | `min_part_load_fischer_tropsch` | number | `0.5` | The minimum unit dispatch (`p_min_pu`) for the Fischer-Tropsch process. |
    | `min_part_load_methanolisation` | number | `0.3` | The minimum unit dispatch (`p_min_pu`) for the methanolisation process. |
    | `min_part_load_methanation` | number | `0.3` | Minimum part load methanation. |
    | `use_fischer_tropsch_waste_heat` | number | `0.25` | Add option for using waste heat of Fischer Tropsch in district heating networks. |
    | `use_haber_bosch_waste_heat` | number | `0.25` | Use Haber-Bosch waste heat. |
    | `use_methanolisation_waste_heat` | number | `0.25` | Use methanolisation waste heat. |
    | `use_methanation_waste_heat` | number | `0.25` | Use methanation waste heat. |
    | `use_fuel_cell_waste_heat` | number | `1` | Add option for using waste heat of fuel cells in district heating networks. |
    | `use_electrolysis_waste_heat` | number | `0.25` | Add option for using waste heat of electrolysis in district heating networks. |
    | `electricity_transmission_grid` | boolean | `true` | Switch for enabling/disabling the electricity transmission grid. |
    | `electricity_distribution_grid` | boolean | `true` | Add a simplified representation of the exchange capacity between transmission and distribution grid level through a link. |
    | `electricity_distribution_grid_cost_factor` | number | `1.0` | Multiplies the investment cost of the electricity distribution grid. |
    | `electricity_grid_connection` | boolean | `true` | Add the cost of electricity grid connection for onshore wind and solar. |
    | `transmission_efficiency` | any |  | Configuration for `sector.transmission_efficiency` settings. |
    |   `enable` | list of string |  | Switch to select the carriers for which transmission efficiency is to be added. Carriers not listed assume lossless transmission. |
    |   `DC` | dict (str -> number) |  | DC transmission efficiency. |
    |   `H2 pipeline` | dict (str -> number) |  | H2 pipeline transmission efficiency. |
    |   `gas pipeline` | dict (str -> number) |  | Gas pipeline transmission efficiency. |
    |   `electricity distribution grid` | dict (str -> number) |  | Electricity distribution grid efficiency. |
    | `H2_network` | boolean | `true` | Add option for new hydrogen pipelines. |
    | `gas_network` | boolean | `true` | Add existing natural gas infrastructure, incl. LNG terminals, production and entry-points. The existing gas network is added with a lossless transport model. A length-weighted [k-edge augmentation algorithm ](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation.html#networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation) can be run to add new candidate gas pipelines such that all regions of the model can be connected to the gas network. When activated, all the gas demands are regionally disaggregated as well. |
    | `H2_retrofit` | boolean | `false` | Add option for retrofiting existing pipelines to transport hydrogen. |
    | `H2_retrofit_capacity_per_CH4` | number | `0.6` | The ratio for H2 capacity per original CH4 capacity of retrofitted pipelines. The [European Hydrogen Backbone (April, 2020) p.15 ](https://gasforclimate2050.eu/wp-content/uploads/2020/07/2020_European-Hydrogen-Backbone_Report.pdf) 60% of original natural gas capacity could be used in cost-optimal case as H2 capacity. |
    | `gas_network_connectivity_upgrade` | number | `1` | The number of desired edge connectivity (k) in the length-weighted [k-edge augmentation algorithm ](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation.html#networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation) used for the gas network. |
    | `gas_distribution_grid` | boolean | `true` | Add a gas distribution grid. |
    | `gas_distribution_grid_cost_factor` | number | `1.0` | Multiplier for the investment cost of the gas distribution grid. |
    | `biomass_spatial` | boolean | `true` | Add option for resolving biomass demand regionally. |
    | `biomass_transport` | boolean | `false` | Add option for transporting solid biomass between nodes. |
    | `biogas_upgrading` | boolean | `true` | Biogas upgrading. |
    | `biogas_upgrading_cc` | boolean | `false` | Add option to capture CO2 from biomass upgrading. |
    | `conventional_generation` | dict (str -> string) |  | Add a more detailed description of conventional carriers. Any power generation requires the consumption of fuel from nodes representing that fuel. |
    | `biomass_to_liquid` | boolean | `true` | Add option for transforming solid biomass into liquid fuel with the same properties as oil. |
    | `biomass_to_liquid_cc` | boolean | `false` | Add option for transforming solid biomass into liquid fuel with the same properties as oil with carbon capture. |
    | `electrobiofuels` | boolean | `true` | Electrobiofuels. |
    | `biosng` | boolean | `false` | Add option for transforming solid biomass into synthesis gas with the same properties as natural gas. |
    | `biosng_cc` | boolean | `false` | Add option for transforming solid biomass into synthesis gas with the same properties as natural gas with carbon capture. |
    | `bioH2` | boolean | `false` | Add option for transforming solid biomass into hydrogen with carbon capture. |
    | `municipal_solid_waste` | boolean | `false` | Add option for municipal solid waste. |
    | `limit_max_growth` | any |  | Configuration for `sector.limit_max_growth` settings. |
    |   `enable` | boolean | `false` | Add option to limit the maximum growth of a carrier. |
    |   `factor` | number | `1.3` | The maximum growth factor of a carrier (e.g. 1.3 allows  30% larger than max historic growth). |
    |   `max_growth` | dict (str -> number) |  | The historic maximum growth of a carrier. |
    |   `max_relative_growth` | dict (str -> number) |  | The historic maximum relative growth of a carrier. |
    | `enhanced_geothermal` | any |  | Configuration for `sector.enhanced_geothermal` settings. |
    |   `enable` | boolean | `false` | Add option to include Enhanced Geothermal Systems. |
    |   `flexible` | boolean | `true` | Add option for flexible operation (see Ricks et al. 2024). |
    |   `max_hours` | integer | `240` | The maximum hours the reservoir can be charged under flexible operation. |
    |   `max_boost` | number | `0.25` | The maximum boost in power output under flexible operation. |
    |   `var_cf` | boolean | `true` | Add option for variable capacity factor (see Ricks et al. 2024). |
    |   `sustainability_factor` | number | `0.0025` | Share of sourced heat that is replenished by the earth's core (see details in [build_egs_potentials.py ](https://github.com/PyPSA/pypsa-eur-sec/blob/master/scripts/build_egs_potentials.py)). |
    | `solid_biomass_import` | any |  | Configuration for `sector.solid_biomass_import` settings. |
    |   `enable` | boolean | `false` | Add option to include solid biomass imports. |
    |   `price` | number | `54` | Price for importing solid biomass (currency/MWh). |
    |   `max_amount` | number | `1390` | Maximum solid biomass import potential (TWh). |
    |   `upstream_emissions_factor` | number | `0.1` | Upstream emissions of solid biomass imports. |
    | `imports` | any |  | Configuration for `sector.imports` settings. |
    |   `enable` | boolean | `false` | Add option to include renewable energy imports. |
    |   `limit` | number |  | Maximum allowed renewable energy imports (TWh). |
    |   `limit_sense` | string | `<=` | Sense of the limit. |
    |   `price` | dict (str -> number) |  | Price for importing renewable energy of carrier. |
    

**YAML Syntax**

```yaml
{{ yaml_section("sector") }}
```


!!! note
    Only used for sector-coupling studies.


## `industry` {#industry_cf}

Only used for sector-coupling studies.

Configuration for `industry` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `St_primary_fraction` | dict (str -> number) |  | The fraction of steel produced via primary route versus secondary route (scrap+EAF). Current fraction is 0.6. |
| `DRI_fraction` | dict (str -> number) |  | The fraction of the primary route DRI + EAF. |
| `H2_DRI` | number | `1.7` | The hydrogen consumption in Direct Reduced Iron (DRI) Mwh_H2 LHV/ton_Steel from 51kgH2/tSt in [Vogl et al (2018) ](https://doi.org/10.1016/j.jclepro.2018.08.279). |
| `elec_DRI` | number | `0.322` | The electricity consumed in Direct Reduced Iron (DRI) shaft. From [HYBRIT brochure ](https://ssabwebsitecdn.azureedge.net/-/media/hybrit/files/hybrit_brochure.pdf). |
| `Al_primary_fraction` | dict (str -> number) |  | The fraction of aluminium produced via the primary route versus scrap. Current fraction is 0.4. |
| `MWh_NH3_per_tNH3` | number | `5.166` | The energy amount per ton of ammonia (LHV). |
| `MWh_CH4_per_tNH3_SMR` | number | `10.8` | The energy amount of methane needed to produce a ton of ammonia using steam methane reforming (SMR). Value derived from 2012's demand from [Center for European Policy Studies (2008) ](https://ec.europa.eu/docsroom/documents/4165/attachments/1/translations/en/renditions/pdf). |
| `MWh_elec_per_tNH3_SMR` | number | `0.7` | The energy amount of electricity needed to produce a ton of ammonia using steam methane reforming (SMR). same source, assuming 94-6% split methane-elec of total energy demand 11.5 MWh/tNH3. |
| `MWh_H2_per_tNH3_electrolysis` | number | `5.93` | The energy amount of hydrogen needed to produce a ton of ammonia using Haber–Bosch process. From [Wang et al (2018) ](https://doi.org/10.1016/j.joule.2018.04.017), Base value assumed around 0.197 tH2/tHN3 (>3/17 since some H2 lost and used for energy). |
| `MWh_elec_per_tNH3_electrolysis` | number | `0.2473` | The energy amount of electricity needed to produce a ton of ammonia using Haber–Bosch process. From [Wang et al (2018) ](https://doi.org/10.1016/j.joule.2018.04.017), Table 13 (air separation and HB). |
| `MWh_NH3_per_MWh_H2_cracker` | number | `1.46` | The energy amount of amonia needed to produce an energy amount hydrogen using ammonia cracker. |
| `NH3_process_emissions` | number | `24.5` | The emission of ammonia production from steam methane reforming (SMR). From UNFCCC for 2015 for EU28. |
| `petrochemical_process_emissions` | number | `25.5` | The emission of petrochemical production. From UNFCCC for 2015 for EU28. |
| `HVC_primary_fraction` | dict (str -> number) |  | The fraction of high value chemicals (HVC) produced via primary route. |
| `HVC_mechanical_recycling_fraction` | dict (str -> number) |  | The fraction of high value chemicals (HVC) produced using mechanical recycling. |
| `HVC_chemical_recycling_fraction` | dict (str -> number) |  | The fraction of high value chemicals (HVC) produced using chemical recycling. |
| `HVC_environment_sequestration_fraction` | number | `0.0` | The fraction of high value chemicals (HVC) put into landfill resulting in additional carbon sequestration. The default value is 0. |
| `waste_to_energy` | boolean | `false` | Switch to enable expansion of waste to energy CHPs for conversion of plastics. Default is false. |
| `waste_to_energy_cc` | boolean | `false` | Switch to enable expansion of waste to energy CHPs for conversion of plastics with carbon capture. Default is false. |
| `sector_ratios_fraction_future` | dict (str -> number) |  | The fraction of total progress in fuel and process switching achieved in the industry sector. |
| `basic_chemicals_without_NH3_production_today` | number | `69.0` | The amount of basic chemicals produced without ammonia (= 86 Mtethylene-equiv - 17 MtNH3). |
| `HVC_production_today` | number | `52.0` | The amount of high value chemicals (HVC) produced. This includes ethylene, propylene and BTX. From [DECHEMA (2017) ](https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf), Figure 16, page 107. |
| `MWh_elec_per_tHVC_mechanical_recycling` | number | `0.547` | The energy amount of electricity needed to produce a ton of high value chemical (HVC) using mechanical recycling. From SI of [Meys et al (2020) ](https://doi.org/10.1016/j.resconrec.2020.105010), Table S5, for HDPE, PP, PS, PET. LDPE would be 0.756. |
| `MWh_elec_per_tHVC_chemical_recycling` | number | `6.9` | The energy amount of electricity needed to produce a ton of high value chemical (HVC) using chemical recycling. The default value is based on pyrolysis and electric steam cracking. From [Material Economics (2019) ](https://materialeconomics.com/latest-updates/industrial-transformation-2050), page 125. |
| `chlorine_production_today` | number | `9.58` | The amount of chlorine produced. From [DECHEMA (2017) ](https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf), Table 7, page 43. |
| `MWh_elec_per_tCl` | number | `3.6` | The energy amount of electricity needed to produce a ton of chlorine. From [DECHEMA (2017) ](https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf), Table 6 page 43. |
| `MWh_H2_per_tCl` | number | `-0.9372` | The energy amount of hydrogen needed to produce a ton of chlorine. The value is negative since hydrogen produced in chloralkali process. From [DECHEMA (2017) ](https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf), page 43. |
| `methanol_production_today` | number | `1.5` | The amount of methanol produced. From [DECHEMA (2017) ](https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf), page 62. |
| `MWh_elec_per_tMeOH` | number | `0.167` | The energy amount of electricity needed to produce a ton of methanol from fossil gas. From [DECHEMA (2017) ](https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf), Table 14, page 65. |
| `MWh_CH4_per_tMeOH` | number | `10.25` | The energy amount of methane needed to produce a ton of methanol from fossil gas. From [DECHEMA (2017) ](https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf), Table 14, page 65. |
| `MWh_MeOH_per_tMeOH` | number | `5.528` | The energy amount per ton of methanol (LHV). From [DECHEMA (2017) ](https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf), page 74. |
| `hotmaps_locate_missing` | boolean | `false` | Locate industrial sites without valid locations based on city and countries. |
| `reference_year` | integer | `2023` | The year used as the baseline for industrial energy demand and production. Data extracted from [JRC-IDEES 2015 ](https://data.jrc.ec.europa.eu/dataset/jrc-10110-10001). |
| `oil_refining_emissions` | number | `0.013` | The emissions from oil fuel processing (e.g. oil in petrochemical refinieries). The default value of 0.013 tCO2/MWh is based on DE statistics for 2019; the EU value is very similar. |

**YAML Syntax**

```yaml
{{ yaml_section("industry") }}
```

!!! note
    Only used for sector-coupling studies.


## `costs` {#costs_cf}

Configuration for `costs` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `year` | integer | `2050` | Year for which to retrieve cost assumptions of `data/costs/primary/<version>/costs_<year>.csv`. |
| `social_discountrate` | number | `0.02` | Social discount rate to compare costs in different investment periods. 0.02 corresponds to a social discount rate of 2%. |
| `fill_values` | any |  | Configuration for `costs.fill_values` settings. |
| ↳ `FOM` | number | `0` | Default fixed operation and maintenance cost. |
| ↳ `VOM` | number | `0` | Default variable operation and maintenance cost. |
| ↳ `efficiency` | number | `1` | Default efficiency. |
| ↳ `fuel` | number | `0` | Default fuel cost. |
| ↳ `investment` | number | `0` | Default investment cost. |
| ↳ `lifetime` | integer | `25` | Default lifetime in years. |
| ↳ `CO2 intensity` | number | `0` | Default CO2 intensity. |
| ↳ `discount rate` | number | `0.07` | Default discount rate. |
| ↳ `standing losses` | number | `0` | Default standing losses. |
| `custom_cost_fn` | string \| null | `data/custom_costs.csv` | Path to the custom costs file. None if it should not be used. Default `data/custom_costs.csv` contains minor adjustments for stabilising the optimisation results. |
| `overwrites` | dict (str -> dict (str -> number)) |  | For the given parameters and technologies, assumptions about their parameter are overwritten the corresponding value of the technology. |
| `capital_cost` | dict (str -> number) |  | For the given technologies, assumptions about their capital investment costs are set to the corresponding value. Optional; overwrites cost assumptions from `resources/costs.csv`. |
| `marginal_cost` | dict (str -> number) |  | For the given technologies, assumptions about their marginal operating costs are set to the corresponding value. Optional; overwrites cost assumptions from `resources/costs.csv`. |
| `emission_prices` | any |  | Configuration for `costs.emission_prices` settings. |
| ↳ `enable` | boolean | `false` | Add cost for a carbon-dioxide price configured in `costs: emission_prices: co2` to `marginal_cost` of generators. Config setting can also be enabled with the keyword `Ep` in the `{opts}` wildcard for electricity-only runs. |
| ↳ `co2` | number \| dict (str -> number) | `0.0` | Exogenous price of carbon-dioxide. In electricity-only runs it is added to the marginal costs of fossil-fuelled generators according to their carbon intensity, while for sector networks it applies to emissions ending up in CO2 atmosphere. |
| ↳ `dynamic` | boolean | `false` | Add time-varying cost for a carbon-dioxide price based on historical values built by the rule `build_co2_prices`. |
| ↳ `rolling_window` | integer | `90` | Rolling window (in days) for smoothing the historical CO2 prices when `dynamic` is set to True. |

**YAML Syntax**

```yaml
{{ yaml_section("costs") }}
```


## `clustering` {#clustering_cf}

use `min` in `p_nom_max:` for more conservative assumptions.

Configuration for `clustering` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `mode` | enum (`busmap`, `custom_busmap`, `administrative`, `custom_busshapes`) | `busmap` | 'busmap': Default. 'custom_busmap': Enable the use of custom busmaps in rule `cluster_network`. If activated the rule looks for provided busmaps at `data/busmaps/base_s_{clusters}_{base_network}.csv` which should have the same format as `resources/busmap_base_s_{clusters}.csv`, i.e. the index should contain the buses of `networks/base_s.nc`. {base_network} is the name of the selected base_network in electricity, e.g. `gridkit`, `osm-prebuilt`, or `osm-raw`. 'administrative': Clusters and indexes the network based on the administrative regions of the countries based on `nuts3_shapes.geojson` (level: 1, 2, 3, bz). To activate this, additionally set the `clusters` wildcard in `scenario` to 'adm'. 'custom_busshapes': Enable the use of custom shapes in rule `cluster_network`. If activated the rule looks for provided busshapes at `data/busshapes/base_s_{clusters}_{base_network}.geojson`. |
| `administrative` | any |  | Configuration for `clustering.administrative` settings. |
| ↳ `level` | enum (`0`, `1`, `2`, `3`, `bz`) | `1` | Level of administrative regions to cluster the network. 0: Country level, 1: NUTS1 level, 2: NUTS2 level, 3: NUTS3 level, 'bz': Bidding zones. Only applies when mode is set to `administrative`. Note that non-NUTS countries 'BA', 'MD', 'UA', and 'XK' can only be clustered to level 0 and 1. |
| ↳ `countries` | dict (str -> integer) |  | Optionally include dictionary of individual country codes and their individual NUTS levels. Overwrites country-specific `level`. For example: `{'DE': 1, 'FR': 2}`. Only applies when mode is set to `administrative`. |
| `focus_weights` | boolean \| dict (str -> number) | `false` | Optionally specify the focus weights for the clustering of countries. For instance: `DE: 0.8` will distribute 80% of all nodes to Germany and 20% to the rest of the countries. Only applies when mode is set to `busmap`. |
| `copperplate_regions` | list of list of string |  | Optionally specify the regions to copperplate as a list of groups. Each group is a list of region codes that will be connected with infinite capacity lines. |
| `build_bidding_zones` | any |  | Configuration for `clustering.build_bidding_zones` settings. |
| ↳ `remove_islands` | boolean | `false` | Exclude from the shape file the Balearic Islands, Bornholm, the Canary Islands, the Orkney Islands, the Shetland Islands, the Azores Islands and Madeira. |
| ↳ `aggregate_to_tyndp` | boolean | `false` | Adjust the shape file to the TYNDP topology. Aggregate the Southern Norwegian bidding zones and extract Crete as a separate zone from the Greek shape. |
| `simplify_network` | any |  | Configuration for `clustering.simplify_network` settings. |
| ↳ `to_substations` | boolean | `false` | Aggregates all nodes without power injection (positive or negative, i.e. demand or generation) to electrically closest ones. |
| ↳ `remove_stubs` | boolean | `true` | Controls whether radial parts of the network should be recursively aggregated. Defaults to true. |
| ↳ `remove_stubs_across_borders` | boolean | `false` | Controls whether radial parts of the network should be recursively aggregated across borders. Defaults to true. |
| `cluster_network` | any |  | Configuration for `clustering.cluster_network` settings. |
| ↳ `algorithm` | enum (`kmeans`, `hac`) | `kmeans` | Clustering algorithm to use. |
| ↳ `hac_features` | list of string |  | List of meteorological variables contained in the weather data cutout that should be considered for hierarchical clustering. |
| `exclude_carriers` | list of string |  | List of carriers which will not be aggregated. If empty, all carriers will be aggregated. |
| `consider_efficiency_classes` | boolean \| list of number | `false` | Aggregate each carrier into efficiency classes defined by quantile boundaries. If True, uses [0.1, 0.9] as default quantiles (labels: Q0, Q10, Q90). If a list of floats, defines custom quantile boundaries, e.g. [0.1, 0.5, 0.9]. |
| `aggregation_strategies` | any |  | Configuration for `clustering.aggregation_strategies` settings. |
| ↳ `generators` | dict (str -> string) |  | Aggregates the component according to the given strategy. For example, if sum, then all values within each cluster are summed to represent the new generator. |
| ↳ `buses` | dict (str -> string) |  | Aggregates the component according to the given strategy. For example, if sum, then all values within each cluster are summed to represent the new bus. |
| `temporal` | any |  | Configuration for `clustering.temporal` settings. |
| ↳ `resolution_elec` | boolean \| string | `false` | Resample the time-resolution by averaging over every `n` snapshots in `prepare_network`. **Warning:** This option should currently only be used with electricity-only networks, not for sector-coupled networks. |
| ↳ `resolution_sector` | boolean \| string | `false` | Resample the time-resolution by averaging over every `n` snapshots in `prepare_sector_network`. |

**YAML Syntax**

```yaml
{{ yaml_section("clustering") }}
```

!!! tip
    use `min` in `p_nom_max:` for more conservative assumptions.


## `adjustments` {#adjustments_cf}

Configuration for top-level adjustments key.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `electricity` | boolean \| any | `false` | Parameter adjustments applied in `prepare_network`. |
| ↳ `factor` | boolean \| dict (str -> dict (str -> dict (str -> number \| dict (str -> any)))) | `false` | Multiply original value with given factor |
| ↳ `absolute` | boolean \| dict (str -> dict (str -> dict (str -> number \| dict (str -> any)))) | `false` | Set attribute to absolute value. Can be also a dictionary with planning horizons as keys. |
| `sector` | boolean \| any |  | Parameter adjustments applied in `prepare_sector_network`. |
| ↳ `factor` | boolean \| dict (str -> dict (str -> dict (str -> number \| dict (str -> any)))) | `false` | Multiply original value with given factor |
| ↳ `absolute` | boolean \| dict (str -> dict (str -> dict (str -> number \| dict (str -> any)))) | `false` | Set attribute to absolute value. Can be also a dictionary with planning horizons as keys. |

**YAML Syntax**

```yaml
{{ yaml_section("adjustments") }}
```


## `solving` {#solving_cf}

Configuration for `solving` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `options` | any |  | Configuration for `solving.options` settings. |
| ↳ `clip_p_max_pu` | number | `0.01` | To avoid too small values in the renewables` per-unit availability time series values below this threshold are set to zero. |
| ↳ `load_shedding` | any |  | Configuration for `solving.options.load_shedding` settings. |
| ↳↳ `enable` | boolean | `false` | Enable load shedding by adding high-cost generators to avoid infeasibilities. Requires either all_carriers: true or at least one entry in carriers. |
| ↳↳ `default_cost` | number | `100000` | The default cost for load-shedding in the unit of the bus carrier (e.g. EUR/MWh for electricity, EUR/t_CO2 for CO2). Must be positive. |
| ↳↳ `all_carriers` | boolean | `true` | Switch to apply load shedding to all carriers. Otherwise, load shedding will be applied to listed carriers only. |
| ↳↳ `carriers` | dict (str -> number) | `{}` | Dictionary of carriers and their specific load shedding cost in the unit of the bus carrier (e.g. EUR/MWh for electricity, EUR/t_CO2 for CO2). If load shedding is enabled for all carriers, the default cost is assumed for non-listed carriers. |
| ↳ `load_sinks` | any |  | Configuration for `solving.options.load_sinks` settings. |
| ↳↳ `enable` | boolean | `false` | Add load sinks by adding negative-cost, energy consuming generators to avoid infeasibilities by absorbing excess energy. Requires either all_carriers: true or at least one entry in carriers. |
| ↳↳ `default_cost` | number | `100000` | The default cost for load sinks in the unit of the bus carrier (e.g. EUR/MWh for electricity, EUR/t_CO2 for CO2). Must be positive. |
| ↳↳ `all_carriers` | boolean | `false` | Switch to add load sinks for all carriers. Otherwise, load sinks will be added for listed carriers only. |
| ↳↳ `carriers` | dict (str -> number) | `{}` | Dictionary of carriers and their specific load sink cost in the unit of the bus carrier (e.g. EUR/MWh for electricity, EUR/t_CO2 for CO2). If load sinks are added for all carriers, the default cost is assumed for non-listed carriers. |
| ↳ `curtailment_mode` | boolean | `false` | Fixes the dispatch profiles of generators with time-varying p_max_pu by setting `p_min_pu = p_max_pu` and adds an auxiliary curtailment generator (with negative sign to absorb excess power) at every AC bus. This can speed up the solving process as the curtailment decision is aggregated into a single generator per region. Defaults to `false`. |
| ↳ `noisy_costs` | boolean | `true` | Add random noise to marginal cost of generators by `\mathcal{U}(0.009,0,011)` and capital cost of lines and links by `\mathcal{U}(0.09,0,11)`. |
| ↳ `skip_iterations` | boolean | `true` | Skip iterating, do not update impedances of branches. Defaults to true. |
| ↳ `rolling_horizon` | boolean | `false` | Switch for rule `solve_operations_network` whether to optimize the network in a rolling horizon manner, where the snapshot range is split into slices of size `horizon` which are solved consecutively. This setting has currently no effect on sector-coupled networks. |
| ↳ `seed` | integer | `123` | Random seed for increased deterministic behaviour. |
| ↳ `custom_extra_functionality` | string \| null | `../data/custom_extra_functionality.py` | Path to a Python file with custom extra functionality code to be injected into the solving rules of the workflow relative to `rules` directory. |
| ↳ `io_api` | string \| null |  | Passed to linopy and determines the API used to communicate with the solver. With the `'lp'` and `'mps'` options linopy passes a file to the solver; with the `'direct'` option (only supported for HIGHS and Gurobi) linopy uses an in-memory python API resulting in better performance. |
| ↳ `track_iterations` | boolean | `false` | Flag whether to store the intermediate branch capacities and objective function values are recorded for each iteration in `network.lines['s_nom_opt_X']` (where `X` labels the iteration) |
| ↳ `min_iterations` | integer | `2` | Minimum number of solving iterations in between which resistance and reactence (`x/r`) are updated for branches according to `s_nom_opt` of the previous run. |
| ↳ `max_iterations` | integer | `3` | Maximum number of solving iterations in between which resistance and reactence (`x/r`) are updated for branches according to `s_nom_opt` of the previous run. |
| ↳ `transmission_losses` | integer | `2` | Add piecewise linear approximation of transmission losses based on n tangents. Defaults to 0, which means losses are ignored. |
| ↳ `linearized_unit_commitment` | boolean | `true` | Whether to optimise using the linearized unit commitment formulation. |
| ↳ `horizon` | integer | `365` | Number of snapshots to consider in each iteration. Defaults to 100. |
| ↳ `overlap` | integer | `0` | Number of overlapping snapshots between consecutive iterations in rolling horizon optimization. Defaults to 0, which means no overlap. |
| ↳ `post_discretization` | any |  | Configuration for `solving.options.post_discretization` settings. |
| ↳↳ `enable` | boolean | `false` | Switch to enable post-discretization of the network. Disabled by default. |
| ↳↳ `line_unit_size` | number | `1700` | Discrete unit size of lines in MW. |
| ↳↳ `line_threshold` | number | `0.3` | The threshold relative to the discrete line unit size beyond which to round up to the next unit. |
| ↳↳ `link_unit_size` | dict (str -> number) |  | Discrete unit size of links in MW by carrier (given in dictionary style). |
| ↳↳ `link_threshold` | dict (str -> number) |  | The threshold relative to the discrete link unit size beyond which to round up to the next unit by carrier (given in dictionary style). |
| ↳↳ `fractional_last_unit_size` | boolean | `false` | When true, links and lines can be built up to p_nom_max. When false, they can only be built up to a multiple of the unit size. |
| ↳ `keep_files` | boolean | `false` | Whether to keep LPs and MPS files after solving. |
| ↳ `store_model` | boolean | `false` | Store the linopy model to a NetCDF file after solving. Not supported with rolling_horizon. Not scenario-aware. |
| ↳ `model_kwargs` | any |  | Configuration for `solving.options.model_kwargs` settings. |
| ↳↳ `solver_dir` | string | `""` | Absolute path to the directory where linopy saves files. |
| `agg_p_nom_limits` | any |  | Configuration for `solving.agg_p_nom_limits` settings. |
| ↳ `agg_offwind` | boolean | `false` | Aggregate together all the types of offwind when writing the constraint (`offwind-all` as a carrier in the `.csv` file). Default is false. |
| ↳ `agg_solar` | boolean | `false` | Aggregate together all the types of electric solar when writing the constraint (`solar-all` as a carrier in the `.csv` file). Default is false. |
| ↳ `include_existing` | boolean | `false` | Take existing capacities into account when writing the constraint. Default is false. |
| ↳ `file` | string | `data/agg_p_nom_minmax.csv` | Reference to `.csv` file specifying per carrier generator nominal capacity constraints for individual countries and planning horizons. Defaults to `data/agg_p_nom_minmax.csv`. |
| `constraints` | any |  | Configuration for `solving.constraints` settings. |
| ↳ `CCL` | boolean | `false` | Add minimum and maximum levels of generator nominal capacity per carrier for individual countries. These can be specified in the file linked at `electricity: agg_p_nom_limits` in the configuration. File defaults to `data/agg_p_nom_minmax.csv`. Does not work with a time resolution resampling. |
| ↳ `EQ` | boolean \| string | `false` | Require each country or node to on average produce a minimal share of its total consumption itself. Example: `EQ0.5c` demands each country to produce on average at least 50% of its consumption; `EQ0.5` demands each node to produce on average at least 50% of its consumption. |
| ↳ `BAU` | boolean | `false` | Add a per-`carrier` minimal overall capacity; i.e. at least `40GW` of `OCGT` in Europe; configured in `electricity: BAU_mincapacities` |
| ↳ `SAFE` | boolean | `false` | Add a capacity reserve margin of a certain fraction above the peak demand to which renewable generators and storage do *not* contribute. Ignores network. |
| `solver` | any |  | Configuration for `solving.solver` settings. |
| ↳ `name` | string | `gurobi` | Solver to use for optimisation problems in the workflow; e.g. clustering and linear optimal power flow. |
| ↳ `options` | string | `gurobi-default` | Link to specific parameter settings. |
| `solver_options` | dict (str -> object) |  | Dictionaries with solver-specific parameter settings. |
| `check_objective` | any |  | Configuration for `solving.check_objective` settings. |
| ↳ `enable` | boolean | `false` | Enable objective value checking. |
| ↳ `expected_value` | number \| null |  | Expected objective value. |
| ↳ `atol` | number | `1000000` | Absolute tolerance. |
| ↳ `rtol` | number | `0.01` | Relative tolerance. |
| `oetc` | any \| null |  | Configuration options for Open Energy Transition Computing (OETC) cluster support. |
| ↳ `name` | string | `pypsa-eur` | Name identifier for the OETC job. |
| ↳ `authentication_server_url` | string | `""` | URL of the OETC authentication server for job submission. |
| ↳ `orchestrator_server_url` | string | `""` | URL of the OETC orchestrator server for job management. |
| ↳ `cpu_cores` | integer | `8` | Number of CPU cores to request for the OETC job. (includes RAM amount at the moment with a factor of 8) |
| ↳ `disk_space_gb` | integer | `50` | Amount of disk space in gigabytes to request for the OETC job. |
| ↳ `delete_worker_on_error` | boolean | `true` | Whether to delete the worker instance when an error occurs during job execution. |
| `mem_mb` | integer | `128000` | Estimated maximum memory requirement for solving networks (MB). |
| `memory_logging_frequency` | integer | `5` | Interval in seconds at which memory usage is logged. |
| `runtime` | string | `48h` | Runtime in humanfriendly style. |

**YAML Syntax**

```yaml
{{ yaml_section("solving") }}
```


## `data` {#data_cf}

Controls which versions of input data are used for building the model.
Versions that are available for each dataset can be found in `data/versions.csv`.
By default, we retrieve the `latest` supported version for each dataset from an archive source.
This means that when upgrading between PyPSA-Eur versions, new versions of input data may also be downloaded and used.
To freeze a model to a specific version of input data, you can set a specific version in the `version` field for each dataset to one specific version as listed in `data/versions.csv`.

Some datasets support `primary` or `build` as a source option, meaning that the data can be retrieved from the original
data source or build it from the latest available data.
See the `data/versions.csv` file for all available datasets and their sources/versions that are supported.

??? note "Details"

    Configuration for `data` settings.
    
    | Property | Type | Default | Description |
    |----------|------|---------|-------------|
    | `hotmaps_industrial_sites` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `enspreso_biomass` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `osm` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `worldbank_urban_population` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `worldbank_commodity_prices` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `gem_europe_gas_tracker` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `gem_gcct` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `instrat_co2_prices` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `co2stop` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `nitrogen_statistics` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `eu_nuts2013` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `eu_nuts2021` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `eurostat_balances` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `eurostat_household_balances` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `wdpa` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `wdpa_marine` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `luisa_land_cover` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `jrc_idees` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `scigrid_gas` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `seawater_temperature` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `swiss_energy_balances` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `synthetic_electricity_demand` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `opsd_electricity_demand` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `entsoe_electricity_demand` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `neso_electricity_demand` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `copernicus_land_cover` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `ship_raster` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `eez` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `nuts3_population` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `gdp_per_capita` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `population_count` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `ghg_emissions` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `gebco` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `attributed_ports` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `corine` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `emobility` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `h2_salt_caverns` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `lau_regions` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `aquifer_data` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `osm_boundaries` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `gem_gspt` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `tyndp` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `powerplants` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `costs` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `country_runoff` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `country_hdd` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `natura` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `bfs_road_vehicle_stock` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `bfs_gdp_and_population` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `mobility_profiles` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `cutout` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `dh_areas` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `geothermal_heat_utilisation_potentials` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `jrc_ardeco` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `jrc_energy_atlas` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `desnz_electricity_consumption` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `ons_lad` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `bidding_zones_electricitymaps` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    | `bidding_zones_entsoepy` | any |  | Configuration for a single data source. |
    |   `source` | enum (`archive`, `primary`, `build`) | `archive` | Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source. |
    |   `version` | string | `latest` | Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source. |
    

**YAML Syntax**

```yaml
{{ yaml_section("data") }}
```


## `overpass_api` {#overpass_api_cf}

Configuration for `overpass_api` settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `url` | string | `https://overpass-api.de/api/interpreter` | Overpass API endpoint URL. See [Overpass API Wiki ](https://wiki.openstreetmap.org/wiki/Overpass_API#Public_Overpass_API_instances) for available public instances. |
| `max_tries` | integer | `5` | Maximum retry attempts for Overpass API requests. Please be respectful to the Overpass API fair use policy of the individual instances. |
| `timeout` | integer | `600` | Timeout in seconds for Overpass API requests. |
| `user_agent` | any |  | Configuration for `overpass_api.user_agent` settings. |
| ↳ `project_name` | string | `PyPSA-Eur` | Project name used to identify the user agent of the Overpass API requests. |
| ↳ `email` | string | `contact@pypsa.org` | Contact email address for the project using the Overpass API. |
| ↳ `website` | string | `https://github.com/PyPSA/pypsa-eur` | Website URL for the project using the Overpass API. |

**YAML Syntax**

```yaml
{{ yaml_section("overpass_api") }}
```



## `plotting` {#plotting_cf}

```yaml
{{ yaml_section("plotting", source="plotting") }}
```

