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

{{ schema_table("logging") }}

**YAML Syntax**

```yaml
{{ yaml_section("logging") }}
```


## `remote` {#remote_cf}

"Remote" indicates the address of a server used for data exchange, often for clusters and data pushing/pulling.

Configuration for top level `remote` settings.

{{ schema_table("remote") }}

**YAML Syntax**

```yaml
{{ yaml_section("remote") }}
```


## `run` {#run_cf}

It is common conduct to analyse energy system optimisation models for **multiple scenarios** for a variety of reasons,
e.g. assessing their sensitivity towards changing the temporal and/or geographical resolution or investigating how
investment changes as more ambitious greenhouse-gas emission reduction targets are applied.

The `run` section is used for running and storing scenarios with different configurations which are not covered by [wildcards](#wildcards).
It determines the path at which resources, networks and results are stored.
Therefore the user can run different configurations within the same directory.

Configuration for top level `run` settings.

{{ schema_table("run") }}

**YAML Syntax**

```yaml
{{ yaml_section("run") }}
```


## `foresight` {#foresight_cf}

[planning_horizons](#planning_horizons) in scenario has to be set.

Configuration for `foresight` settings.

- **Type:** enum (`overnight`, `myopic`, `perfect`)
- **Default:** `overnight`

**YAML Syntax**

```yaml
{{ yaml_section("foresight") }}
```

!!! note
    If you use myopic or perfect foresight, the planning horizon in
    [planning_horizons](#planning_horizons) in scenario has to be set.


## `scenario` {#scenario}

The `scenario` section is an extraordinary section of the config file
that is strongly connected to the [wildcards](#wildcards) and is designed to
facilitate running multiple scenarios through a single command


```console
# for electricity-only studies
   $ snakemake -call solve_elec_networks

   # for sector-coupling studies
   $ snakemake -call solve_sector_networks

For each wildcard, a **list of values** is provided. The rule
```

`solve_all_elec_networks` will trigger the rules for creating
`results/networks/base_s_{clusters}_elec_{opts}.nc` for **all
combinations** of the provided wildcard values as defined by Python's
[itertools.product(...)
](https://docs.python.org/2/library/itertools.html#itertools.product) function
that snakemake's [expand(...) function
](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#targets)
uses.

An exemplary dependency graph (starting from the simplification rules) then looks like this:

Configuration for top level `scenario` settings.

{{ schema_table("scenario") }}

**YAML Syntax**

```yaml
{{ yaml_section("scenario") }}
```


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

{{ schema_table("snapshots") }}

**YAML Syntax**

```yaml
{{ yaml_section("snapshots") }}
```


## `enable` {#enable_cf}

Switches for some rules and optional features.

Configuration for `enable` settings.

{{ schema_table("enable") }}

**YAML Syntax**

```yaml
{{ yaml_section("enable") }}
```


## `co2 budget` {#CO2_budget_cf}

sector_opts.

Configuration for `co2_budget` settings.

- **Type:** dict (str -> number)

**YAML Syntax**

```yaml
{{ yaml_section("co2_budget") }}
```

!!! note
    this parameter is over-ridden if `Co2Lx` or `cb` is set in
    sector_opts.


## `electricity` {#electricity_cf}

Configuration for `electricity` settings.

{{ schema_table("electricity") }}

**YAML Syntax**

```yaml
{{ yaml_section("electricity") }}
```


## `atlite` {#atlite_cf}

Define and specify the `atlite.Cutout` used for calculating renewable potentials and time-series. All options except for `features` are directly used as [cutout parameters ](https://atlite.readthedocs.io/en/latest/ref_api.html#cutout).

Configuration for `atlite` settings.

{{ schema_table("atlite") }}

**YAML Syntax**

```yaml
{{ yaml_section("atlite") }}
```


## `renewable` {#renewable_cf}

### `onwind`

Configuration for onshore wind.

{{ schema_table("renewable.onwind") }}

Configuration for offshore wind.

{{ schema_table("renewable.offwind-ac") }}

Configuration for offshore wind.

{{ schema_table("renewable.offwind-dc") }}

Configuration for offshore wind.

{{ schema_table("renewable.offwind-float") }}

Configuration for solar PV.

{{ schema_table("renewable.solar") }}

Configuration for hydropower.

{{ schema_table("renewable.hydro") }}

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

{{ schema_table("conventional") }}

**YAML Syntax**

```yaml
{{ yaml_section("conventional") }}
```


## `lines` {#lines_cf}

Configuration for `lines` settings.

{{ schema_table("lines") }}

**YAML Syntax**

```yaml
{{ yaml_section("lines") }}
```


## `links` {#links_cf}

Configuration for `links` settings.

{{ schema_table("links") }}

**YAML Syntax**

```yaml
{{ yaml_section("links") }}
```


## `transmission_projects` {#transmission_projects_cf}

Allows to define additional transmission projects that will be added to the base network, e.g., from the TYNDP 2020 dataset. The projects are read in from the CSV files in the subfolder of `data/transmission_projects/`. New transmission projects, e.g. from TYNDP 2024, can be added in a new subfolder of transmission projects, e.g. `data/transmission_projects/tyndp2024` while extending the list of `transmission_projects` in the `config.yaml` by `tyndp2024`. The CSV files in the project folder should have the same columns as the CSV files in the template folder `data/transmission_projects/template`.

Configuration for `transmission_projects` settings.

{{ schema_table("transmission_projects") }}

**YAML Syntax**

```yaml
{{ yaml_section("transmission_projects") }}
```


## `transformers` {#transformers_cf}

Configuration for `transformers` settings.

{{ schema_table("transformers") }}

**YAML Syntax**

```yaml
{{ yaml_section("transformers") }}
```


## `load` {#load_cf}

Configuration for `load` settings.

{{ schema_table("load") }}

**YAML Syntax**

```yaml
{{ yaml_section("load") }}
```


## `energy` {#energy_cf}

Only used for sector-coupling studies.

Configuration for `energy` settings.

{{ schema_table("energy") }}

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

{{ schema_table("biomass") }}

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

{{ schema_table("solar_thermal") }}

**YAML Syntax**

```yaml
{{ yaml_section("solar_thermal") }}
```

!!! note
    Only used for sector-coupling studies.


## `existing_capacities` {#existing_capacities_cf}

Only used for sector-coupling studies. The value for grouping years are only used in myopic or perfect foresight scenarios.

Configuration for `existing_capacities` settings.

{{ schema_table("existing_capacities") }}

**YAML Syntax**

```yaml
{{ yaml_section("existing_capacities") }}
```

!!! note
    Only used for sector-coupling studies. The value for grouping years are only used in myopic or perfect foresight scenarios.


## `sector` {#sector_cf}

Only used for sector-coupling studies.

Configuration for `sector` settings.

{{ schema_table("sector") }}

**YAML Syntax**

```yaml
{{ yaml_section("sector") }}
```


!!! note
    Only used for sector-coupling studies.


## `industry` {#industry_cf}

Only used for sector-coupling studies.

Configuration for `industry` settings.

{{ schema_table("industry") }}

**YAML Syntax**

```yaml
{{ yaml_section("industry") }}
```

!!! note
    Only used for sector-coupling studies.


## `costs` {#costs_cf}

Configuration for `costs` settings.

{{ schema_table("costs") }}

**YAML Syntax**

```yaml
{{ yaml_section("costs") }}
```


## `clustering` {#clustering_cf}

use `min` in `p_nom_max:` for more conservative assumptions.

Configuration for `clustering` settings.

{{ schema_table("clustering") }}

**YAML Syntax**

```yaml
{{ yaml_section("clustering") }}
```

!!! tip
    use `min` in `p_nom_max:` for more conservative assumptions.


## `adjustments` {#adjustments_cf}

Configuration for top-level adjustments key.

{{ schema_table("adjustments") }}

**YAML Syntax**

```yaml
{{ yaml_section("adjustments") }}
```


## `solving` {#solving_cf}

Configuration for `solving` settings.

{{ schema_table("solving") }}

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

{{ schema_table("overpass_api") }}

**YAML Syntax**

```yaml
{{ yaml_section("overpass_api") }}
```



## `plotting` {#plotting_cf}

```yaml
{{ yaml_section("plotting", source="plotting") }}
```

