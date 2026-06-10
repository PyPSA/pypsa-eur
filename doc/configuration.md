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

