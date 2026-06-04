<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!---->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Tutorial: Sector-Coupled {#tutorial_sector}

!!! note
    If you have not done it yet, follow the [installation](installation.md) steps first.

    Also, checkout the tutorial for electricity-only systems first at [Tutorial: Electricity-Only](tutorial.md#tutorial).

In this tutorial, we will add further sectors to the electricity-only model from
[Tutorial: Electricity-Only](tutorial.md#tutorial), namely industry, transport, and buildings. This
requires processing of a few more raw data sources.

The sector-coupling code can be run as an overnight / greenfield scenario or
with multi-horizon investment with myopic foresight. Pathway analysis with
perfect foresight is under development. See also the documentation on
[Foresight](configuration.md#foresight_cf).

## Overnight Scenarios

### Configuration

The default configuration file (``config/config.default.yaml``) is set up for running
overnight scenarios. Running a sector-coupled model unlocks many further
configuration options. In the example below, we say that the gas network should
be added and spatially resolved. We also say that the existing gas network may
be retrofitted to transport hydrogen instead.

```yaml
{{ yaml_section("sector", source="test/config.overnight.yaml") }}
```

Documentation for all options will be added successively to [Configuration](configuration.md).

Scenarios can be defined like for electricity-only studies, but with additional
wildcard options.

```yaml
{{ yaml_section("scenario", source="test/config.overnight.yaml") }}
```

For allowed wildcard values, refer to [Wildcards](wildcards.md).

### Execution

To run an overnight / greenfiled scenario with the specifications above, run

```console
$ snakemake -call all --configfile config/test/config.overnight.yaml
```

which will result in the following jobs ``snakemake`` wants to run, some of
which were already included in the electricity-only tutorial:

```console
job                                                 count
------------------------------------------------  -------
add_electricity                                         1
add_transmission_projects_and_dlr                       1
all                                                     1
base_network                                            1
build_ammonia_production                                1
build_biomass_potentials                                1
build_central_heating_temperature_profiles              1
build_clustered_population_layouts                      1
build_cop_profiles                                      1
build_daily_heat_demand                                 1
build_direct_heat_source_utilisation_profiles           1
build_district_heat_share                               1
build_electricity_demand                                1
build_electricity_demand_base                           1
build_energy_totals                                     1
build_gas_input_locations                               1
build_gas_network                                       1
build_heat_totals                                       1
build_hourly_heat_demand                                1
build_industrial_distribution_key                       1
build_industrial_energy_demand_per_country_today        1
build_industrial_energy_demand_per_node                 1
build_industrial_energy_demand_per_node_today           1
build_industrial_production_per_country                 1
build_industrial_production_per_country_tomorrow        1
build_industrial_production_per_node                    1
build_industry_sector_ratios                            1
build_industry_sector_ratios_intermediate               1
build_osm_boundaries                                    4
build_population_layouts                                1
build_population_weighted_energy_totals                 2
build_powerplants                                       1
build_renewable_profiles                                6
build_salt_cavern_potentials                            1
build_shapes                                            1
build_ship_raster                                       1
build_shipping_demand                                   1
build_temperature_profiles                              1
build_transmission_projects                             1
build_transport_demand                                  1
cluster_gas_network                                     1
cluster_network                                         1
determine_availability_matrix                           6
make_summary                                            1
plot_gas_network                                        1
plot_hydrogen_network                                   1
plot_power_network                                      1
plot_power_network_clustered                            1
plot_summary                                            1
prepare_network                                         1
prepare_sector_network                                  1
retrieve_cost_data                                      1
retrieve_databundle                                     1
retrieve_eez                                            1
retrieve_electricity_demand                             1
retrieve_eurostat_data                                  1
retrieve_eurostat_household_data                        1
retrieve_gas_infrastructure_data                        1
retrieve_gem_europe_gas_tracker                         1
retrieve_gem_steel_plant_tracker                        1
retrieve_hotmaps_industrial_sites                       1
retrieve_jrc_ardeco                                     1
retrieve_jrc_enspreso_biomass                           1
retrieve_jrc_idees                                      1
retrieve_nuts_2013_shapes                               1
retrieve_nuts_2021_shapes                               1
retrieve_osm_boundaries                                 4
retrieve_osm_prebuilt                                   1
retrieve_ship_raster                                    1
retrieve_synthetic_electricity_demand                   1
retrieve_usgs_ammonia_production                        1
retrieve_worldbank_urban_population                     1
simplify_network                                        1
solve_sector_network                                    1
time_aggregation                                        1
total                                                  92
```

This covers the retrieval of additional raw data from online resources and
preprocessing data about the transport, industry, and heating sectors as well as
additional rules about geological storage and sequestration potentials, gas
infrastructure, and biomass potentials. The collection rule ``all`` will also
generate summary CSV files and plots after the network has been solved
successfully.

[![Sector-coupled overnight DAG](img/dag_overnight.svg)](img/dag_overnight.svg)

## Myopic Foresight Scenarios

### Configuration

To activate the myopic foresight mode, set

```yaml
{{ yaml_section("foresight", source="test/config.myopic.yaml") }}
```

Scenarios can be defined like for electricity-only studies, but with additional
wildcard options. For the myopic foresight mode, the ``{planning_horizons}`` wildcard
defines the sequence of investment horizons.

```yaml
{{ yaml_section("scenario", source="test/config.myopic.yaml") }}
```

For allowed wildcard values, refer to [Wildcards](wildcards.md).

In the myopic foresight mode, you can tweak for instance exogenously given transition paths, like the one for
the share of primary steel production we change below:

```yaml
{{ yaml_section("industry", source="test/config.myopic.yaml") }}
```

Documentation for all options will be added successively to [Configuration](configuration.md).

### Execution

To run a myopic foresight scenario with the specifications above, run

```console
$ snakemake -call all --configfile config/test/config.myopic.yaml
```

which will result in additional jobs ``snakemake`` wants to run, which
translates to the following workflow diagram which nicely outlines how the
sequential pathway optimisation with myopic foresight is implemented in the
workflow:

[![Sector-coupled myopic DAG](img/dag_myopic.svg)](img/dag_myopic.svg)

!!! note
    The above DAG is abbreviated. The full DAG for the myopic foresight scenario
    contains all intermediate rules for each planning horizon (2030, 2040, 2050).

## Scaling-Up

If you now feel confident and want to tackle runs with larger temporal, technological and
spatial scope, clean-up the repository and after modifying the ``config/config.yaml`` file
target the collection rule ``all`` again without providing the test
configuration file.

```console
$ snakemake -call purge
$ snakemake -call all
```

!!! note
    It is good practice to perform a dry-run using the option `-n`, before you
    commit to a run:

    ```console
    $ snakemake -call all -n
    ```
