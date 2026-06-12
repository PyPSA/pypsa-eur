<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!---->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Building Electricity Networks

The preparation process of the PyPSA-Eur energy system model consists of a group of `snakemake`
rules which are briefly outlined and explained in detail in the sections below. The pipeline follows
the order `base → simplified → clustered → composed → solved`;
intermediate networks are stored under `resources/{run}/networks`, while solved networks are written to
`results/{run}/networks/solved_{horizon}.nc`.

Not all data dependencies are shipped with the git repository.
Instead we provide separate data bundles which can be obtained
using the `retrieve*` rules ([Retrieving Data](retrieve.md)).
Having downloaded the necessary data, it can build a base PyPSA network with the following rules

- [build_shapes][] generates GeoJSON files with shapes of the countries, exclusive economic zones and [NUTS3](https://en.wikipedia.org/wiki/Nomenclature_of_Territorial_Units_for_Statistics) areas.
- [base_network][] builds and stores the base network with all buses, HVAC lines and HVDC links, and determines [Voronoi cells](https://en.wikipedia.org/wiki/Voronoi_diagram) for all substations.

The network is then simplified by preparing **approximations** of the network model, for which it is computationally viable to co-optimize generation, storage and transmission capacities.


- [simplify_network][] transforms the transmission grid to a 380 kV only equivalent network, while
- [cluster_network][] uses a [k-means](https://en.wikipedia.org/wiki/K-means_clustering) based clustering technique to partition the network into a given number of zones and then reduce the network to a representation with one bus per zone.

The simplification and clustering steps are described in detail in the paper

- Jonas Hörsch and Tom Brown. [The role of spatial scale in joint optimisations of generation and transmission for European highly renewable scenarios](https://arxiv.org/abs/1705.07617)), *14th International Conference on the European Energy Market*, 2017. [arXiv:1705.07617](https://arxiv.org/abs/1705.07617), [doi:10.1109/EEM.2017.7982024](https://doi.org/10.1109/EEM.2017.7982024).

Then, the process continues by calculating conventional power plant capacities, potentials, and per-unit availability time series for variable renewable energy carriers and hydro power plants with the following rules:

- [build_powerplants][] for today's thermal power plant capacities using [powerplantmatching](https://github.com/PyPSA/powerplantmatching) allocating these to the matching clustered region for each powerplant,
- [determine_availability_matrix][] for the land eligibility analysis of each cutout grid cell for PV, onshore and offshore wind,
- [build_renewable_profiles][] for the hourly capacity factors and installation potentials constrained by land-use in each substation's Voronoi cell for PV, onshore and offshore wind, and
- [build_hydro_profile][] for the hourly per-unit hydro power availability time series.

Once `networks/clustered.nc` and the associated bus/line maps exist, the single rule
[compose_network][] stitches electricity, sector, and brownfield inputs together into
`networks/composed_{horizon}.nc`. These per-horizon files are then consumed by
[solve_network][].

## Rule `build_cutout` {#cutout}

::: build_cutout


## Rule `build_osm_boundaries`

::: build_osm_boundaries

## Rule `clean_osm_data`

::: clean_osm_data


## Rule `build_osm_network`

::: build_osm_network

## Rule `build_tyndp_network`

::: build_tyndp_network

## Rule `base_network` {#base}

::: base_network

## Rule `build_natura_raster`

::: build_natura


## Rule `build_transmission_projects`

::: build_transmission_projects

## Rule `build_line_rating`

::: build_line_rating

## Rule `add_transmission_projects_and_dlr`

::: add_transmission_projects_and_dlr

## Rule `build_bidding_zones`

::: build_bidding_zones

## Rule `build_shapes` {#shapes}

::: build_shapes


## Rule `build_electricity_demand_base`

::: build_electricity_demand_base

## Rule `build_electricity_demand` {#electricity_demand}

::: build_electricity_demand

## Rule `build_hac_features`

::: build_hac_features

## Rule `simplify_network` {#simplify}

::: simplify_network

## Rule `cluster_network` {#cluster}

::: cluster_network


## Rule `build_fossil_fuel_prices` {#monthlyprices}

::: build_monthly_prices

## Rule `build_ship_raster` {#ship}

::: build_ship_raster

## Rule `determine_availability_matrix_MD_UA` {#availabilitymatrixmdua}

::: determine_availability_matrix_MD_UA


## Rule `determine_availability_matrix` {#renewableprofiles}

::: determine_availability_matrix


## Rule `build_renewable_profiles`

::: build_renewable_profiles


## Rule `build_hydro_profile` {#hydroprofiles}

::: build_hydro_profile

## Rule `build_powerplants` {#powerplants}

::: build_powerplants

## Rule `compose_network` {#compose}

::: compose_network

## Library modules

The following modules are no longer standalone rules; their functions are
imported and called by [compose_network][] during network assembly.

### `add_electricity`

::: add_electricity

### `prepare_network`

::: prepare_network
