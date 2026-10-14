<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Transmission Grid

The rules on this page build the transmission network from OpenStreetMap or the
ENTSO-E reference grid, add planned projects and dynamic line rating, and
reduce the network to the model regions. They follow the order `base →
simplified → clustered`; the clustered network is the starting point of every
model run. The simplification and clustering steps are described in
[@horschRoleSpatial2017] and [@frysztackiComparisonClustering2022].

## Base network

{{ rules("clean_osm_data", "build_osm_network", "build_tyndp_network", "base_network") }}

## Transmission projects and line rating

{{ rules("build_transmission_projects", "build_line_rating", "add_transmission_projects_and_dlr") }}

## Simplification and clustering

{{ rules("build_hac_features", "simplify_network", "cluster_network", "chain_busmaps") }}
