<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Regions and Weather

The rules on this page define the geography of the model: the country and
administrative shapes, the offshore zones and the weather cutouts from which
all time series are derived. Everything downstream, from the base network to
the population layouts, is built on these shapes.

## Weather data

{{ rules("build_cutout") }}

## Shapes

{{ rules("build_shapes", "build_nuts3_shapes", "build_offshore_shapes", "build_osm_boundaries", "build_bidding_zones") }}
