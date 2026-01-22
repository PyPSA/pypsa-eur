# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

if (
    config["electricity"]["base_network"] == "osm"
    and config["data"]["osm"]["source"] == "build"
):

    rule prepare_osm_network_release:
        params:
            line_types=config["lines"]["types"],
        input:
            base_network=resources("networks/base.nc"),
            stations_polygon=resources("osm/geojson/stations_polygon.geojson"),
            buses_polygon=resources("osm/geojson/buses_polygon.geojson"),
        output:
            buses=resources("osm/upstream/release/buses.csv"),
            converters=resources("osm/upstream/release/converters.csv"),
            lines=resources("osm/upstream/release/lines.csv"),
            links=resources("osm/upstream/release/links.csv"),
            transformers=resources("osm/upstream/release/transformers.csv"),
            map=resources("osm/upstream/release/map.html"),
        log:
            logs("prepare_osm_network_release.log"),
        benchmark:
            benchmarks("prepare_osm_network_release")
        threads: 1
        resources:
            mem_mb=1000,
        script:
            scripts("/prepare_osm_network_release.py")
