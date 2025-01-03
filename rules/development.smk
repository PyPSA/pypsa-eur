# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

if config["electricity"]["base_network"] == "osm-raw":

    rule prepare_osm_network_release:
        params:
            line_types=config["lines"]["types"],
        input:
            base_network=resources("networks/base.nc"),
            stations_polygon=resources("osm-raw/build/geojson/stations_polygon.geojson"),
            buses_polygon=resources("osm-raw/build/geojson/buses_polygon.geojson"),
        output:
            buses=resources("osm-raw/release/buses.csv"),
            converters=resources("osm-raw/release/converters.csv"),
            lines=resources("osm-raw/release/lines.csv"),
            links=resources("osm-raw/release/links.csv"),
            transformers=resources("osm-raw/release/transformers.csv"),
            map=resources("osm-raw/release/map.html"),
        log:
            logs("prepare_osm_network_release.log"),
        benchmark:
            benchmarks("prepare_osm_network_release")
        threads: 1
        resources:
            mem_mb=1000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/prepare_osm_network_release.py"
