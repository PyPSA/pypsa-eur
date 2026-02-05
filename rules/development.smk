# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule base_network_incumbent:
    message:
        "Building base network to which to compare against."
    params:
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        lines=config_provider("lines"),
        links=config_provider("links"),
        transformers=config_provider("transformers"),
        clustering=config_provider("clustering", "mode"),
        admin_levels=config_provider("clustering", "administrative"),
    input:
        unpack(input_base_network_incumbent),
        nuts3_shapes=resources("nuts3_shapes.geojson"),
        country_shapes=resources("country_shapes.geojson"),
        offshore_shapes=resources("offshore_shapes.geojson"),
        europe_shape=resources("europe_shape.geojson"),
    output:
        base_network=resources("osm-network/comparison/incumbent/networks/base.nc"),
        regions_onshore=resources(
            "osm-network/comparison/incumbent/regions_onshore.geojson"
        ),
        regions_offshore=resources(
            "osm-network/comparison/incumbent/regions_offshore.geojson"
        ),
        admin_shapes=resources("osm-network/comparison/incumbent/admin_shapes.geojson"),
    log:
        logs("base_network_incumbent.log"),
    benchmark:
        benchmarks("base_network_incumbent")
    threads: 4
    resources:
        mem_mb=2000,
    script:
        "../scripts/base_network.py"


rule make_network_comparison:
    message:
        "Create network comparison between two PyPSA networks."
    params:
        countries=config_provider("countries"),
        base_network=config_provider("electricity", "base_network"),
        compare_to_version=config_provider(
            "osm_network_release", "compare_to", "version"
        ),
        voltages=config_provider("electricity", "voltages"),
    input:
        n_release=resources("networks/base.nc"),
        n_incumbent=resources("osm-network/comparison/incumbent/networks/base.nc"),
    output:
        lengths=resources("osm-network/comparison/lengths.pdf"),
    log:
        logs("make_network_comparison.log"),
    benchmark:
        benchmarks("make_network_comparison")
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "../scripts/make_network_comparison.py"


rule prepare_osm_network_release:
    params:
        line_types=config["lines"]["types"],
        release_version=config_provider("osm_network_release", "release_version"),
    input:
        base_network=resources("networks/base.nc"),
        stations_polygon=resources("osm-network/build/geojson/stations_polygon.geojson"),
        buses_polygon=resources("osm-network/build/geojson/buses_polygon.geojson"),
    output:
        buses=resources("osm-network/release/buses.csv"),
        converters=resources("osm-network/release/converters.csv"),
        lines=resources("osm-network/release/lines.csv"),
        links=resources("osm-network/release/links.csv"),
        transformers=resources("osm-network/release/transformers.csv"),
        map=resources("osm-network/release/map.html"),
    log:
        logs("prepare_osm_network_release.log"),
    benchmark:
        benchmarks("prepare_osm_network_release")
    threads: 1
    resources:
        mem_mb=1000,
    script:
        "../scripts/prepare_osm_network_release.py"
