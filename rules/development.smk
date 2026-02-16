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
        base_network=resources("osm/comparison/incumbent/networks/base.nc"),
        regions_onshore=resources("osm/comparison/incumbent/regions_onshore.geojson"),
        regions_offshore=resources("osm/comparison/incumbent/regions_offshore.geojson"),
        admin_shapes=resources("osm/comparison/incumbent/admin_shapes.geojson"),
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
        n_incumbent=resources("osm/comparison/incumbent/networks/base.nc"),
        country_shapes=resources("country_shapes.geojson"),
        regions_offshore=resources("regions_offshore.geojson"),
    output:
        lengths=resources("osm/comparison/lengths.pdf"),
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
    message:
        "Preparing OSM network release files and map."
    params:
        line_types=config["lines"]["types"],
        release_version=config_provider("osm_network_release", "release_version"),
        include_polygons=True,
        export=True,
    input:
        base_network=resources("networks/base.nc"),
        stations_polygon=resources("osm/build/geojson/stations_polygon.geojson"),
        buses_polygon=resources("osm/build/geojson/buses_polygon.geojson"),
    output:
        buses=resources("osm/release/buses.csv"),
        converters=resources("osm/release/converters.csv"),
        lines=resources("osm/release/lines.csv"),
        links=resources("osm/release/links.csv"),
        transformers=resources("osm/release/transformers.csv"),
        map=resources("osm/release/map.html"),
    log:
        logs("prepare_osm_network_release.log"),
    benchmark:
        benchmarks("prepare_osm_network_release")
    threads: 1
    resources:
        mem_mb=1000,
    script:
        "../scripts/prepare_osm_network_release.py"


rule map_incumbent:
    message:
        "Preparing map of incumbent network for comparison with OSM release."
    params:
        line_types=config["lines"]["types"],
        release_version="Incumbent",
        include_polygons=False,
        export=False,
    input:
        base_network=resources("osm/comparison/incumbent/networks/base.nc"),
    output:
        map=resources("osm/comparison/map_incumbent.html"),
    log:
        logs("prepare_osm_network_release.log"),
    benchmark:
        benchmarks("prepare_osm_network_release")
    threads: 1
    resources:
        mem_mb=1000,
    script:
        "../scripts/prepare_osm_network_release.py"


rule osm_release:
    message:
        "Creating OSM network release files, map and comparison with incumbent network."
    input:
        resources("osm/release/buses.csv"),
        resources("osm/release/converters.csv"),
        resources("osm/release/lines.csv"),
        resources("osm/release/links.csv"),
        resources("osm/release/transformers.csv"),
        resources("osm/release/map.html"),
        resources("osm/comparison/map_incumbent.html"),
        resources("osm/comparison/lengths.pdf"),
