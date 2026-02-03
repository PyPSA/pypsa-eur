# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


def get_osm_network_incumbent(
    version: str,
    float=0.6,
    source: str = "archive",
) -> pd.Series:

    fp = workflow.source_path("../data/versions.csv")
    data_versions = load_data_versions(fp)
    name = "osm"

    dataset = data_versions.loc[
        (data_versions["dataset"] == name)
        & (data_versions["source"] == source)
        & (data_versions["supported"])  # Limit to supported versions only
        & (data_versions["version"] == version if "latest" != version else True)
        & (data_versions["latest"] if "latest" == version else True)
    ]

    if dataset.empty:
        raise ValueError(
            f"OSM network for version '{version}' not found in data/versions.csv."
        )

    # Return single-row DataFrame as a Series
    dataset = dataset.squeeze()

    # Generate output folder path in the `data` directory
    dataset["folder"] = Path(
        "data", name, dataset["source"], dataset["version"]
    ).as_posix()

    return dataset


def input_base_network_incumbent(w):
    version = config_provider("osm_network_release", "compare_to", "version")(w)
    source = config_provider("osm_network_release", "compare_to", "source")(w)
    osm_dataset = get_osm_network_incumbent(version, source)
    osm_path = osm_dataset["folder"]
    components = {"buses", "lines", "links", "converters", "transformers"}

    inputs = {c: f"{osm_path}/{c}.csv" for c in components}

    return inputs


OSM_DATASET_INCUMBENT = get_osm_network_incumbent(
    version=config.get("osm_network_release", {})
    .get("compare_to", {})
    .get("version", "0.6"),
    source=config.get("osm_network_release", {})
    .get("compare_to", {})
    .get("source", "archive"),
)

OSM_ARCHIVE_FILES_INCUMBENT = [
    "buses.csv",
    "converters.csv",
    "lines.csv",
    "links.csv",
    "transformers.csv",
    # Newer versions include the additional map.html file for visualisation
    *(["map.html"] if float(OSM_DATASET_INCUMBENT["version"]) >= 0.6 else []),
]


rule retrieve_osm_archive_incumbent:
    message:
        "Retrieving OSM archive data"
    input:
        **{
            file: storage(f"{OSM_DATASET_INCUMBENT['url']}/{file}")
            for file in OSM_ARCHIVE_FILES_INCUMBENT
        },
    output:
        **{
            file: f"{OSM_DATASET_INCUMBENT['folder']}/{file}"
            for file in OSM_ARCHIVE_FILES_INCUMBENT
        },
    log:
        "logs/retrieve_osm_archive.log",
    threads: 1
    resources:
        mem_mb=500,
    run:
        for key in input.keys():
            copy2(input[key], output[key])


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
        regions_onshore=resources("osm-network/comparison/incumbent/regions_onshore.geojson"),
        regions_offshore=resources("osm-network/comparison/incumbent/regions_offshore.geojson"),
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
