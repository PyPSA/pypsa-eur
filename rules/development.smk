# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule prepare_osm_network_release:
    input:
        base_network=resources("networks/base.nc"),
    output:
        buses=resources("osm/release/buses.csv"),
        converters=resources("osm/release/converters.csv"),
        lines=resources("osm/release/lines.csv"),
        links=resources("osm/release/links.csv"),
        transformers=resources("osm/release/transformers.csv"),
    log:
        logs("prepare_osm_network_release.log"),
    benchmark:
        benchmarks("prepare_osm_network_release")
    script:
        "../scripts/prepare_osm_network_release.py"
