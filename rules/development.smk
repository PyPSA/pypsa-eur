# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

if config["electricity"]["base_network"] == "osm-raw":

    rule prepare_osm_network_release:
        input:
            base_network=resources("networks/base.nc"),
        output:
            buses=resources("osm-raw/release/buses.csv"),
            converters=resources("osm-raw/release/converters.csv"),
            lines=resources("osm-raw/release/lines.csv"),
            links=resources("osm-raw/release/links.csv"),
            transformers=resources("osm-raw/release/transformers.csv"),
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
