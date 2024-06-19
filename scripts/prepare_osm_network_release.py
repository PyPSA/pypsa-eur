# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import logging

import pypsa
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def prepare_osm_network_release(network):
    return None


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_osm_network_release")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    network = pypsa.Network(snakemake.input.base_network)

    buses_columns = [
        "bus_id",
        "voltage",
        "dc",
        "symbol",
        "under_construction",
        "x",
        "y",
        "country",
        "geometry",
    ]

    lines_columns = [
        "line_id",
        "bus0",
        "bus1",
        "voltage",
        "circuits",
        "length",
        "underground",
        "under_construction",
        "geometry",
    ]

    links_columns = [
        "link_id",
        "bus0",
        "bus1",
        "voltage",
        "p_nom",
        "length",
        "underground",
        "under_construction",
        "geometry",
    ]

    transformers_columns = [
        "transformer_id",
        "bus0",
        "bus1",
        "voltage_bus0",
        "voltage_bus1",
        "geometry",
    ]

    converters_columns = []
