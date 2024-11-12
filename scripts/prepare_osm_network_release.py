# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import logging
import os

import pandas as pd
import pypsa
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


BUSES_COLUMNS = [
    "bus_id",
    "voltage",
    "dc",
    "symbol",
    "under_construction",
    "tags",
    "x",
    "y",
    "country",
    "geometry",
]
LINES_COLUMNS = [
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
LINKS_COLUMNS = [
    "link_id",
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "length",
    "underground",
    "under_construction",
    "tags",
    "geometry",
]
TRANSFORMERS_COLUMNS = [
    "transformer_id",
    "bus0",
    "bus1",
    "voltage_bus0",
    "voltage_bus1",
    "s_nom",
    "geometry",
]
CONVERTERS_COLUMNS = [
    "converter_id",
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "geometry",
]


def export_clean_csv(df, columns, output_file):
    """
    Export a cleaned DataFrame to a CSV file.

    Args:
        df (pandas.DataFrame): The DataFrame to be exported.
        columns (list): A list of column names to include in the exported CSV file.
        output_file (str): The path to the output CSV file.

    Returns:
        None
    """
    rename_dict = {
        "Bus": "bus_id",
        "Line": "line_id",
        "Link": "link_id",
        "Transformer": "transformer_id",
        "v_nom": "voltage",
        "num_parallel": "circuits",
    }

    if "converter_id" in columns:
        rename_dict["Link"] = "converter_id"

    df.reset_index().rename(columns=rename_dict).loc[:, columns].replace(
        {True: "t", False: "f"}
    ).to_csv(output_file, index=False, quotechar="'")

    return None


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_osm_network_release")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    network = pypsa.Network(snakemake.input.base_network)

    network.buses["dc"] = network.buses.pop("carrier").map({"DC": "t", "AC": "f"})
    network.lines.length = network.lines.length * 1e3
    network.links.length = network.links.length * 1e3

    # Sort alphabetically
    network.buses.sort_index(inplace=True)
    network.transformers.sort_index(inplace=True)
    network.lines.sort_index(inplace=True)
    network.links.sort_index(inplace=True)

    # Export to clean csv for release
    logger.info(f"Exporting {len(network.buses)} buses to %s", snakemake.output.buses)
    export_clean_csv(network.buses, BUSES_COLUMNS, snakemake.output.buses)

    logger.info(
        f"Exporting {len(network.transformers)} transformers to %s",
        snakemake.output.transformers,
    )
    export_clean_csv(
        network.transformers, TRANSFORMERS_COLUMNS, snakemake.output.transformers
    )

    logger.info(f"Exporting {len(network.lines)} lines to %s", snakemake.output.lines)
    export_clean_csv(network.lines, LINES_COLUMNS, snakemake.output.lines)

    # Boolean that specifies if link element is a converter
    is_converter = network.links.index.str.startswith("conv") == True

    logger.info(
        f"Exporting {len(network.links[~is_converter])} links to %s",
        snakemake.output.links,
    )
    export_clean_csv(
        network.links[~is_converter], LINKS_COLUMNS, snakemake.output.links
    )

    logger.info(
        f"Exporting {len(network.links[is_converter])} converters to %s",
        snakemake.output.converters,
    )
    export_clean_csv(
        network.links[is_converter], CONVERTERS_COLUMNS, snakemake.output.converters
    )

    logger.info("Export of OSM network for release complete.")
