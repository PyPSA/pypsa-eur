# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
TODO To fill later
"""

# import geojson
import geopandas as gpd
import logging
# import numpy as np
# import os
import pandas as pd
# import pypsa
# import requests
import tqdm.auto as tqdm

from _helpers import configure_logging
logger = logging.getLogger(__name__)

def clean_osm_data(output):
    with open(output, "w") as file:
        file.write("Hello, world!\n")


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("clean_osm_data")
    
    configure_logging(snakemake)
    logger.info("Dummy log: clean_osm_data()")

    output = str(snakemake.output)
    clean_osm_data(output)

    # Create df by iterating over lines_way and append them to df_lines_way
    gdf1 = gpd.read_file(snakemake.input["lines_way"])
    
    
    snakemake.wildcards
    snakemake.input["lines_way"].keys()

