# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
TODO To fill later
"""

# import geojson
import logging
# import numpy as np
# import overpass as op
# import os
# import pandas as pd
# import pypsa
# import requests

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


