# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build country-level heating degree days in Europe for each country. Used for rescaling heat demand in weather years not covered by energy balance statistics.

Outputs
-------

- ``data/country_runoff/build/unknown/era5-hdd-per-country.csv``:

    ===================  ==========  ===========  =========================================================
    Field                Dimensions  Unit         Description
    ===================  ==========  ===========  =========================================================
    index/time           time        day          Datestamp, YYYY-MM-DD
    -------------------  ----------  -----------  ---------------------------------------------------------
    <columns>            country     ISO-3166 A2  Aggregated HDDs per country
    ===================  ==========  ===========  =========================================================

"""

import logging

import atlite
import geopandas as gpd
import pandas as pd

from scripts._helpers import (
    configure_logging,
    load_cutout,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_country_hdd")
    configure_logging(snakemake)

    cutout: atlite.Cutout = load_cutout(snakemake.input.cutouts)
    country_shapes = gpd.read_file(snakemake.input.country_shapes).set_index("name")[
        "geometry"
    ]

    ds = cutout.heat_demand(shapes=country_shapes)

    df: pd.DataFrame = ds.to_pandas().astype(int)

    df.to_csv(snakemake.output.era5_hdd)
