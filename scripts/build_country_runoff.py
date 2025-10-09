# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build daily hydro runoff for each country to fill missing EIA statistics for hydro generation.

Outputs
-------

- ``data/country_runoff/build/unknown/era5-runoff-per-country.csv``:

    ===================  ==========  ===========  =========================================================
    Field                Dimensions  Unit         Description
    ===================  ==========  ===========  =========================================================
    index/time           time        day          Datestamp, YYYY-MM-DD
    -------------------  ----------  -----------  ---------------------------------------------------------
    <columns>            country     ISO-3166 A2  Daily total runoff (volume per area) per country
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

        snakemake = mock_snakemake("build_country_runoff")
    configure_logging(snakemake)

    cutout: atlite.Cutout = load_cutout(snakemake.input.cutouts)
    country_shapes = gpd.read_file(snakemake.input.country_shapes).set_index("name")[
        "geometry"
    ]

    ds = cutout.runoff(shapes=country_shapes)

    df: pd.DataFrame = ds.to_pandas()
    df = df.resample("D").sum().astype(int)

    df.to_csv(snakemake.output.era5_runoff)
