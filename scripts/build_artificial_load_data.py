# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2022 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This rule downloads the load data.
"""

import logging

logger = logging.getLogger(__name__)
import pandas as pd
from _helpers import configure_logging

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_artificial_load_data", weather_year="")

    configure_logging(snakemake)

    weather_year = snakemake.wildcards.weather_year
    if weather_year:
        snapshots = dict(
            start=weather_year, end=str(int(weather_year) + 1), inclusive="left"
        )
    else:
        snapshots = snakemake.config["snapshots"]
    snapshots = pd.date_range(freq="h", **snapshots)

    fixed_year = snakemake.config["load"].get("fixed_year", False)
    years = (
        slice(str(fixed_year), str(fixed_year))
        if fixed_year
        else slice(snapshots[0], snapshots[-1])
    )
    countries = snakemake.config["countries"]

    load = pd.read_csv(snakemake.input[0], index_col=0, parse_dates=True).loc[
        snapshots, countries
    ]

    assert not load.isna().any().any(), "Load data contains nans."

    if fixed_year:
        load.index = load.index.map(lambda t: t.replace(year=snapshots.year[0]))

    load.to_csv(snakemake.output[0])
