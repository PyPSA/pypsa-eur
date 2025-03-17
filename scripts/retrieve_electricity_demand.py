# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve electricity prices from OPSD.
"""

import logging

import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_electricity_demand")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    url = "https://data.open-power-system-data.org/time_series/{version}/time_series_60min_singleindex.csv"

    df1, df2 = [
        pd.read_csv(url.format(version=version), index_col=0)
        for version in snakemake.params.versions
    ]
    combined = pd.concat([df1, df2[df2.index > df1.index[-1]]])

    pattern = "_load_actual_entsoe_transparency"
    transparency = combined.filter(like=pattern).rename(
        columns=lambda x: x.replace(pattern, "")
    )
    pattern = "_load_actual_entsoe_power_statistics"
    powerstatistics = combined.filter(like=pattern).rename(
        columns=lambda x: x.replace(pattern, "")
    )

    res = transparency.fillna(powerstatistics)

    res.to_csv(snakemake.output[0])
