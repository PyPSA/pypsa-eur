# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build CO2 price time series.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_co2_prices")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    df = pd.read_csv(snakemake.input["csv"], parse_dates=True, index_col=0)["price"]

    rolling_window = snakemake.params["rolling_window"]

    df_smoothed = (
        df.rolling(window=rolling_window, center=True, min_periods=1)
        .mean()
        .resample("1h")
        .mean()
        .ffill()
    )

    df_smoothed.index = df_smoothed.index.tz_localize(None)

    df_smoothed.to_csv(snakemake.output["csv"])
