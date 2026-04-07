# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Preprocess energy balances from Eurostat.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def build_eurostat(
    input_eurostat: str,
) -> pd.DataFrame:
    """
    Return long-format energy balances for all countries in TWh/a.

    Parameters
    ----------
    input_eurostat : str
        Path to the Eurostat database.

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame containing energy data for all countries in TWh/a.
    """

    # Python engine to allow multiple delimiters, but slower
    raw = pd.read_csv(
        input_eurostat,
        sep=",|\t| [^ ]?\t",
        engine="python",
        na_values=[":", ": m", 0.0],
    ).rename(columns={"geo\\TIME_PERIOD": "country"})

    df = raw.melt(
        id_vars=["freq", "nrg_bal", "siec", "unit", "country"],
        var_name="year",
        value_name="value",
    ).query("unit == 'GWH'")

    # Clean up formatting
    df["year"] = df["year"].astype(int)
    df.drop(["freq", "unit"], axis=1, inplace=True)

    # For total and fossil energy, fill in missing values with
    # closest non-missing value in year index
    for country in df.country.unique():
        mask = (
            (df["country"] == country)
            & (df["nrg_bal"] == "FC_TRA_DAVI_E")
            & (df["siec"].isin(["TOTAL", "FE"]))
        )

        df.loc[mask, "value"] = (
            df[mask].groupby("siec")["value"].transform(lambda x: x.ffill().bfill())
        )

    # Follow JRC-IDEES country code convention
    df["country"] = df["country"].replace({"UK": "GB", "EL": "GR"})

    # Convert from GWh to TWh
    df["value"] = df["value"] / 1000

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_eurostat_balances")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    eurostat = build_eurostat(snakemake.input.tsv_gz)

    eurostat.to_csv(snakemake.output.csv, index=False)
