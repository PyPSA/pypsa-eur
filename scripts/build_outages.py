# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve power plant outage data from ENTSOE.
"""

import ast
import logging

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def outage_as_dense(df: pd.DataFrame) -> pd.DataFrame:
    """
    Expand outage data to a dense DataFrame with hourly time index and
    production_resource_id columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing ENTSOE outage data

    Returns
    -------
    pd.DataFrame
    """

    if df.empty:
        return pd.DataFrame()

    df["start_hour"] = df.start.dt.tz_convert("UTC").dt.floor("h")
    df["end_hour"] = df.end.dt.tz_convert("UTC").dt.floor("h")
    df["avail_qty"] = df["avail_qty"].astype(float)
    df["unavail_mw"] = df["nominal_power"] - df["avail_qty"]

    nominal = (
        df.drop_duplicates("production_resource_id")
        .set_index("production_resource_id")
        .nominal_power
    )

    def _expand_outage(r: pd.Series) -> pd.Series:
        return pd.Series(
            r.unavail_mw,
            index=pd.date_range(r.start_hour, r.end_hour, freq="h", inclusive="both"),
        )

    outages_t = df.apply(_expand_outage, axis=1)
    outages_t.index = df.production_resource_id

    outages_t = (
        outages_t.groupby(level=0).sum(min_count=1).T.clip(upper=nominal, axis=1)
    )

    return outages_t


def entsoe_to_powerplantmatching_id(c: str, ppl: pd.DataFrame) -> float | int:
    """
    Map ENTSOE IDs (EIC) to Powerplantmatching IDs.

    Parameters
    ----------
    c : str
        ENTSOE ID (EIC)
    ppl : pd.DataFrame
        Powerplantmatching DataFrame

    Returns
    -------
    float | int
        Powerplantmatching ID or NaN if not found
    """
    idx = ppl.loc[ppl.ENTSOE.str.contains(c, na=False)].index
    return idx[0] if len(idx) == 1 else np.nan


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_outages")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    units = pd.read_csv(snakemake.input.outages, index_col=0)

    ppl = pd.read_csv(snakemake.input.powerplants, index_col=0)
    ppl["ENTSOE"] = (
        ppl["projectID"]
        .map(lambda x: ast.literal_eval(x).get("ENTSOE", np.nan))
        .combine_first(ppl["EIC"])
        .astype(str)
    )

    outages = []

    planned = units.query(
        "businesstype == 'Planned maintenance' and docstatus != 'Cancelled'"
    ).copy()
    outages.append(outage_as_dense(planned))

    unplanned = units.query("businesstype == 'Unplanned outage'").copy()
    outages.append(outage_as_dense(unplanned))

    outages_t = (
        pd.concat(outages, axis=1)
        .T.groupby(level=0)
        .sum()
        .T.replace(0, np.nan)
        .dropna(how="all")
        .dropna(how="all", axis=1)
    )
    mapped_columns = outages_t.columns.map(
        lambda c: entsoe_to_powerplantmatching_id(c, ppl)
    )
    outages_t_mapped = outages_t.loc[:, mapped_columns.notna()]
    outages_t_unmapped = outages_t.loc[:, mapped_columns.isna()]
    outages_t_mapped.columns = mapped_columns.dropna()
    outages_t_mapped.columns.name = "powerplantmatching_id"

    outages_t_mapped.to_csv(snakemake.output.mapped)
    outages_t_unmapped.to_csv(snakemake.output.unmapped)
