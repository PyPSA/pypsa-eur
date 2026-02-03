# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
This rule downloads the load data from `Open Power System Data Time series
<https://data.open-power-system-data.org/time_series/>`_. For all countries in
the network, the per country load timeseries are extracted from the dataset.
After filling small gaps linearly and large gaps by copying time-slice of a
given period, the load data is exported to a ``.csv`` file.
"""

import logging

import numpy as np
import pandas as pd
from pandas import Timedelta as Delta

from scripts._helpers import configure_logging, get_snapshots, set_scenario_config

logger = logging.getLogger(__name__)

OPSD_DATE_FORMAT = "%Y-%m-%dT%H:%M:%SZ"


def load_timeseries(fn, years, countries, date_format=None):
    """
    Read and pre-filter load data.

    Parameters
    ----------
    years : None or slice()
        Years for which to read load data (defaults to
        slice("2018","2019"))
    fn : str
        File name or url location (file format .csv)
    countries : listlike
        Countries for which to read load data.

    Returns
    -------
    load : pd.DataFrame
        Load time-series with UTC timestamps x ISO-2 countries
    """
    return (
        pd.read_csv(fn, index_col=0, parse_dates=[0], date_format=date_format)
        .tz_localize(None)
        .dropna(how="all", axis=0)
        .filter(items=countries)
        .loc[years]
    )


def consecutive_nans(ds):
    return (
        ds.isnull()
        .astype(int)
        .groupby(ds.notnull().astype(int).cumsum()[ds.isnull()])
        .transform("sum")
        .fillna(0)
    )


def fill_large_gaps(ds, shift):
    """
    Fill up large gaps with load data from the previous week.

    This function fills gaps ragning from 3 to 168 hours (one week).
    """
    shift = Delta(shift)
    nhours = shift / np.timedelta64(1, "h")
    if (consecutive_nans(ds) > nhours).any():
        logger.warning(
            "There exist gaps larger then the time shift used for copying time slices."
        )
    time_shift = pd.Series(ds.values, ds.index + shift)
    return ds.where(ds.notnull(), time_shift.reindex_like(ds))


def nan_statistics(df):
    def max_consecutive_nans(ds):
        return (
            ds.isnull()
            .astype(int)
            .groupby(ds.notnull().astype(int).cumsum())
            .sum()
            .max()
        )

    consecutive = df.apply(max_consecutive_nans)
    total = df.isnull().sum()
    max_total_per_month = df.isnull().resample("m").sum().max()
    return pd.concat(
        [total, consecutive, max_total_per_month],
        keys=["total", "consecutive", "max_total_per_month"],
        axis=1,
    )


def copy_timeslice(load, cntry, start, stop, delta, fn_load=None):
    start = pd.Timestamp(start)
    stop = pd.Timestamp(stop)
    if start in load.index and stop in load.index:
        if start - delta in load.index and stop - delta in load.index and cntry in load:
            load.loc[start:stop, cntry] = load.loc[
                start - delta : stop - delta, cntry
            ].values
        elif fn_load is not None and cntry in load:
            duration = pd.date_range(freq="h", start=start - delta, end=stop - delta)
            load_raw = load_timeseries(fn_load, duration, [cntry])
            load.loc[start:stop, cntry] = load_raw.loc[
                start - delta : stop - delta, cntry
            ].values


def manual_adjustment(load, fn_load, countries):
    """
    Adjust gaps manual for load data.

    Parameters
    ----------
     load : pd.DataFrame
         Load time-series with UTC timestamps x ISO-2 countries
    load_fn: str
         File name or url location (file format .csv)

    Returns
    -------
     load : pd.DataFrame
         Manual adjusted and interpolated load time-series with UTC
         timestamps x ISO-2 countries
    """

    copy_timeslice(load, "UA", "2010-01-01 00:00", "2010-01-01 01:00", Delta(days=-1))
    copy_timeslice(load, "XK", "2010-01-01 00:00", "2010-01-01 01:00", Delta(days=-1))
    copy_timeslice(load, "MD", "2010-01-01 00:00", "2010-01-01 01:00", Delta(days=-1))
    copy_timeslice(load, "BA", "2010-01-01 00:00", "2010-01-01 01:00", Delta(days=-1))

    copy_timeslice(load, "CH", "2010-01-19 07:00", "2010-01-19 22:00", Delta(days=1))
    copy_timeslice(load, "CH", "2010-03-28 00:00", "2010-03-28 21:00", Delta(days=1))
    copy_timeslice(load, "CH", "2010-11-04 04:00", "2010-11-04 22:00", Delta(days=1))
    copy_timeslice(load, "CH", "2025-04-30 00:00", "2025-04-30 23:00", Delta(days=1))
    copy_timeslice(load, "CH", "2023-02-16 23:00", "2025-02-17 09:00", Delta(days=1))
    copy_timeslice(load, "NO", "2010-12-09 11:00", "2010-12-09 18:00", Delta(days=1))
    copy_timeslice(load, "PL", "2014-12-08 00:00", "2014-12-08 23:00", Delta(days=1))
    copy_timeslice(load, "PL", "2025-04-28 00:00", "2025-04-28 23:00", Delta(days=1))

    copy_timeslice(load, "AT", "2015-11-29 00:00", "2015-11-29 23:00", Delta(days=-1))
    copy_timeslice(load, "AT", "2014-12-17 23:00", "2014-12-18 07:00", Delta(days=-1))
    copy_timeslice(load, "CH", "2023-12-05 00:00", "2025-12-05 23:00", Delta(days=-1))
    copy_timeslice(load, "SI", "2014-12-30 12:00", "2014-12-31 23:00", Delta(days=-2))

    copy_timeslice(load, "LU", "2018-12-07 00:00", "2018-12-09 23:00", Delta(days=-2))
    copy_timeslice(load, "LU", "2018-12-04 00:00", "2018-12-04 23:00", Delta(days=-1))
    copy_timeslice(load, "LU", "2018-11-15 00:00", "2018-11-15 23:00", Delta(days=-1))

    copy_timeslice(load, "MK", "2019-09-01 00:00", "2019-09-02 23:00", Delta(weeks=-1))
    copy_timeslice(load, "MK", "2024-09-02 00:00", "2024-09-02 23:00", Delta(days=-1))

    copy_timeslice(load, "GB", "2018-12-14 21:00", "2018-12-15 02:00", Delta(days=1))
    copy_timeslice(load, "IE", "2025-11-17 23:00", "2025-11-19 09:00", Delta(weeks=-1))

    # is a WE, so take WE before
    copy_timeslice(load, "CH", "2010-10-08 13:00", "2010-10-10 21:00", Delta(weeks=1))

    # longer stretch
    copy_timeslice(load, "NL", "2014-12-01 00:00", "2014-12-19 23:00", Delta(weeks=-3))

    def _safe_where(df, col, pred):
        if col in df.columns:
            df[col] = df[col].where(pred(df[col]), np.nan)

    def _safe_setna(df, idx, col):
        if col in df.columns:
            df.loc[idx, col] = np.nan

    # remove erroneous peaks
    _safe_where(load, "BA", lambda s: s < 2650)
    _safe_where(load, "DK", lambda s: s < 6400)
    _safe_where(load, "MK", lambda s: s < 2000)
    _safe_where(load, "PT", lambda s: s < 11000)

    # remove erroneous drops
    _safe_where(load, "EE", lambda s: s > 300)
    _safe_where(load, "GB", lambda s: s > 10000)
    _safe_where(load, "GR", lambda s: s > 2300)
    _safe_where(load, "LT", lambda s: s > 600)
    _safe_where(load, "LV", lambda s: s > 300)
    _safe_where(load, "ME", lambda s: s > 90)
    _safe_where(load, "NO", lambda s: s > 6000)
    _safe_where(load, "SE", lambda s: s > 6000)
    _safe_where(load, "SI", lambda s: s > 400)
    _safe_where(load, "XK", lambda s: s > 200)
    _safe_where(load, "PT", lambda s: s > 3000)
    _safe_where(load, "AT", lambda s: s > 3000)
    _safe_where(load, "BG", lambda s: s > 2000)

    _safe_setna(load, "2019-04-03 01:00", "BG")
    _safe_setna(load, "2019-10-29 01:00", "BG")
    _safe_setna(load, "2019-10-29 03:00", "BG")
    _safe_setna(load, "2018-09-12 16:00", "GB")
    _safe_setna(load, "2019-10-28 11:00", "SK")

    copy_timeslice(load, "MK", "2019-04-30 00:00", "2019-05-01 23:00", Delta(weeks=-1))
    for year in [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2024]:
        copy_timeslice(
            load, "BA", f"{year}-07-05 00:00", f"{year}-07-08 00:00", Delta(weeks=-1)
        )
        copy_timeslice(
            load, "BA", f"{year}-11-22 19:00", f"{year}-11-22 21:00", Delta(days=-1)
        )

    return load


def repeat_years(s: pd.Series, years: list) -> pd.Series:
    s = s[~((s.index.month == 2) & (s.index.day == 29))]  # drop leap day
    return pd.concat(
        [s.set_axis(s.index.map(lambda t: t.replace(year=y))) for y in years]
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_electricity_demand")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    # supported year ranges
    years = slice("2010", "2025")

    interpolate_limit = snakemake.params.load["fill_gaps"]["interpolate_limit"]
    countries = snakemake.params.countries

    time_shift = snakemake.params.load["fill_gaps"]["time_shift_for_large_gaps"]

    opsd = load_timeseries(snakemake.input.opsd, years, countries, OPSD_DATE_FORMAT)

    entsoe = load_timeseries(snakemake.input.entsoe, years, countries)

    neso = load_timeseries(snakemake.input.neso, years, countries)

    load = entsoe.combine_first(opsd).combine_first(neso)

    # zero values are considered missing values
    load.replace(0, np.nan, inplace=True)

    if "UA" in countries:
        # use 2021 as template for filling missing years
        s = load.loc["2021", "UA"]
        fill_values = repeat_years(s, range(2010, 2026))
        load["UA"] = load["UA"].combine_first(fill_values)

    if "MD" in countries:
        # use 2022 as template for filling missing years
        s = load.loc["2022", "MD"]
        fill_values = repeat_years(s, range(2010, 2025))
        load["MD"] = load["MD"].combine_first(fill_values)

    if "AL" in countries:
        # use 2024 as template for filling missing and erroneous years
        s = load.loc["2024", "AL"]
        fill_values = repeat_years(s, range(2010, 2026))
        load.loc[fill_values.index, "AL"] = fill_values

    if "BA" in countries:
        # use 2024 as template for filling missing and erroneous years
        s = load.loc["2024", "BA"]
        fill_values = repeat_years(s, range(2010, 2026))
        load["BA"] = load["BA"].combine_first(fill_values)

    if "XK" in countries:
        # use 2024 as template for filling missing and erroneous years
        s = load.loc["2024", "XK"]
        fill_values = repeat_years(s, range(2010, 2024))
        load["XK"] = load["XK"].combine_first(fill_values)

    if "MK" in countries:
        # use 2024 as template for filling missing and erroneous years
        s = load.loc["2024", "MK"]
        fill_values = repeat_years(s, range(2025, 2026))
        load["MK"] = load["MK"].combine_first(fill_values)

    if "CY" in countries:
        # use 2021 as template for filling missing years
        s = load.loc["2021", "CY"]
        fill_values = repeat_years(s, range(2010, 2026))
        load["CY"] = load["CY"].combine_first(fill_values)

    if snakemake.params.load["manual_adjustments"]:
        load = manual_adjustment(load, snakemake.input[0], countries)

    if snakemake.params.load["fill_gaps"]["enable"]:
        logger.info(f"Linearly interpolate gaps of size {interpolate_limit} and less.")
        load = load.interpolate(method="linear", limit=interpolate_limit)

        logger.info(
            f"Filling larger gaps by copying time-slices of period '{time_shift}'."
        )
        load = load.apply(fill_large_gaps, shift=time_shift)

    if snakemake.params.load["supplement_synthetic"]:
        logger.info("Supplement missing data with synthetic data.")
        fn = snakemake.input.synthetic
        synthetic_load = pd.read_csv(fn, index_col=0, parse_dates=True)
        # UA, MD, XK, CY, MT do not appear in synthetic load data
        countries = list(set(countries) - set(["UA", "MD", "XK", "CY", "MT"]))
        synthetic_load = synthetic_load.loc[snapshots, countries]
        load = load.combine_first(synthetic_load)

    assert not load.isna().any().any(), (
        "Load data contains nans. Adjust the parameters "
        "`time_shift_for_large_gaps` or modify the `manual_adjustment` function "
        "for implementing the needed load data modifications."
    )

    fixed_year = snakemake.params["load"].get("fixed_year", False)
    years = (
        slice(str(fixed_year), str(fixed_year))
        if fixed_year
        else slice(snapshots[0], snapshots[-1])
    )

    load = load.loc[years].reindex(index=snapshots)

    # need to reindex load time series to target year
    if fixed_year:
        load.index = load.index.map(lambda t: t.replace(year=snapshots.year[0]))

    load.to_csv(snakemake.output[0])
