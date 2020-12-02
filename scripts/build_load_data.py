# coding: utf-8
"""

This rule downloads the load data from `Open Power System Data Time series <https://data.open-power-system-data.org/time_series/>`_. For all countries in the network, the per country load timeseries with suffix ``_load_actual_entsoe_transparency`` are extracted from the dataset. After filling small gaps linearly and large gaps by copying time-slice of a given period, the load data is exported to a ``.csv`` file.

Relevant Settings
-----------------

.. code:: yaml

    snapshots:

    load:
        url:
        interpolate_limit:
        time_shift_for_large_gaps:
        manual_adjustments: true


.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`load_cf`

Inputs
------


Outputs
-------

- ``resource/load.csv``:


"""

import logging
logger = logging.getLogger(__name__)
from _helpers import configure_logging

import pandas as pd
import numpy as np


def load_timeseries(fn, years, countries):
    """
    Read load data from OPSD time-series package version 2020-10-06.

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
    logger.info(f"Retrieving load data from '{fn}'.")

    rename = lambda s: s[:-len('_load_actual_entsoe_transparency')]
    return (pd.read_csv(fn, index_col=0, parse_dates=True)
            .filter(like='_load_actual_entsoe_transparency')
            .rename(columns=rename)
            .dropna(how="all", axis=0)
            .rename(columns={'GB_UKM' : 'GB'})
            .filter(items=countries)
            .loc[years])


def consecutive_nans(ds):
    return (ds.isnull().astype(int)
            .groupby(ds.notnull().astype(int).cumsum()[ds.isnull()])
            .transform('sum').fillna(0))


def fill_large_gaps(ds, shift):
    """
    Fill up large gaps with load data from the previous week.

    This function fills gaps ragning from 3 to 168 hours (one week).
    """
    shift = pd.Timedelta(shift)
    nhours = shift / np.timedelta64(1, 'h')
    if (consecutive_nans(ds) > nhours).any():
        logger.warning('There exist gaps larger then the time shift used for '
                       'copying time slices.')
    time_shift = pd.Series(ds.values, ds.index + shift)
    return ds.where(ds.notnull(), time_shift.reindex_like(ds))


def nan_statistics(df):
    def max_consecutive_nans(ds):
        return (ds.isnull().astype(int)
                  .groupby(ds.notnull().astype(int).cumsum())
                  .sum().max())
    consecutive = df.apply(max_consecutive_nans)
    total = df.isnull().sum()
    max_total_per_month = df.isnull().resample('m').sum().max()
    return pd.concat([total, consecutive, max_total_per_month],
                 keys=['total', 'consecutive', 'max_total_per_month'], axis=1)


def copy_timeslice(load, cntry, start, stop, delta):
    start = pd.Timestamp(start)
    stop = pd.Timestamp(stop)
    if start in load.index and stop in load.index:
        load.loc[start:stop, cntry] = load.loc[start-delta:stop-delta, cntry].values
    return load


def manual_adjustment(load, source="ENTSOE_power_statistics"):
    """
    Adjust gaps manual for load data from OPSD time-series package.


    Albania and Macedonia do not exist in the data set. Both get the same load
    curve as Montenegro,  scaled by correspoding ratio of total energy
    consumptions reported by  <IEA Data browser `https://www.iea.org/data-and-statistics?country=WORLD&fuel=Electricity%20and%20heat&indicator=TotElecCons>`_ for the
    year 2016.


    Parameters
    ----------
    load : pd.DataFrame
        Load time-series with UTC timestamps x ISO-2 countries

    Returns
    -------
    load : pd.DataFrame
        Manual adjusted and interpolated load time-series with UTC
        timestamps x ISO-2 countries
    """

    if 'ME' in load:
        if 'AL' not in load and 'AL' in countries:
            load['AL'] = load.ME * (5.7/2.9)
        if 'MK' not in load and 'MK' in countries:
            load['MK'] = load.ME * (6.7/2.9)

    if 'BG' in load:
        load = copy_timeslice(load, 'BG', '2018-10-27 21:00', '2018-10-28 22:00',
                              pd.Timedelta(weeks=1))

    return load


if __name__ == "__main__":

    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_load_data')

    configure_logging(snakemake)

    config = snakemake.config
    url = config['load']['url']
    interpolate_limit = config['load']['interpolate_limit']
    countries = config['countries']
    snapshots = pd.date_range(freq='h', **config['snapshots'])
    years = slice(snapshots[0], snapshots[-1])
    time_shift = config['load']['time_shift_for_large_gaps']

    load = load_timeseries(url, years, countries)

    if config['load']['manual_adjustments']:
        load = manual_adjustment(load)

    logger.info(f"Interpolate gaps of size {interpolate_limit} or less linearly.")
    load = load.interpolate(method='linear', limit=interpolate_limit)

    logger.info("Filling larger gaps by copying time-slices of period "
                f"'{time_shift}'.")
    load = load.apply(fill_large_gaps, shift=time_shift)

    assert not load.isna().any().any(), (
        'Load data contains nans. Adjust the parameters '
        '`time_shift_for_large_gaps` or modify the `manual_adjustment` function '
        'for implementing the needed load data modifications.')

    load.to_csv(snakemake.output[0])

