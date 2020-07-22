# coding: utf-8
"""

This rule downloads the load data from `Open Power System Data Time series <https://data.open-power-system-data.org/time_series/>` and interpolate the date to fill gaps. The resulting data saved as ``.csv`` file. The user can choose between two different data sources ("ENTSOE_transparency" or "ENTSOE_power_statistics"). 

Relevant Settings
-----------------

.. code:: yaml

    snapshots:

    load:
        url:
        source:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`load`

Inputs
------


Outputs
-------

- ``resource/load.csv``:


Description
-----------

The configuration options ``load: source:`` can be used to control which load data should be used in the model "ENTSOE_transparency" or "ENTSOE_power_statistics".


"""

import logging
logger = logging.getLogger(__name__)
from _helpers import configure_logging

import pandas as pd


def load_timeseries_opsd(years=None, fn=None, countries=None, source="ENTSOE_power_statistics"):
    """
    Read load data from OPSD time-series package version 2019-06-05.

    Parameters
    ----------
    years : None or slice()
        Years for which to read load data (defaults to
        slice("2018","2019"))
        
    fn : file name or url location (file format .csv)
    
    countries : Countries for which to read load data.
        
    source : "ENTSOE_transparency" or "ENTSOE_power_statistics"

    Returns
    -------
    load : pd.DataFrame
        Load time-series with UTC timestamps x ISO-2 countries
    """

    if countries is None:
        countries = snakemake.config['countries']
        
    logger.info(f"Retrieving load data from '{fn}'.")
     
    if source == 'ENTSOE_transparency':
        load = (pd.read_csv(fn, index_col=0, parse_dates=True)
                .filter(like='_load_actual_entsoe_transparency')
                .rename(columns=lambda s: s[:-len('_load_actual_entsoe_transparency')])
                .dropna(how="all", axis=0))
        
    elif source == 'ENTSOE_power_statistics':
        load = (pd.read_csv(fn, index_col=0, parse_dates=True)
            .filter(like='_load_actual_entsoe_power_statistics')
            .rename(columns=lambda s: s[:-len('_load_actual_entsoe_power_statistics')])
            .dropna(how="all", axis=0))
    else:
        raise NotImplementedError(f"Data for source `{source}` not available.")
    
    
    load = load.rename(columns={'GB_UKM' : 'GB'}).filter(items=countries)

    if years is not None:
        load = load.loc[years]
        
    
    return load

def consecutive_nans(ds):
    return (ds.isnull().astype(int)
            .groupby(ds.notnull().astype(int).cumsum()[ds.isnull()])
            .transform('sum').fillna(0))

def fill_large_gaps(ds, lower_bound=3, upper_bound=168):
    """Fill up large gaps (longer than 3 hours and up to one week (168 hours)) with load data from the previous week."""
    week_shift = pd.Series(ds.values, ds.index + pd.Timedelta('1w'))
    return ds.where((consecutive_nans(ds) < lower_bound) | (consecutive_nans(ds) > upper_bound), week_shift.reindex_like(ds))

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


def manual_adjustment(load, source="ENTSOE_power_statistics"):
    """
    Adjust gaps manual for load data from OPSD time-series package.

    Parameters
    ----------
    load : pd.DataFrame
        Load time-series with UTC timestamps x ISO-2 countries
        
    source : "ENTSOE_transparency" or "ENTSOE_power_statistics"

    Returns
    -------
    load : pd.DataFrame
        Manual adjusted and interpolated load time-series with UTC timestamps x ISO-2 countries
    """
    def copy_timeslice(load, cntry, start, stop, delta):
        start = pd.Timestamp(start)
        stop = pd.Timestamp(stop)
        if start in load.index and stop in load.index:
            load.loc[start:stop, cntry] = load.loc[start-delta:stop-delta, cntry].values
        return load
    
    # Kosovo (KV) and Albania (AL) do not exist in the data set 
    
    # Kosovo gets the same load curve as Serbia
    # scale parameter selected by energy consumption ratio from IEA Data browser
    # IEA 2012 -> scaling parameter = 4.8 / 27.
    # IEA 2017 -> scaling parameter = 5. / 33.
    # source https://www.iea.org/data-and-statistics?country=KOSOVO&fuel=Electricity%20and%20heat&indicator=Electricity%20final%20consumption
    
    if 'RS' in load.columns:
        if 'KV' not in load.columns or load.KV.isnull().values.all():
            load['KV'] = load['RS'] * (4.8 / 27.)
    
    # Albania gets the same load curve as Macedonia
    # scale parameter selected by energy consumption ratio from IEA Data browser
    # IEA 2012 -> scaling parameter = 4.1 / 7.4
    # IEA 2017 -> scaling parameter = 6.0 / 7.0
    # source https://www.iea.org/data-and-statistics?country=ALBANIA&fuel=Electricity%20and%20heat&indicator=Electricity%20final%20consumption

    if 'MK' in load.columns:
        if 'AL' not in load.columns or load.AL.isnull().values.all():
            load['AL'] = load['MK'] * (4.1 / 7.4)

    
    if source == 'ENTSOE_power_statistics':

        # To fill periods of load-gaps, we copy a period before into it
        load = copy_timeslice(load, 'GR', '2015-08-11 21:00', '2015-08-15 20:00', pd.Timedelta(weeks=1))
        load = copy_timeslice(load, 'AT', '2018-12-31 22:00', '2019-01-01 22:00', pd.Timedelta(days=2))
        load = copy_timeslice(load, 'CH', '2010-01-19 07:00', '2010-01-19 22:00', pd.Timedelta(days=1))
        load = copy_timeslice(load, 'CH', '2010-03-28 00:00', '2010-03-28 21:00', pd.Timedelta(days=1))
        load = copy_timeslice(load, 'CH', '2010-10-08 13:00', '2010-10-10 21:00', pd.Timedelta(weeks=1)) #is a WE, so take WE before
        load = copy_timeslice(load, 'CH', '2010-11-04 04:00', '2010-11-04 22:00', pd.Timedelta(days=1))
        load = copy_timeslice(load, 'NO', '2010-12-09 11:00', '2010-12-09 18:00', pd.Timedelta(days=1))
        load = copy_timeslice(load, 'GB', '2009-12-31 23:00', '2010-01-31 23:00', pd.Timedelta(days=-364)) #whole january missing
        

    if source == 'ENTSOE_transparency':
        
        # To fill periods of load-gaps, we copy a period before into it
        load = copy_timeslice(load, 'BG', '2018-10-27 21:00', '2018-10-28 22:00', pd.Timedelta(weeks=1))

    return load


if __name__ == "__main__":

    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_load_data')
        rootpath = '..'
    else:
        rootpath = '.'
    
    configure_logging(snakemake)
    
    url = snakemake.config['load']['url']
    
    opsd_load = (load_timeseries_opsd(years = slice(*pd.date_range(freq='y', **snakemake.config['snapshots'])[[0,-1]].year.astype(str)),
                                 fn=url,
                                 countries = snakemake.config['countries'],
                                 source = snakemake.config['load']['source']))

    # Convert to naive UTC (has to be explicit since pandas 0.24)
    opsd_load.index = opsd_load.index.tz_localize(None)
    
    # check the number and lenght of gaps
    nan_stats = nan_statistics(opsd_load)
    
    gap_filling_threshold = snakemake.config['load']['gap_filling_threshold']
      
    if nan_stats.consecutive.max() > gap_filling_threshold:        
        logger.warning(f"Load data contains consecutive gaps of longer than '{gap_filling_threshold}' hours! Check dataset carefully!")

    # adjust gaps manuel
    if snakemake.config['load']['adjust_gaps_manuel']:
        logger.info(f"Load data are adjusted manual. See adjustes in the manual_adjustment() function")
        opsd_load = manual_adjustment(load=opsd_load, source=snakemake.config['load']['source'])

    # Adjust gaps and interpolate load data with predefined heuristics
    logger.info(f"Gaps up to a length of one week and longer than {gap_filling_threshold} hours filled with data from previous week. Smaler gaps interpolated linearly.")
    # Longer gaps filled with data from previous week lower_bound is gap_filling_threshold and upper bound is one week.
    opsd_load = opsd_load.apply(fill_large_gaps, lower_bound=gap_filling_threshold)
    # Interpolate over the data with filling limit of gap_filling_threshold
    opsd_load = opsd_load.interpolate(method='linear', limit=gap_filling_threshold)

   
    # check the number and lenght of gaps after adjustment and interpolating
    nan_stats = nan_statistics(opsd_load)
    
    if nan_stats.consecutive.max() > 0:
        raise NotImplementedError(f'Load data contains gaps after adjustments. Modify manual_adjustment() function for implementing the needed load data modifications!')

    opsd_load.to_csv(snakemake.output[0])
    
