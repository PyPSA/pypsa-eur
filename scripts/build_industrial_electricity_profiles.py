# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build hourly industrial electricity demand profiles per model region for a planning horizon in MW.

The electricity demand per subsector and region is the specific electricity consumption from [build_industry_sector_ratios_intermediate][] (MWh/t) multiplied by the regional production from [build_industrial_production_per_node][]. The hourly electricity profile of a region is the sum of normalised profiles per industry branch, weighted by the region's electricity demand per subsector. The profiles are synthetic profiles for Germany in 2017 from FfE and are applied to all countries. They are rotated to match the weekdays of the modelled year, follow the public holidays and the local time of each country, and keep the annual mean load of a flat profile. FfE provides no profiles for chemicals and non-ferrous metals. These subsectors use the almost flat profiles of paper and steel, because their processes mostly run continuously at full load.

References
----------
- Ganz et al. (2021), [Wie können europäische Branchen-Lastgänge die Energiewende im Industriesektor unterstützen?](https://www.ffe.de/wp-content/uploads/2021/11/Wie-koennen-europaeische-Branchen-Lastgaenge-die-Energiewende-im-Industriesektor-unterstuetzen.pdf)
- FfE (2021), [Normalized industrial electrical load profiles (Germany)](https://opendata.ffe.de/dataset/normalized-industrial-electrical-load-profiles-germany/)
"""

import json
import logging

import holidays
import numpy as np
import pandas as pd
import pytz

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

FFE_YEAR = 2017

FFE_PROFILES = {
    1: "Iron & steel industry",
    4: "Non-metallic Minerals",
    5: "Transport Equipment",
    6: "Machinery",
    7: "Mining and Quarrying",
    8: "Food and Tobacco",
    9: "Paper, Pulp and Print",
    10: "Wood and Wood Products",
    11: "Construction",
    12: "Textile and Leather",
    13: "Non-specified (Industry)",
}

INDUSTRY_CATEGORY_TO_PROFILE = {
    "Electric arc": "Iron & steel industry",
    "DRI + Electric arc": "Iron & steel industry",
    "Integrated steelworks": "Iron & steel industry",
    "HVC": "Paper, Pulp and Print",
    "HVC (mechanical recycling)": "Paper, Pulp and Print",
    "HVC (chemical recycling)": "Paper, Pulp and Print",
    "Ammonia": "Paper, Pulp and Print",
    "Chlorine": "Paper, Pulp and Print",
    "Methanol": "Paper, Pulp and Print",
    "Other chemicals": "Paper, Pulp and Print",
    "Pharmaceutical products etc.": "Food and Tobacco",
    "Cement": "Non-metallic Minerals",
    "Ceramics & other NMM": "Non-metallic Minerals",
    "Glass production": "Non-metallic Minerals",
    "Pulp production": "Paper, Pulp and Print",
    "Paper production": "Paper, Pulp and Print",
    "Printing and media reproduction": "Paper, Pulp and Print",
    "Food, beverages and tobacco": "Food and Tobacco",
    "Alumina production": "Iron & steel industry",
    "Aluminium - primary production": "Iron & steel industry",
    "Aluminium - secondary production": "Iron & steel industry",
    "Other non-ferrous metals": "Iron & steel industry",
    "Transport equipment": "Transport Equipment",
    "Machinery equipment": "Machinery",
    "Textiles and leather": "Textile and Leather",
    "Wood and wood products": "Wood and Wood Products",
    "Other industrial sectors": "Non-specified (Industry)",
}


def load_ffe_load_profiles(json_file: str) -> pd.DataFrame:
    """
    Load the normalised FfE load profiles per industry branch.

    Returns
    -------
    pd.DataFrame
        Hourly profiles in German local time with one column per branch.
    """
    with open(json_file) as f:
        data = json.load(f)["data"]

    index = pd.date_range(f"{FFE_YEAR}-01-01", f"{FFE_YEAR}-12-31 23:00", freq="h")
    profiles = {FFE_PROFILES[d["internal_id"][0]]: d["values"] for d in data}

    return pd.DataFrame(profiles, index=index)


def shift_profiles(profiles: pd.DataFrame, year: int, country: str) -> pd.DataFrame:
    """
    Shift the FfE profiles to the weekdays and holidays of another year and country.

    German holidays in the profiles are replaced by the average day of the same
    weekday. The profiles are then rotated to align the weekdays, and the
    holidays of the country are replaced by the average Sunday.

    Parameters
    ----------
    profiles : pd.DataFrame
        Hourly profiles for Germany in the year of the FfE data.
    year : int
        Year to shift the profiles to.
    country : str
        Two-letter country code.

    Returns
    -------
    pd.DataFrame
        Hourly profiles in local time for the whole year, relative to their mean.
    """
    df = profiles.copy()
    average_days = df.groupby([df.index.dayofweek, df.index.hour]).mean()

    def average_day(index: pd.DatetimeIndex, sunday: bool = False) -> pd.DataFrame:
        dayofweek = [6] * len(index) if sunday else index.dayofweek
        return average_days.loc[list(zip(dayofweek, index.hour))].set_axis(index)

    def holiday_hours(index: pd.DatetimeIndex, country: str) -> pd.DatetimeIndex:
        calendar = holidays.country_holidays(country)
        return index[[hour in calendar for hour in index]]

    source_holidays = holiday_hours(df.index, "DE")
    df.loc[source_holidays] = average_day(source_holidays)

    hours = pd.date_range(f"{year}-01-01", f"{year}-12-31 23:00", freq="h")
    day_drift = df.index[0].dayofweek - hours[0].dayofweek
    df = pd.DataFrame(
        np.roll(df, 24 * day_drift, axis=0), index=hours[: len(df)], columns=df.columns
    )
    # leap years end with an additional average day
    df = pd.concat([df, average_day(hours[len(df) :])])

    target_holidays = holiday_hours(df.index, country)
    df.loc[target_holidays] = average_day(target_holidays, sunday=True)

    return df / df.mean()


def build_electricity_profiles(
    demand: pd.DataFrame, ffe_profiles: pd.DataFrame, snapshots: pd.DatetimeIndex
) -> pd.DataFrame:
    """
    Build hourly industrial electricity demand profiles per node.

    Parameters
    ----------
    demand : pd.DataFrame
        Annual electricity demand in TWh/a per subsector (index) and node (columns).
    ffe_profiles : pd.DataFrame
        Hourly profiles per industry branch from `load_ffe_load_profiles`.
    snapshots : pd.DatetimeIndex
        Snapshots in UTC.

    Returns
    -------
    pd.DataFrame
        Electricity demand in MW per snapshot (index) and node (columns).
    """
    demand = demand.groupby(INDUSTRY_CATEGORY_TO_PROFILE.__getitem__).sum()
    nodes_by_country = demand.columns.groupby(demand.columns.str[:2])

    profiles = []
    for country, nodes in nodes_by_country.items():
        timezone = pytz.country_timezones["RS" if country == "XK" else country][0]
        local = snapshots.tz_localize("UTC").tz_convert(timezone).tz_localize(None)
        shifted = pd.concat(
            shift_profiles(ffe_profiles, year, country) for year in local.year.unique()
        )
        shifted = shifted.loc[local, demand.index].set_axis(snapshots)
        profiles.append(shifted @ demand[nodes])

    return pd.concat(profiles, axis=1)[demand.columns] * 1e6 / 8760


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_electricity_profiles", horizon=2030
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    fn = snakemake.input.industry_sector_ratios
    sector_ratios = pd.read_csv(fn, header=[0, 1], index_col=0).loc["elec"].unstack(0)

    # material demand per industry and node (Mt/a)
    fn = snakemake.input.industrial_production_per_node
    production = pd.read_csv(fn, index_col=0).T / 1e3

    # electricity demand per industry and node (TWh/a)
    countries = production.columns.str[:2]
    demand = production * sector_ratios[countries].set_axis(production.columns, axis=1)

    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    ffe_profiles = load_ffe_load_profiles(snakemake.input.ffe_profiles)

    profiles = build_electricity_profiles(demand, ffe_profiles, snapshots)

    profiles.to_csv(snakemake.output.industrial_electricity_profiles)
