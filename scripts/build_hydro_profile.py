# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build hydroelectric inflow time-series for each country.

Outputs
-------

- ``resources/profile_hydro.nc``:

    ===================  ================  =========================================================
    Field                Dimensions        Description
    ===================  ================  =========================================================
    inflow               countries, time   Inflow to the state of charge (in MW),
                                           e.g. due to river inflow in hydro reservoir.
    ===================  ================  =========================================================

    .. image:: img/inflow-ts.png
        :scale: 33 %

    .. image:: img/inflow-box.png
        :scale: 33 %
"""

import logging

import atlite
import country_converter as coco
import geopandas as gpd
import pandas as pd
from _helpers import configure_logging, get_snapshots, set_scenario_config
from numpy.polynomial import Polynomial

cc = coco.CountryConverter()


def get_eia_annual_hydro_generation(fn, countries, capacities=False):
    # in billion kWh/a = TWh/a
    df = pd.read_csv(fn, skiprows=2, index_col=1, na_values=[" ", "--"]).iloc[1:, 1:]
    df.index = df.index.str.strip()
    df.columns = df.columns.astype(int)

    former_countries = {
        "Former Czechoslovakia": dict(
            countries=["Czechia", "Slovakia"], start=1980, end=1992
        ),
        "Former Serbia and Montenegro": dict(
            countries=["Serbia", "Montenegro", "Kosovo"], start=1992, end=2005
        ),
        "Former Yugoslavia": dict(
            countries=[
                "Slovenia",
                "Croatia",
                "Bosnia and Herzegovina",
                "Serbia",
                "Kosovo",
                "Montenegro",
                "North Macedonia",
            ],
            start=1980,
            end=1991,
        ),
    }

    for k, v in former_countries.items():
        period = [i for i in range(v["start"], v["end"] + 1)]
        ratio = df.loc[v["countries"]].T.dropna().sum()
        ratio /= ratio.sum()
        for country in v["countries"]:
            df.loc[country, period] = df.loc[k, period] * ratio[country]

    baltic_states = ["Latvia", "Estonia", "Lithuania"]
    df.loc[baltic_states] = (
        df.loc[baltic_states].T.fillna(df.loc[baltic_states].mean(axis=1)).T
    )

    df.loc["Germany"] = df.filter(like="Germany", axis=0).sum()
    df = df.loc[~df.index.str.contains("Former")]
    df.drop(["Europe", "Germany, West", "Germany, East"], inplace=True)

    df.index = cc.convert(df.index, to="iso2")
    df.index.name = "countries"

    # convert to MW of MWh/a
    factor = 1e3 if capacities else 1e6
    df = df.T[countries] * factor

    df.ffill(axis=0, inplace=True)

    return df


def correct_eia_stats_by_capacity(eia_stats, fn, countries, baseyear=2019):
    cap = get_eia_annual_hydro_generation(fn, countries, capacities=True)
    ratio = cap / cap.loc[baseyear]
    eia_stats_corrected = eia_stats / ratio
    to_keep = ["AL", "AT", "CH", "DE", "GB", "NL", "RS", "XK", "RO", "SK"]
    to_correct = eia_stats_corrected.columns.difference(to_keep)
    eia_stats.loc[:, to_correct] = eia_stats_corrected.loc[:, to_correct]


def approximate_missing_eia_stats(eia_stats, runoff_fn, countries):
    runoff = pd.read_csv(runoff_fn, index_col=0).T[countries]
    runoff.index = runoff.index.astype(int)

    # fix outliers; exceptional floods in 1977-1979 in ES & PT
    if "ES" in runoff:
        runoff.loc[1978, "ES"] = runoff.loc[1979, "ES"]
    if "PT" in runoff:
        runoff.loc[1978, "PT"] = runoff.loc[1979, "PT"]

    runoff_eia = runoff.loc[eia_stats.index]

    eia_stats_approximated = {}

    for c in countries:
        X = runoff_eia[c]
        Y = eia_stats[c]

        to_predict = runoff.index.difference(eia_stats.index)
        X_pred = runoff.loc[to_predict, c]

        p = Polynomial.fit(X, Y, 1)
        Y_pred = p(X_pred)

        eia_stats_approximated[c] = pd.Series(Y_pred, index=to_predict)

    eia_stats_approximated = pd.DataFrame(eia_stats_approximated)
    return pd.concat([eia_stats, eia_stats_approximated]).sort_index()


logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_hydro_profile")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params_hydro = snakemake.params.hydro

    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)

    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=time)

    countries = snakemake.params.countries
    country_shapes = (
        gpd.read_file(snakemake.input.country_shapes)
        .set_index("name")["geometry"]
        .reindex(countries)
    )
    country_shapes.index.name = "countries"

    fn = snakemake.input.eia_hydro_generation
    eia_stats = get_eia_annual_hydro_generation(fn, countries)

    config_hydro = snakemake.config["renewable"]["hydro"]

    if config_hydro.get("eia_correct_by_capacity"):
        fn = snakemake.input.eia_hydro_capacity
        correct_eia_stats_by_capacity(eia_stats, fn, countries)

    if config_hydro.get("eia_approximate_missing"):
        fn = snakemake.input.era5_runoff
        eia_stats = approximate_missing_eia_stats(eia_stats, fn, countries)

    contained_years = pd.date_range(freq="YE", **snakemake.params.snapshots).year
    norm_year = config_hydro.get("eia_norm_year")
    missing_years = contained_years.difference(eia_stats.index)
    if norm_year:
        eia_stats.loc[contained_years] = eia_stats.loc[norm_year]
    elif missing_years.any():
        eia_stats.loc[missing_years] = eia_stats.median()

    inflow = cutout.runoff(
        shapes=country_shapes,
        smooth=True,
        lower_threshold_quantile=True,
        normalize_using_yearly=eia_stats,
    )

    if "clip_min_inflow" in params_hydro:
        inflow = inflow.where(inflow > params_hydro["clip_min_inflow"], 0)

    inflow.to_netcdf(snakemake.output.profile)
