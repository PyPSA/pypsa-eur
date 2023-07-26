#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build hydroelectric inflow time-series for each country.

Relevant Settings
-----------------

.. code:: yaml

    countries:

    renewable:
        hydro:
            cutout:
            clip_min_inflow:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`toplevel_cf`, :ref:`renewable_cf`

Inputs
------

- ``data/bundle/EIA_hydro_generation_2000_2014.csv``: Hydroelectricity net generation per country and year (`EIA <https://www.eia.gov/beta/international/data/browser/#/?pa=000000000000000000000000000000g&c=1028i008006gg6168g80a4k000e0ag00gg0004g800ho00g8&ct=0&ug=8&tl_id=2-A&vs=INTL.33-12-ALB-BKWH.A&cy=2014&vo=0&v=H&start=2000&end=2016>`_)

    .. image:: img/hydrogeneration.png
        :scale: 33 %

- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``"cutouts/" + config["renewable"]['hydro']['cutout']``: confer :ref:`cutout`

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

Description
-----------

.. seealso::
    :mod:`build_renewable_profiles`
"""

import logging

import atlite
import country_converter as coco
import geopandas as gpd
import pandas as pd
from _helpers import configure_logging
from numpy.polynomial import Polynomial

cc = coco.CountryConverter()


def get_eia_annual_hydro_generation(fn, countries, capacities=False):
    # in billion kWh/a = TWh/a
    # capacity: in billion kW = GW
    df = pd.read_csv(fn, skiprows=2, index_col=1, na_values=[" ", "--"]).iloc[1:, 1:]
    df.index = df.index.str.strip()

    former_countries = {
        "Former Czechoslovakia": dict(
            countries=["Czech Republic", "Slovakia"], start=1980, end=1992
        ),
        "Former Serbia and Montenegro": dict(
            countries=["Serbia", "Montenegro"], start=1992, end=2005
        ),
        "Former Yugoslavia": dict(
            countries=[
                "Slovenia",
                "Croatia",
                "Bosnia and Herzegovina",
                "Serbia",
                "Montenegro",
                "North Macedonia",
            ],
            start=1980,
            end=1991,
        ),
    }

    for k, v in former_countries.items():
        period = [str(i) for i in range(v["start"], v["end"] + 1)]
        ratio = df.loc[v["countries"]].T.dropna().sum()
        ratio /= ratio.sum()
        for country in v["countries"]:
            df.loc[country, period] = df.loc[k, period] * ratio[country]

    baltic_states = ["Latvia", "Estonia", "Lithuania"]
    df.loc[baltic_states] = (
        df.loc[baltic_states].T.fillna(df.loc[baltic_states].mean(axis=1)).T
    )

    df.loc["Germany"] = df.filter(like="Germany", axis=0).sum()
    df.loc["Serbia"] += df.loc["Kosovo"].fillna(0.0)
    df = df.loc[~df.index.str.contains("Former")]
    df.drop(["Europe", "Germany, West", "Germany, East", "Kosovo"], inplace=True)

    df.index = cc.convert(df.index, to="iso2")
    df.index.name = "countries"

    # convert to MW of MWh/a
    factor = 1e3 if capacities else 1e6
    df = df.T[countries] * factor

    return df


def correct_eia_stats_by_capacity(eia_stats, fn, countries, baseyear=2019):
    cap = get_eia_annual_hydro_generation(fn, countries, capacities=True)
    ratio = cap / cap.loc[str(baseyear)]
    eia_stats_corrected = eia_stats / ratio
    to_keep = ["AL", "AT", "CH", "DE", "GB", "NL", "RS", "RO", "SK"]
    to_correct = eia_stats_corrected.columns.difference(to_keep)
    eia_stats.loc[:, to_correct] = eia_stats_corrected.loc[:, to_correct]


def approximate_missing_eia_stats(eia_stats, runoff_fn, countries):
    runoff = pd.read_csv(runoff_fn, index_col=0).T[countries]

    # fix ES, PT data points
    runoff.loc["1978", ["ES", "PT"]] = runoff.loc["1979", ["ES", "PT"]]

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

        snakemake = mock_snakemake("build_hydro_profile", weather_year="")
    configure_logging(snakemake)

    params_hydro = snakemake.params.hydro
    cutout = atlite.Cutout(snakemake.input.cutout)

    countries = snakemake.params.countries
    country_shapes = (
        gpd.read_file(snakemake.input.country_shapes)
        .set_index("name")["geometry"]
        .reindex(countries)
    )
    country_shapes.index.name = "countries"

    fn = snakemake.input.eia_hydro_generation
    eia_stats = get_eia_annual_hydro_generation(fn, countries)

    if config_hydro.get("eia_correct_by_capacity"):
        fn = snakemake.input.eia_hydro_capacity
        correct_eia_stats_by_capacity(eia_stats, fn, countries)

    if config_hydro.get("eia_approximate_missing"):
        fn = snakemake.input.era5_runoff
        eia_stats = approximate_missing_eia_stats(eia_stats, fn, countries)

    eia_stats.to_csv(snakemake.output.eia_hydro)

    weather_year = snakemake.wildcards.weather_year
    norm_year = config_hydro.get("eia_norm_year")
    if norm_year:
        eia_stats.loc[weather_year] = eia_stats.loc[norm_year]
    elif weather_year and weather_year not in eia_stats.index:
        eia_stats.loc[weather_year] = eia_stats.median()

    inflow = cutout.runoff(
        shapes=country_shapes,
        smooth=True,
        lower_threshold_quantile=True,
        normalize_using_yearly=eia_stats,
    )

    if "clip_min_inflow" in params_hydro:
        inflow = inflow.where(inflow > params_hydro["clip_min_inflow"], 0)

    inflow.to_netcdf(snakemake.output.profile)
