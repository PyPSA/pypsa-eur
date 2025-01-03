# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import pandas as pd
from _helpers import configure_logging, set_scenario_config
from entsoe import EntsoePandasClient
from entsoe.exceptions import NoMatchingDataError

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_cross_border_flows")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    api_key = snakemake.config["private"]["keys"]["entsoe_api"]
    client = EntsoePandasClient(api_key=api_key)

    start = pd.Timestamp(snakemake.params.snapshots["start"], tz="Europe/Brussels")
    end = pd.Timestamp(snakemake.params.snapshots["end"], tz="Europe/Brussels")

    countries = snakemake.params.countries

    prices = []
    unavailable_countries = []

    for country in countries:
        country_code = country

        try:
            gen = client.query_day_ahead_prices(country, start=start, end=end)
            gen = gen.tz_localize(None).resample("1h").mean()
            gen = gen.loc[start.tz_localize(None) : end.tz_localize(None)]
            prices.append(gen)
        except NoMatchingDataError:
            unavailable_countries.append(country)

    if unavailable_countries:
        logger.warning(
            f"Historical electricity prices for countries {', '.join(unavailable_countries)} not available."
        )

    keys = [c for c in countries if c not in unavailable_countries]
    prices = pd.concat(prices, keys=keys, axis=1)
    prices.to_csv(snakemake.output[0])
