# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import pandas as pd
from _helpers import configure_logging, set_scenario_config
from entsoe import EntsoePandasClient
from entsoe.exceptions import NoMatchingDataError

logger = logging.getLogger(__name__)


carrier_grouper = {
    "Waste": "Biomass",
    "Hydro Pumped Storage": "Hydro",
    "Hydro Water Reservoir": "Hydro",
    "Hydro Run-of-river and poundage": "Run of River",
    "Fossil Coal-derived gas": "Gas",
    "Fossil Gas": "Gas",
    "Fossil Oil": "Oil",
    "Fossil Oil shale": "Oil",
    "Fossil Brown coal/Lignite": "Lignite",
    "Fossil Peat": "Lignite",
    "Fossil Hard coal": "Coal",
    "Wind Onshore": "Onshore Wind",
    "Wind Offshore": "Offshore Wind",
    "Other renewable": "Other",
    "Marine": "Other",
}


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_electricity_production")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    api_key = snakemake.config["private"]["keys"]["entsoe_api"]
    client = EntsoePandasClient(api_key=api_key)

    start = pd.Timestamp(snakemake.params.snapshots["start"], tz="Europe/Brussels")
    end = pd.Timestamp(snakemake.params.snapshots["end"], tz="Europe/Brussels")

    countries = snakemake.params.countries

    generation = []
    unavailable_countries = []

    for country in countries:
        country_code = country

        try:
            gen = client.query_generation(country, start=start, end=end, nett=True)
            gen = gen.tz_localize(None).resample("1h").mean()
            gen = gen.loc[start.tz_localize(None) : end.tz_localize(None)]
            gen = gen.rename(columns=carrier_grouper).T.groupby(level=0).sum().T
            generation.append(gen)
        except NoMatchingDataError:
            unavailable_countries.append(country)

    if unavailable_countries:
        logger.warning(
            f"Historical electricity production for countries {', '.join(unavailable_countries)} not available."
        )

    keys = [c for c in countries if c not in unavailable_countries]
    generation = pd.concat(generation, keys=keys, axis=1)
    generation.to_csv(snakemake.output[0])
