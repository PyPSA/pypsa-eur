#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Created on Mon Jul  3 11:19:54 2023.

@author: fabian
"""


import logging

import pandas as pd
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
}


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_historical_electricity_generation")


api_key = "aeff3346-a240-40df-bd12-692772b845d0"
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
        gen = gen.rename(columns=carrier_grouper).groupby(level=0, axis=1).sum()
        generation.append(gen)
    except NoMatchingDataError:
        unavailable_countries.append(country)


if unavailable_countries:
    logger.warning(
        f"Historical electricity production for countries {', '.join(unavailable_countries)} not available."
    )

keys = [c for c in countries if c not in unavailable_countries]
generation = pd.concat(generation, keys=keys, axis=1)
generation = generation.loc[start.tz_localize(None) : end.tz_localize(None)]
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
generation.to_csv(snakemake.output[0])
