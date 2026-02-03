# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve electricity demand for United Kingdom from NESO.

https://www.neso.energy/data-portal/historic-demand-data
https://www.neso.energy/data-portal/api-guidance
"""

import logging

import pandas as pd
import requests
from tqdm import tqdm

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

# Resource IDs for yearly demand data
RESOURCE_IDS = [
    "e8608e9a-f56c-457f-b9e7-bfffcfd19731",  # 2001
    "4daaac31-ae56-461e-8efe-6b33e147225a",  # 2002
    "30650247-acbf-414d-8ae2-2aa4aac57537",  # 2003
    "5f5c076c-d3a3-4f33-a137-501cc3d2ca9b",  # 2004
    "42a41acd-53b6-450d-af1b-4a92a352895c",  # 2005
    "949bb4a4-8374-4730-89ef-302d82428d2c",  # 2006
    "c0da767c-447d-48d1-93fd-ca3955d078ae",  # 2007
    "ec63adc8-e8a1-45cc-a593-5953b1ede7b1",  # 2008
    "ed8a37cb-65ac-4581-8dbc-a3130780da3a",  # 2009
    "b3eae4a5-8c3c-4df1-b9de-7db243ac3a09",  # 2010
    "01522076-2691-4140-bfb8-c62284752efd",  # 2011
    "4bf713a2-ea0c-44d3-a09a-63fc6a634b00",  # 2012
    "2ff7aaff-8b42-4c1b-b234-9446573a1e27",  # 2013
    "b9005225-49d3-40d1-921c-03ee2d83a2ff",  # 2014
    "cc505e45-65ae-4819-9b90-1fbb06880293",  # 2015
    "3bb75a28-ab44-4a0b-9b1c-9be9715d3c44",  # 2016
    "2f0f75b8-39c5-46ff-a914-ae38088ed022",  # 2017
    "fcb12133-0db0-4f27-a4a5-1669fd9f6d33",  # 2018
    "dd9de980-d724-415a-b344-d8ae11321432",  # 2019
    "33ba6857-2a55-479f-9308-e5c4c53d4381",  # 2020
    "18c69c42-f20d-46f0-84e9-e279045befc6",  # 2021
    "bb44a1b5-75b1-4db2-8491-257f23385006",  # 2022
    "bf5ab335-9b40-4ea4-b93a-ab4af7bce003",  # 2023
    "f6d02c0f-957b-48cb-82ee-09003f2ba759",  # 2024
    "b2bde559-3455-4021-b179-dfe60c0337b0",  # 2025
    "8a4a771c-3929-4e56-93ad-cdf13219dea5",  # 2026
]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_electricity_demand_neso")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    dfs = []
    for resource_id in tqdm(
        RESOURCE_IDS,
        ascii=False,
        unit="year",
        desc="Retrieving UK electricity demand from NESO",
    ):
        sql = (
            f'SELECT "SETTLEMENT_DATE", "SETTLEMENT_PERIOD", "ND" FROM "{resource_id}"'
        )
        r = requests.get(
            "https://api.neso.energy/api/3/action/datastore_search_sql",
            params={"sql": sql},
            timeout=90,
        )
        r.raise_for_status()

        data = r.json()["result"]["records"]
        df = pd.DataFrame(data)
        df.index = pd.date_range(
            start=df["SETTLEMENT_DATE"].iloc[0],
            periods=len(df),
            freq="30min",
            tz="UTC",
            inclusive="left",
        )
        df = df["ND"].astype(float).rename("GB")
        dfs.append(df)

    df = pd.concat(dfs)
    df = df.resample("1h").mean()
    df.to_csv(snakemake.output.csv)
