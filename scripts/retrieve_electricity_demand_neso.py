# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve electricity demand for United Kingdom from NESO.

https://www.neso.energy/data-portal/historic-demand-data
"""

import logging

import pandas as pd
from tqdm import tqdm

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

URLS = [
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/e8608e9a-f56c-457f-b9e7-bfffcfd19731/download/demanddata_2001.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/4daaac31-ae56-461e-8efe-6b33e147225a/download/demanddata_2002.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/30650247-acbf-414d-8ae2-2aa4aac57537/download/demanddata_2003.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/5f5c076c-d3a3-4f33-a137-501cc3d2ca9b/download/demanddata_2004.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/42a41acd-53b6-450d-af1b-4a92a352895c/download/demanddata_2005.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/949bb4a4-8374-4730-89ef-302d82428d2c/download/demanddata_2006.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/c0da767c-447d-48d1-93fd-ca3955d078ae/download/demanddata_2007.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/ec63adc8-e8a1-45cc-a593-5953b1ede7b1/download/demanddata_2008.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/ed8a37cb-65ac-4581-8dbc-a3130780da3a/download/demanddata_2009.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/b3eae4a5-8c3c-4df1-b9de-7db243ac3a09/download/demanddata_2010.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/01522076-2691-4140-bfb8-c62284752efd/download/demanddata_2011.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/4bf713a2-ea0c-44d3-a09a-63fc6a634b00/download/demanddata_2012.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/2ff7aaff-8b42-4c1b-b234-9446573a1e27/download/demanddata_2013.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/b9005225-49d3-40d1-921c-03ee2d83a2ff/download/demanddata_2014.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/cc505e45-65ae-4819-9b90-1fbb06880293/download/demanddata_2015.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/3bb75a28-ab44-4a0b-9b1c-9be9715d3c44/download/demanddata_2016.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/2f0f75b8-39c5-46ff-a914-ae38088ed022/download/demanddata_2017.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/fcb12133-0db0-4f27-a4a5-1669fd9f6d33/download/demanddata_2018.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/dd9de980-d724-415a-b344-d8ae11321432/download/demanddata_2019.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/33ba6857-2a55-479f-9308-e5c4c53d4381/download/demanddata_2020.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/18c69c42-f20d-46f0-84e9-e279045befc6/download/demanddata_2021.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/bb44a1b5-75b1-4db2-8491-257f23385006/download/demanddata_2022.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/bf5ab335-9b40-4ea4-b93a-ab4af7bce003/download/demanddata_2023.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/f6d02c0f-957b-48cb-82ee-09003f2ba759/download/demanddata_2024.csv",
    "https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/b2bde559-3455-4021-b179-dfe60c0337b0/download/demanddata_2025.csv",
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

    tqdm_kwargs = dict(
        ascii=False,
        unit=" years",
        total=len(URLS),
        desc="Retrieving UK electricity demand from NESO",
    )
    dfs = []
    for url in tqdm(URLS, **tqdm_kwargs):
        df = pd.read_csv(url, usecols=["SETTLEMENT_DATE", "SETTLEMENT_PERIOD", "ND"])
        df.index = pd.date_range(
            start=df.SETTLEMENT_DATE.iloc[0],
            periods=len(df),
            freq="30min",
            tz="UTC",
            inclusive="left",
        )
        # for checking with settlement date and period column
        # df = df.tz_convert("Europe/London")
        df = df["ND"].rename("GB")
        dfs.append(df)

    df = pd.concat(dfs)
    df = df.resample("1h").mean()

    df.to_csv(snakemake.output.csv)
