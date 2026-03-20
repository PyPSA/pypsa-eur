# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Extract historic Swiss energy balances in TWh/a from Excel file.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def load_sheet(xlsx: pd.ExcelFile, sheet_name: str, index_col: str) -> pd.DataFrame:
    df = xlsx.parse(sheet_name, skiprows=5, index_col=index_col)
    year_columns = [col for col in df.columns if isinstance(col, int)]
    df = df[year_columns].dropna(how="all", axis=0)
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_swiss_energy_balances")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    xlsx = pd.ExcelFile(snakemake.input.xlsx)

    tab1 = load_sheet(xlsx, "Tabelle1", index_col="Verwendungszweck")

    tab17 = load_sheet(xlsx, "Tabelle17", index_col="Verwendungszweck")

    tab18 = load_sheet(xlsx, "Tabelle18", index_col="Verwendungszweck")

    tab26 = load_sheet(xlsx, "Tabelle26", index_col="Verwendungszweck")

    tab28 = load_sheet(xlsx, "Tabelle28", index_col="Verwendungszweck")

    tab35 = load_sheet(xlsx, "Tabelle35", index_col="Verkehrsträger")

    tab38 = load_sheet(xlsx, "Tabelle38", index_col="Verkehrsträger")

    data = {
        "total residential": tab17.loc["Total Endenergieverbrauch"],
        "total residential space": tab17.loc["Raumwärme"],
        "total residential water": tab17.loc["Warmwasser"],
        "total residential cooking": tab17.loc["Kochen / Geschirrspüler"],
        "electricity residential": tab18.loc["Total"],
        "electricity residential space": tab18.loc["Raumwärme"],
        "electricity residential water": tab18.loc["Warmwasser"],
        "electricity residential cooking": tab18.loc["Kochherde"],
        "total services": tab26.loc["Total Endenergie"],
        "total services space": tab26.loc["Raumwärme"],
        "total services water": tab26.loc["Warmwasser"],
        "total services cooking": tab26.loc["Prozesswärme"],
        "electricity services": tab28.loc["Total Elektrizität"],
        "electricity services space": tab28.loc["Raumwärme"],
        "electricity services water": tab28.loc["Warmwasser"],
        "electricity services cooking": tab28.loc["Prozesswärme"],
        "total rail": tab35.loc["Schiene"],
        "total road": tab35.loc["Strasse"],
        "electricity road": tab38.loc["Strasse"],
        "electricity rail": tab38.loc["Schiene"],
        "total domestic aviation": tab35.loc["Luft (Inland)"],
        "total international aviation": tab1.loc["int. Flugverkehr"],
        "total domestic navigation": tab35.loc["Wasser"],
        "total international navigation": 0,
    }

    df = pd.DataFrame(data)
    df.columns.name = None
    df.index = pd.MultiIndex.from_product([["CH"], df.index], names=["country", "year"])

    # convert PJ/a to TWh/a
    df /= 3.6

    df.to_csv(snakemake.output.csv)
