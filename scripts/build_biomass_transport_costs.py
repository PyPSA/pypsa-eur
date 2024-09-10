# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Reads biomass transport costs for different countries of the JRC report.

    "The JRC-EU-TIMES model.
    Bioenergy potentials
    for EU and neighbouring countries."
    (2015)

converts them from units 'EUR per km/ton' -> 'EUR/ (km MWh)'

assuming as an approximation energy content of wood pellets

@author: bw0928
"""
# try:
#     import tabula as tbl
#     import platform
# except:
#     ImportError("Please install tabula-py and platform")

# system = platform.system()
# encoding = "cp1252" if system == "Windows" else "utf-8"


import pandas as pd

ENERGY_CONTENT = 4.8  # unit MWh/t (wood pellets)


def get_countries():
    pandas_options = dict(
        skiprows=range(6), header=None, index_col=0, encoding=encoding
    )

    return tbl.read_pdf(
        str(snakemake.input.transport_cost_data),
        pages="145",
        multiple_tables=False,
        pandas_options=pandas_options,
        encoding=encoding,
    )[0].index


def get_cost_per_tkm(page, countries):
    pandas_options = dict(
        skiprows=range(6),
        header=0,
        sep=" |,",
        engine="python",
        index_col=False,
        encoding=encoding,
    )

    sc = tbl.read_pdf(
        str(snakemake.input.transport_cost_data),
        pages=page,
        multiple_tables=False,
        pandas_options=pandas_options,
        encoding=encoding,
    )[0]
    sc.index = countries
    sc.columns = sc.columns.str.replace("€", "EUR")

    return sc


def build_biomass_transport_costs():
    # Optional build from JRC report pdf, requires tabula and java dependencies.
    # countries = get_countries()
    # sc1 = get_cost_per_tkm(146, countries)
    # sc2 = get_cost_per_tkm(147, countries)

    # Use extracted csv from JRC report
    # https://publications.jrc.ec.europa.eu/repository/bitstream/JRC98626/biomass%20potentials%20in%20europe_web%20rev.pdf
    # Pages 146 (144) for supply chain 1 and 147 (145) for supply chain 2
    sc1 = pd.read_csv(snakemake.input.sc1, index_col=0, skiprows=2)
    sc2 = pd.read_csv(snakemake.input.sc2, index_col=0, skiprows=2)

    # take mean of both supply chains
    to_concat = [sc1["EUR/km/ton"], sc2["EUR/km/ton"]]
    transport_costs = pd.concat(to_concat, axis=1).mean(axis=1)

    # convert tonnes to MWh
    transport_costs /= ENERGY_CONTENT
    transport_costs.name = "EUR/km/MWh"

    # rename country names
    to_rename = {"UK": "GB", "EL": "GR"}
    transport_costs.rename(to_rename, inplace=True)

    # add missing Norway with data from Sweden
    transport_costs["NO"] = transport_costs["SE"]

    transport_costs.to_csv(snakemake.output[0])


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_biomass_transport_costs")

    build_biomass_transport_costs()
