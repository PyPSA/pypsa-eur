"""
Reads biomass transport costs for different countries of the JRC report

    "The JRC-EU-TIMES model.
    Bioenergy potentials
    for EU and neighbouring countries."
    (2015)

converts them from units 'EUR per km/ton' -> 'EUR/ (km MWh)'

assuming as an approximation energy content of wood pellets

@author: bw0928
"""

import pandas as pd
import tabula as tbl

ENERGY_CONTENT = 4.8  # unit MWh/t (wood pellets)

def get_countries():

    pandas_options = dict(
        skiprows=range(6),
        header=None,
        index_col=0
    )

    return tbl.read_pdf(
        str(snakemake.input.transport_cost_data),
        pages="145",
        multiple_tables=False,
        pandas_options=pandas_options
    )[0].index


def get_cost_per_tkm(page, countries):
    
    pandas_options = dict(
        skiprows=range(6),
        header=0,
        sep=' |,',
        engine='python',
        index_col=False,
    )

    sc = tbl.read_pdf(
        str(snakemake.input.transport_cost_data),
        pages=page,
        multiple_tables=False,
        pandas_options=pandas_options
    )[0]
    sc.index = countries
    sc.columns = sc.columns.str.replace("€", "EUR")
    
    return sc


def build_biomass_transport_costs():

    countries = get_countries()

    sc1 = get_cost_per_tkm(146, countries)
    sc2 = get_cost_per_tkm(147, countries)

    # take mean of both supply chains
    to_concat = [sc1["EUR/km/ton"], sc2["EUR/km/ton"]]
    transport_costs = pd.concat(to_concat, axis=1).mean(axis=1)

    # convert tonnes to MWh
    transport_costs /= ENERGY_CONTENT
    transport_costs.name = "EUR/km/MWh"

    # rename country names
    to_rename = {
        "UK": "GB",
        "XK": "KO",
        "EL": "GR"
    }
    transport_costs.rename(to_rename, inplace=True)

    # add missing Norway with data from Sweden
    transport_costs["NO"] = transport_costs["SE"]

    transport_costs.to_csv(snakemake.output[0])


if __name__ == "__main__":

    build_biomass_transport_costs()
