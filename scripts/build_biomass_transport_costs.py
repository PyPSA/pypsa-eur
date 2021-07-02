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

ENERGY_CONTENT = 4.8  # unit MWh/tonne (assuming wood pellets)


def build_biomass_transport_costs():

    df_list = tbl.read_pdf(
        snakemake.input[0],
        pages="145-147",
        multiple_tables=True,
    )
    
    countries = df_list[0][0].iloc[6:].rename(index=lambda x: x + 1)

    # supply chain 1
    df = df_list[1].copy().rename(index=countries.to_dict())
    df.rename(
        columns=df.iloc[:6].apply(lambda col: col.str.cat(sep=" "), axis=0).to_dict(),
        inplace=True,
    )
    df = df.iloc[6:]
    df.loc[6] = df.loc[6].str.replace("€", "EUR")

    # supply chain 2
    df2 = df_list[2].copy().rename(index=countries.to_dict())
    df2.rename(
        columns=df2.iloc[:6].apply(lambda col: col.str.cat(sep=" "), axis=0).to_dict(),
        inplace=True,
    )
    df2 = df2.iloc[6:]
    df2.loc[6] = df2.loc[6].str.replace("€", "EUR")

    df.to_csv(snakemake.output.supply_chain1)
    df2.to_csv(snakemake.output.supply_chain1)

    transport_costs = pd.concat([df["per km/ton"], df2["per km/ton"]], axis=1).drop(6)
    transport_costs = transport_costs.astype(float, errors="ignore").mean(axis=1)

    # convert unit to EUR/MWh
    transport_costs /= ENERGY_CONTENT
    transport_costs = pd.DataFrame(transport_costs, columns=["cost [EUR/(km MWh)]"])

    # rename
    transport_costs.rename({"UK": "GB", "XK": "KO", "EL": "GR"}, inplace=True)

    # add missing Norway
    transport_costs.loc["NO"] = transport_costs.loc["SE"]
    transport_costs.to_csv(snakemake.output.transport_costs)


if __name__ == "__main__":

    prepare_biomass_transport_costs()
