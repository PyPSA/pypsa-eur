# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

idx = pd.IndexSlice

# maps processes to temperature bands distributions from Fleiter et al. 2025
process_temperature_band_mapper = {
    "Aluminium - secondary production": idx[
        "Non-ferrous metals", "NE metals", ["Aluminium, secondary"]
    ],
    "Aluminium - primary production": idx[
        "Non-ferrous metals", "NE metals", ["Aluminium, primary"]
    ],
    "Other non-ferrous metals": idx[
        "Non-ferrous metals",
        "NE metals",
        [
            "Copper, primary",
            "Copper, secondary",
            "Copper further treatment",
            "Zinc, primary",
            "Zinc, secondary",
        ],
    ],
    "HVC": idx[
        :,
        :,
        [
            "Adipic acid",
            "Calcium carbide",
            "Carbon black",
            "Nitric acid",
            "Soda ash",
            "TDI",
            "Titanium dioxide",
        ],
    ],
    "Cement": idx["Non-metallic mineral products", "Clinker", :],
    "Ceramics & other NMM": idx["Non-metallic mineral products", "Ceramics", :],
    "Glass production": idx["Non-metallic mineral products", "Glass", :],
    "Pulp production": idx[
        "Paper and printing", "Paper products", ["Paper", "Paper Electrical"]
    ],
    "Paper production": idx[
        "Paper and printing", "Paper products", ["Chemical pulp", "Mechanical pulp"]
    ],
    "Printing and media reproduction": idx[
        "Paper and printing", "Paper products", ["Chemical pulp", "Mechanical pulp"]
    ],
    "Food, beverages and tobacco": idx["Food, drink and tobacco", :, :],
}

backup_temperature_band_shares = {
    # from JRC IDEES EU26 Industry
    "Alumina production": {
        "heat<100": 0.0,
        "heat100-200": 0.0,
        "heat200-500": 0.0,
        "heat>500": 1.0,
    },
    "Pharmaceutical products etc.": {
        "heat<100": 0.3,
        "heat100-200": 0.6,
        "heat200-500": 0.1,
        "heat>500": 1.0,
    },
    # from https://energyinnovation.org/data-explorer/overcoming-all-barriers-to-industrial-electrification
    "Other industrial sectors": {
        "heat<100": 0.1,
        "heat100-200": 0.65,
        "heat200-500": 0.25,
        "heat>500": 0.0,
    },
    "Transport equipment": {
        "heat<100": 0.25,
        "heat100-200": 0.6,
        "heat200-500": 0.15,
        "heat>500": 0.0,
    },
    "Machinery equipment": {
        "heat<100": 0.25,
        "heat100-200": 0.6,
        "heat200-500": 0.15,
        "heat>500": 0.0,
    },
    # from https://www.fpl.fs.usda.gov/documents/fplgtr/fplgtr118.pdf
    "Wood and wood products": {
        "heat<100": 0.65,
        "heat100-200": 0.35,
        "heat200-500": 0.0,
        "heat>500": 0.0,
    },
    # from https://www.ifc.org/content/dam/ifc/doc/2000/2007-textiles-manufacturing-ehs-guidelines-en.pdf
    "Textiles and leather": {
        "heat<100": 0.7,
        "heat100-200": 0.2,
        "heat200-500": 0.1,
        "heat>500": 0.0,
    },
}


def build_industry_sector_ratios_endogenous():
    # in TWh/a
    demand = pd.read_csv(
        snakemake.input.industrial_energy_demand_per_country_today,
        header=[0, 1],
        index_col=0,
    )

    # in Mt/a
    production = (
        pd.read_csv(snakemake.input.industrial_production_per_country, index_col=0)
        / 1e3
    ).stack()
    production.index.names = [None, None]

    today_sector_ratios = demand.div(production, axis=1).replace([np.inf, -np.inf], 0)

    today_sector_ratios.dropna(how="all", axis=1, inplace=True)

    rename = {
        "waste": "biomass",
        "electricity": "elec",
        "solid": "coke",
        "gas": "methane",
        "other": "biomass",
        "liquid": "naphtha",
    }
    today_sector_ratios = today_sector_ratios.rename(rename).groupby(level=0).sum()

    process_temperature_bands = pd.read_excel(
        snakemake.input.process_temperature_bands,
        sheet_name="Industry_ESC_FEC",
        index_col=[0, 1, 2],
        header=[0, 1],
    )
    process_temperature_bands.columns = process_temperature_bands.columns.droplevel(0)
    process_temperature_bands.rename(
        columns={
            "<\xa0100\xa0°C": "heat<100",
            "<\xa0100–200\xa0°C": "heat100-200",
            "200–500\xa0°C": "heat200-500",
            "500–1000\xa0°C": "heat500-1000",
            ">\xa01000\xa0°C": "heat>1000",
        },
        inplace=True,
    )
    process_temperature_bands.loc[:, "heat>500"] = (
        process_temperature_bands.loc[:, "heat500-1000"]
        + process_temperature_bands.loc[:, "heat>1000"]
    )
    process_temperature_bands = process_temperature_bands[
        ["heat<100", "heat100-200", "heat200-500", "heat>500"]
    ]

    heat_carriers = ["heat", "biomass", "methane"]
    bands = ["heat<100", "heat100-200", "heat200-500", "heat>500"]

    endogenous_sector_ratios = today_sector_ratios.copy()

    endogenous_sector_ratios = endogenous_sector_ratios.reindex(
        endogenous_sector_ratios.index.union(bands)
    ).replace(np.nan, 0)

    # for each heat-endogenous process, split energy in heat carriers into the respective temperature bands
    for process, key in process_temperature_band_mapper.items():
        as_fuels = endogenous_sector_ratios.loc[heat_carriers, idx[:, process]]

        as_heat = pd.DataFrame(
            np.outer(
                process_temperature_bands.loc[key].mean().values, as_fuels.sum().values
            ),
            index=bands,
            columns=as_fuels.sum().index,
        )

        endogenous_sector_ratios.loc[bands, idx[:, process]] = as_heat
        endogenous_sector_ratios.loc[heat_carriers, idx[:, process]] = 0.0

    for process, band_shares in backup_temperature_band_shares.items():
        band_shares = pd.Series(band_shares)

        as_fuels = endogenous_sector_ratios.loc[heat_carriers, idx[:, process]]

        as_heat = pd.DataFrame(
            np.outer(band_shares.values, as_fuels.sum().values),
            index=band_shares.index,
            columns=as_fuels.sum().index,
        )

        endogenous_sector_ratios.loc[bands, idx[:, process]] = as_heat
        endogenous_sector_ratios.loc[heat_carriers, idx[:, process]] = 0.0

    # heat is only nonzero in processes that are not endogenous, so this does not produce an error
    endogenous_sector_ratios.loc["heat<100"] += endogenous_sector_ratios.loc["heat"]
    endogenous_sector_ratios.drop("heat", inplace=True)

    endogenous_sector_ratios.index.name = "MWh/t"

    endogenous_sector_ratios.to_csv(snakemake.output.industry_sector_ratios_endogenous)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_sector_ratios_endogenous",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    build_industry_sector_ratios_endogenous()
