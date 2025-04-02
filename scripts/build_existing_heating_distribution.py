# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Builds table of existing heat generation capacities for initial planning
horizon.

Existing heat generation capacities are distributed to nodes based on population.
Within the nodes, the capacities are distributed to sectors (residential and services) based on sectoral consumption and urban/rural based population distribution.

Outputs:
--------
- Existing heat generation capacities distributed to nodes: `resources/{run_name}/existing_heating_distribution_base_s_{clusters}_{planning_horizons}.csv`


Notes
-----
- Data for Albania, Montenegro and Macedonia is not included in input database and assumed 0.
- Coal and oil boilers are assimilated to oil boilers.
- All ground-source heat pumps are assumed in rural areas and all air-source heat pumps are assumed to be in urban areas.

References
----------
- "Mapping and analyses of the current and future (2020 - 2030) heating/cooling fuel deployment (fossil/renewables)" (https://energy.ec.europa.eu/publications/mapping-and-analyses-current-and-future-2020-2030-heatingcooling-fuel-deployment-fossilrenewables-1_en)
"""

import logging

import country_converter as coco
import numpy as np
import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()


def build_existing_heating():
    """
    Retrieve and clean existing heating capacities for the myopic code.
    Data comes from the study "Mapping and analyses of the current and
    future (2020 - 2030) heating/cooling fuel deployment (fossil/renewables)".

    Source
    ------
    https://energy.ec.europa.eu/publications/mapping-and-analyses-current-and-future-2020-2030-heatingcooling-fuel-deployment-fossilrenewables-1_en

    File
    ----
    "WP2_DataAnnex_1_BuildingTechs_ForPublication_201603.xls" -> "existing_heating_raw.csv".

    Notes
    -----
    Data is for buildings only (i.e. NOT district heating) and represents the year 2012.
    """
    # TODO start from original file

    existing_heating = pd.read_csv(
        snakemake.input.existing_heating, index_col=0, header=0
    )

    # data for Albania, Montenegro and Macedonia not included in database
    existing_heating.loc["Albania"] = np.nan
    existing_heating.loc["Montenegro"] = np.nan
    existing_heating.loc["Macedonia"] = np.nan

    existing_heating.fillna(0.0, inplace=True)

    # convert GW to MW
    existing_heating *= 1e3

    existing_heating.index = cc.convert(existing_heating.index, to="iso2")

    # coal and oil boilers are assimilated to oil boilers
    existing_heating["oil boiler"] = (
        existing_heating["oil boiler"] + existing_heating["coal boiler"]
    )
    existing_heating.drop(["coal boiler"], axis=1, inplace=True)

    # distribute technologies to nodes by population
    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    nodal_heating = existing_heating.loc[pop_layout.ct]
    nodal_heating.index = pop_layout.index
    nodal_heating = nodal_heating.multiply(pop_layout.fraction, axis=0)

    district_heat_info = pd.read_csv(snakemake.input.district_heat_share, index_col=0)
    urban_fraction = district_heat_info["urban fraction"]

    energy_layout = pd.read_csv(
        snakemake.input.clustered_pop_energy_layout, index_col=0
    )

    uses = ["space", "water"]
    sectors = ["residential", "services"]

    nodal_sectoral_totals = pd.DataFrame(dtype=float)

    for sector in sectors:
        nodal_sectoral_totals[sector] = energy_layout[
            [f"total {sector} {use}" for use in uses]
        ].sum(axis=1)

    nodal_sectoral_fraction = nodal_sectoral_totals.div(
        nodal_sectoral_totals.sum(axis=1), axis=0
    )

    nodal_heat_name_fraction = pd.DataFrame(index=district_heat_info.index, dtype=float)

    nodal_heat_name_fraction["urban central"] = 0.0

    for sector in sectors:
        nodal_heat_name_fraction[f"{sector} rural"] = nodal_sectoral_fraction[
            sector
        ] * (1 - urban_fraction)
        nodal_heat_name_fraction[f"{sector} urban decentral"] = (
            nodal_sectoral_fraction[sector] * urban_fraction
        )

    nodal_heat_name_tech = pd.concat(
        {
            name: nodal_heating.multiply(nodal_heat_name_fraction[name], axis=0)
            for name in nodal_heat_name_fraction.columns
        },
        axis=1,
        names=["heat name", "technology"],
    )

    # move all ground HPs to rural, all air to urban

    for sector in sectors:
        nodal_heat_name_tech[(f"{sector} rural", "ground heat pump")] += (
            nodal_heat_name_tech[("urban central", "ground heat pump")]
            * nodal_sectoral_fraction[sector]
            + nodal_heat_name_tech[(f"{sector} urban decentral", "ground heat pump")]
        )
        nodal_heat_name_tech[(f"{sector} urban decentral", "ground heat pump")] = 0.0

        nodal_heat_name_tech[(f"{sector} urban decentral", "air heat pump")] += (
            nodal_heat_name_tech[(f"{sector} rural", "air heat pump")]
        )
        nodal_heat_name_tech[(f"{sector} rural", "air heat pump")] = 0.0

    # add large-scale heat pump sources as columns for district heating with 0 capacity

    for heat_pump_source in snakemake.params.sector["heat_pump_sources"][
        "urban central"
    ]:
        nodal_heat_name_tech[("urban central", f"{heat_pump_source} heat pump")] = 0.0

    nodal_heat_name_tech.to_csv(snakemake.output.existing_heating_distribution)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_existing_heating_distribution",
            clusters=48,
            planning_horizons=2050,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    build_existing_heating()
