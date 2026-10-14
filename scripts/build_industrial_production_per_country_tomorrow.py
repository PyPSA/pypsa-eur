# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build industrial production per country and subsector for a planning horizon in kt/a.

Starting from today's production from [build_industrial_production_per_country][], the steel, aluminium and high-value chemical (HVC) subsectors are re-split between primary and secondary routes according to the `industry` configuration for the horizon year; values given per year are interpolated. Total steel is preserved: a share `St_primary_fraction` stays primary, of which `DRI_fraction` moves from integrated steelworks to the new "DRI + Electric arc" route, and the rest becomes electric arc (scrap) steel. Aluminium is likewise split by `Al_primary_fraction`. HVC production is split into virgin production (`HVC_primary_fraction`) and the new subsectors "HVC (mechanical recycling)" and "HVC (chemical recycling)".

!!! note "Primary shares are European"
    The primary fractions apply to the European total. Each country's primary production is scaled by a common factor so that the sum across countries meets the target share, which preserves the existing spatial pattern of primary plants.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.prepare_sector_network import get

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_industrial_production_per_country_tomorrow")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params.industry

    investment_year = int(snakemake.wildcards.horizon)

    fn = snakemake.input.industrial_production_per_country
    production = pd.read_csv(fn, index_col=0)

    keys = ["Integrated steelworks", "Electric arc"]
    total_steel = production[keys].sum(axis=1)

    st_primary_fraction = get(params["St_primary_fraction"], investment_year)
    dri_fraction = get(params["DRI_fraction"], investment_year)
    int_steel = production["Integrated steelworks"].sum()
    fraction_persistent_primary = st_primary_fraction * total_steel.sum() / int_steel

    dri = (
        dri_fraction * fraction_persistent_primary * production["Integrated steelworks"]
    )
    production.insert(2, "DRI + Electric arc", dri)

    not_dri = 1 - dri_fraction
    production["Integrated steelworks"] = (
        not_dri * fraction_persistent_primary * production["Integrated steelworks"]
    )
    production["Electric arc"] = (
        total_steel
        - production["DRI + Electric arc"]
        - production["Integrated steelworks"]
    )

    keys = ["Aluminium - primary production", "Aluminium - secondary production"]
    total_aluminium = production[keys].sum(axis=1)

    key_pri = "Aluminium - primary production"
    key_sec = "Aluminium - secondary production"

    al_primary_fraction = get(params["Al_primary_fraction"], investment_year)
    fraction_persistent_primary = (
        al_primary_fraction * total_aluminium.sum() / (production[key_pri].sum() or 1)
    )

    production[key_pri] = fraction_persistent_primary * production[key_pri]
    production[key_sec] = total_aluminium - production[key_pri]

    production["HVC (mechanical recycling)"] = (
        get(params["HVC_mechanical_recycling_fraction"], investment_year)
        * production["HVC"]
    )
    production["HVC (chemical recycling)"] = (
        get(params["HVC_chemical_recycling_fraction"], investment_year)
        * production["HVC"]
    )

    production["HVC"] *= get(params["HVC_primary_fraction"], investment_year)

    fn = snakemake.output.industrial_production_per_country_tomorrow
    production.to_csv(fn, float_format="%.2f")
