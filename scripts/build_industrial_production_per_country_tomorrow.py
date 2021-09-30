"""Build future industrial production per country."""

import pandas as pd

from prepare_sector_network import get

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('build_industrial_production_per_country_tomorrow')

    config = snakemake.config["industry"]

    investment_year = int(snakemake.wildcards.planning_horizons)

    fn = snakemake.input.industrial_production_per_country
    production = pd.read_csv(fn, index_col=0)

    keys = ["Integrated steelworks", "Electric arc"]
    total_steel = production[keys].sum(axis=1)

    st_primary_fraction = get(config["St_primary_fraction"], investment_year)
    dri_fraction = get(config["DRI_fraction"], investment_year)
    int_steel = production["Integrated steelworks"].sum()
    fraction_persistent_primary = st_primary_fraction * total_steel.sum() / int_steel

    dri = dri_fraction * fraction_persistent_primary * production["Integrated steelworks"]
    production.insert(2, "DRI + Electric arc", dri)

    not_dri = (1 - dri_fraction)
    production["Integrated steelworks"] = not_dri * fraction_persistent_primary * production["Integrated steelworks"]
    production["Electric arc"] = total_steel - production["DRI + Electric arc"] - production["Integrated steelworks"]

    keys = ["Aluminium - primary production", "Aluminium - secondary production"]
    total_aluminium = production[keys].sum(axis=1)

    key_pri = "Aluminium - primary production"
    key_sec = "Aluminium - secondary production"

    al_primary_fraction = get(config["Al_primary_fraction"], investment_year)
    fraction_persistent_primary = al_primary_fraction * total_aluminium.sum() / production[key_pri].sum()

    production[key_pri] = fraction_persistent_primary * production[key_pri]
    production[key_sec] = total_aluminium - production[key_pri]

    production["HVC (mechanical recycling)"] = get(config["HVC_mechanical_recycling_fraction"], investment_year) * production["HVC"]
    production["HVC (chemical recycling)"] = get(config["HVC_chemical_recycling_fraction"], investment_year) * production["HVC"]

    production["HVC"] *= get(config['HVC_primary_fraction'], investment_year)

    fn = snakemake.output.industrial_production_per_country_tomorrow
    production.to_csv(fn, float_format='%.2f')
