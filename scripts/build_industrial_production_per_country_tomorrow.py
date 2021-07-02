"""Build future industrial production per country."""

import pandas as pd

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('build_industrial_production_per_country_tomorrow')

    config = snakemake.config["industry"]

    fn = snakemake.input.industrial_production_per_country
    production = pd.read_csv(fn, index_col=0)

    keys = ["Integrated steelworks", "Electric arc"]
    total_steel = production[keys].sum(axis=1)

    int_steel = production["Integrated steelworks"].sum()
    fraction_persistent_primary = config["St_primary_fraction"] * total_steel.sum() / int_steel

    dri = fraction_persistent_primary * production["Integrated steelworks"]
    production.insert(2, "DRI + Electric arc", dri)

    production["Electric arc"] = total_steel - production["DRI + Electric arc"]
    production["Integrated steelworks"] = 0.

    keys = ["Aluminium - primary production", "Aluminium - secondary production"]
    total_aluminium = production[keys].sum(axis=1)

    key_pri = "Aluminium - primary production"
    key_sec = "Aluminium - secondary production"
    fraction_persistent_primary = config["Al_primary_fraction"] * total_aluminium.sum() / production[key_pri].sum()
    production[key_pri] = fraction_persistent_primary * production[key_pri]
    production[key_sec] = total_aluminium - production[key_pri]

    production["Basic chemicals (without ammonia)"] *= config['HVC_primary_fraction']

    fn = snakemake.output.industrial_production_per_country_tomorrow
    production.to_csv(fn, float_format='%.2f')
