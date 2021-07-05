"""Build ammonia production."""

import pandas as pd

country_to_alpha2 = {
    "Austriae": "AT",
    "Bulgaria": "BG",
    "Belgiume": "BE",
    "Croatia": "HR",
    "Czechia": "CZ",
    "Estonia": "EE",
    "Finland": "FI",
    "France": "FR",
    "Germany": "DE",
    "Greece": "GR",
    "Hungarye": "HU",
    "Italye": "IT",
    "Lithuania": "LT",
    "Netherlands": "NL",
    "Norwaye": "NO",
    "Poland": "PL",
    "Romania": "RO",
    "Serbia": "RS",
    "Slovakia": "SK",
    "Spain": "ES",
    "Switzerland": "CH",
    "United Kingdom": "GB",
}

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('build_ammonia_production')

    ammonia = pd.read_excel(snakemake.input.usgs,
                            sheet_name="T12",
                            skiprows=5,
                            header=0,
                            index_col=0,
                            skipfooter=19)

    ammonia.rename(country_to_alpha2, inplace=True)

    years = [str(i) for i in range(2013, 2018)]
    countries = country_to_alpha2.values()
    ammonia = ammonia.loc[countries, years].astype(float)

    # convert from ktonN to ktonNH3
    ammonia *= 17 / 14

    ammonia.index.name = "ktonNH3/a"

    ammonia.to_csv(snakemake.output.ammonia_production)
