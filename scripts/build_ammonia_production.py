

import pandas as pd

ammonia = pd.read_excel(snakemake.input.usgs,
                        sheet_name="T12",
                        skiprows=5,
                        header=0,
                        index_col=0,
                        skipfooter=19)

rename = {"Austriae" : "AT",
          "Bulgaria" : "BG",
          "Belgiume" : "BE",
          "Croatia" : "HR",
          "Czechia" : "CZ",
          "Estonia" : "EE",
          "Finland" : "FI",
          "France" : "FR",
          "Germany" : "DE",
          "Greece" : "GR",
          "Hungarye" : "HU",
          "Italye" : "IT",
          "Lithuania" : "LT",
          "Netherlands" : "NL",
          "Norwaye" : "NO",
          "Poland" : "PL",
          "Romania" : "RO",
          "Serbia" : "RS",
          "Slovakia" : "SK",
          "Spain" : "ES",
          "Switzerland" : "CH",
          "United Kingdom" : "GB",
}

ammonia = ammonia.rename(rename)

ammonia = ammonia.loc[rename.values(),[str(i) for i in range(2013,2018)]].astype(float)

#convert from ktonN to ktonNH3
ammonia = ammonia*17/14

ammonia.index.name = "ktonNH3/a"

ammonia.to_csv(snakemake.output.ammonia_production)
