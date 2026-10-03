# SPDX-FileCopyrightText: : 2025 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build historical CO2 or greenhouse gas emissions per country, year and sector from national inventories.

The national inventories reported to the UNFCCC are combined from three sources:
the EEA dataset for EU and EEA member countries, the EEA dataset for Energy
Community Contracting Parties, and the UK statistics for Great Britain. The IPCC
source categories are mapped to the sectors of PyPSA-Eur. Energy emissions not
assigned to another sector (e.g. refineries, manufacturing, fugitive emissions)
form "industrial non-elec". Serbia's inventory includes Kosovo, which is split by
population. Bosnia and Herzegovina and North Macedonia have no inventory and get
the per-capita emissions of Albania, Montenegro, Serbia and Kosovo.

Outputs
-------
- ``resources/<run_name>/co2_totals.csv``: Emissions in Mt per country, year and sector.

References
----------
- `EEA national emissions reported to the UNFCCC <https://doi.org/10.2909/a1407927-71d3-49da-9895-29a5e4cbcaf1>`_
- `EEA Energy Community greenhouse gas inventories <https://sdi.eea.europa.eu/catalogue/srv/api/records/b2f0b8a6-ae5f-44ff-ade4-7054a5cf5daa>`_
- `DESNZ final UK greenhouse gas emissions <https://www.gov.uk/government/collections/final-uk-greenhouse-gas-emissions-national-statistics>`_
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

to_ipcc = {
    "electricity": "1.A.1.a",
    "residential non-elec": "1.A.4.b",
    "services non-elec": "1.A.4.a",
    "rail non-elec": "1.A.3.c",
    "road non-elec": "1.A.3.b",
    "domestic navigation": "1.A.3.d",
    "international navigation": "1.D.1.b",
    "domestic aviation": "1.A.3.a",
    "international aviation": "1.D.1.a",
    "total energy": "1",
    "industrial processes": "2",
    "agriculture": "3",
    "agriculture, forestry and fishing": "1.A.4.c",
    "LULUCF": "4",
    "waste management": "5",
    "other": "6",
    "indirect": "ind_CO2",
}

# Population on 1 January 2019 from Eurostat (demo_gind)
POPULATION = {
    "AL": 2862427,
    "BA": 3492018,
    "ME": 622182,
    "MK": 2077132,
    "RS": 6963764,
    "XK": 1795666,
}


def read_eea(fn: str, scope: str) -> pd.DataFrame:
    """
    Read emissions in Mt by sector from an EEA inventory dataset.

    Parameters
    ----------
    fn : str
        Path to the EEA CSV file.
    scope : str
        Pollutant name, e.g. "CO2".

    Returns
    -------
    pd.DataFrame
        Emissions with a (country, year) index and sectors as columns.
    """
    df = pd.read_csv(fn, encoding="utf-8-sig", low_memory=False)
    df["Year"] = pd.to_numeric(df.Year, errors="coerce")
    return (
        df.query("Pollutant_name == @scope and Sector_code in @to_ipcc.values()")
        .dropna(subset="Year")
        .astype({"Year": int})
        .pivot_table(
            index=["Country_code", "Year"],
            columns="Sector_code",
            values="emissions",
            aggfunc="sum",
        )
        .rename_axis(index=["country", "year"], columns=None)
        .rename(index={"UK": "GB"}, columns={v: k for k, v in to_ipcc.items()})
        .div(1e3)
    )


def read_uk(fn: str, scope: str) -> pd.DataFrame:
    """
    Read UK emissions in Mt by sector from the DESNZ dataset of emissions by source.

    Parameters
    ----------
    fn : str
        Path to the DESNZ Excel file.
    scope : str
        "CO2" or "All greenhouse gases - (CO2 equivalent)".

    Returns
    -------
    pd.DataFrame
        Emissions with a (country, year) index and sectors as columns.
    """
    df = pd.read_excel(fn, sheet_name="UK_by_source")
    if scope == "CO2":
        df = df.query("GHG == 'CO2'")
    code = df["CRT category"].str.replace(".", "")

    def sector(c: str) -> pd.Series:
        # CRT category 1 excludes international bunkers (memo item 1.D)
        mask = code.str.startswith(c.replace(".", ""))
        if c == "1":
            mask &= ~code.str.startswith("1D")
        return df[mask].groupby("Year")["Emissions (MtCO2e)"].sum()

    uk = pd.DataFrame({k: sector(c) for k, c in to_ipcc.items()})
    uk.index = pd.MultiIndex.from_product([["GB"], uk.index], names=["country", "year"])
    return uk


def build_co2_totals(
    countries: list[str], eea: str, energy_community: str, uk: str, scope: str
) -> pd.DataFrame:
    """
    Combine national inventories and aggregate them to the sectors of PyPSA-Eur.

    Parameters
    ----------
    countries : list[str]
        ISO2 country codes to include.
    eea : str
        Path to the EEA dataset for EU and EEA member countries.
    energy_community : str
        Path to the EEA dataset for Energy Community Contracting Parties.
    uk : str
        Path to the DESNZ dataset of UK emissions by source.
    scope : str
        Pollutant name, e.g. "CO2" or "All greenhouse gases - (CO2 equivalent)".

    Returns
    -------
    pd.DataFrame
        Emissions in Mt with a (country, year) index and sectors as columns.
    """
    co2 = pd.concat(
        [read_eea(eea, scope), read_eea(energy_community, scope), read_uk(uk, scope)]
    ).fillna(0.0)

    rs = co2.loc["RS"] * POPULATION["RS"] / (POPULATION["RS"] + POPULATION["XK"])
    balkans = co2.loc[["AL", "ME", "RS"]].groupby("year").sum()
    per_capita = balkans / sum(POPULATION[c] for c in ["AL", "ME", "RS", "XK"])
    co2 = pd.concat(
        [
            co2.drop("RS", level="country"),
            pd.concat({"RS": rs, "XK": co2.loc["RS"] - rs}, names=["country"]),
            pd.concat(
                {c: per_capita * POPULATION[c] for c in ["BA", "MK"]}, names=["country"]
            ),
        ]
    )

    missing = pd.Index(countries).difference(co2.index.unique("country"))
    if not missing.empty:
        raise ValueError(f"No emissions inventory for countries {missing.tolist()}.")
    co2 = co2.loc[countries].sort_index()

    not_industry = [
        "electricity",
        "services non-elec",
        "residential non-elec",
        "road non-elec",
        "rail non-elec",
        "domestic aviation",
        "domestic navigation",
        "agriculture, forestry and fishing",
    ]
    co2["industrial non-elec"] = co2["total energy"] - co2[not_industry].sum(axis=1)
    co2["agriculture"] += co2["agriculture, forestry and fishing"]

    return co2.drop(columns=["total energy", "agriculture, forestry and fishing"])


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_co2_totals")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    co2 = build_co2_totals(
        snakemake.params.countries,
        snakemake.input.eea,
        snakemake.input.energy_community,
        snakemake.input.uk,
        snakemake.params.emissions_scope,
    )
    co2.to_csv(snakemake.output.co2_totals)
