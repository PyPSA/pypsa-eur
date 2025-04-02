# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build total energy demands and carbon emissions per country using JRC IDEES,
eurostat, and EEA data.

- Country-specific data is read in :func:`build_eurostat`, :func:`build_idees` and `build_swiss`.
- :func:`build_energy_totals` then combines energy data from Eurostat, Swiss, and IDEES data and :func:`rescale_idees_from_eurostat` rescales IDEES data to match Eurostat data.
- :func:`build_district_heat_share` calculates the share of district heating for each country from IDEES data.
- Historical CO2 emissions are calculated in :func:`build_eea_co2` and :func:`build_eurostat_co2` and combined in :func:`build_co2_totals`.

Outputs
-------
- `resources/<run_name>/energy_totals.csv`: Energy totals per country, sector and year.
- `resources/<run_name>/co2_totals.csv`: CO2 emissions per country, sector and year.
- `resources/<run_name>/transport_data.csv`: Transport data per country and year.
- `resources/<run_name>/district_heat_share.csv`: District heating share per by country and year.
"""

import logging
import multiprocessing as mp
from functools import partial

import country_converter as coco
import geopandas as gpd
import numpy as np
import pandas as pd
from _helpers import configure_logging, mute_print, set_scenario_config
from tqdm import tqdm

cc = coco.CountryConverter()
logger = logging.getLogger(__name__)
idx = pd.IndexSlice


def cartesian(s1: pd.Series, s2: pd.Series) -> pd.DataFrame:
    """
    Compute the Cartesian product of two pandas Series.

    Parameters
    ----------
        s1: pd.Series
            The first pandas Series
        s2: pd.Series:
            The second pandas Series.

    Returns
    -------
    pd.DataFrame
        A DataFrame representing the Cartesian product of s1 and s2.

    Examples
    --------
    >>> s1 = pd.Series([1, 2, 3], index=["a", "b", "c"])
    >>> s2 = pd.Series([4, 5, 6], index=["d", "e", "f"])
    >>> cartesian(s1, s2)
       d  e   f
    a  4  5   6
    b  8 10  12
    c 12 15  18
    """
    return pd.DataFrame(np.outer(s1, s2), index=s1.index, columns=s2.index)


def reverse(dictionary: dict) -> dict:
    """
    Reverses the keys and values of a dictionary.

    Parameters
    ----------
    dictionary : dict
        The dictionary to be reversed.

    Returns
    -------
    dict
        A new dictionary with the keys and values reversed.

    Examples
    --------
    >>> d = {"a": 1, "b": 2, "c": 3}
    >>> reverse(d)
    {1: 'a', 2: 'b', 3: 'c'}
    """
    return {v: k for k, v in dictionary.items()}


idees_rename = {"GR": "EL", "GB": "UK"}

eu28 = cc.EU28as("ISO2").ISO2.tolist()
eu27 = cc.EU27as("ISO2").ISO2.tolist()
eu28_eea = eu28.copy()
eu28_eea.remove("GB")
eu28_eea.append("UK")


to_ipcc = {
    "electricity": "1.A.1.a - Public Electricity and Heat Production",
    "residential non-elec": "1.A.4.b - Residential",
    "services non-elec": "1.A.4.a - Commercial/Institutional",
    "rail non-elec": "1.A.3.c - Railways",
    "road non-elec": "1.A.3.b - Road Transportation",
    "domestic navigation": "1.A.3.d - Domestic Navigation",
    "international navigation": "1.D.1.b - International Navigation",
    "domestic aviation": "1.A.3.a - Domestic Aviation",
    "international aviation": "1.D.1.a - International Aviation",
    "total energy": "1 - Energy",
    "industrial processes": "2 - Industrial Processes and Product Use",
    "agriculture": "3 - Agriculture",
    "agriculture, forestry and fishing": "1.A.4.c - Agriculture/Forestry/Fishing",
    "LULUCF": "4 - Land Use, Land-Use Change and Forestry",
    "waste management": "5 - Waste management",
    "other": "6 - Other Sector",
    "indirect": "ind_CO2 - Indirect CO2",
    "total wL": "Total (with LULUCF)",
    "total woL": "Total (without LULUCF)",
}


def eurostat_per_country(input_eurostat: str, country: str) -> pd.DataFrame:
    """
    Read energy balance data for a specific country from Eurostat.

    Parameters
    ----------
    input_eurostat : str
        Path to the directory containing Eurostat data files.
    country : str
        Country code for the specific country.

    Returns
    -------
    pd.DataFrame
        Concatenated energy balance data for the specified country.

    Notes
    -----
    - The function reads `<input_eurostat>/<country>.-Energy-balance-sheets-April-2023-edition.xlsb`
    - It removes the "Cover" sheet from the data and concatenates all the remaining sheets into a single DataFrame.
    """

    filename = (
        f"{input_eurostat}/{country}-Energy-balance-sheets-April-2023-edition.xlsb"
    )
    sheet = pd.read_excel(
        filename,
        engine="pyxlsb",
        sheet_name=None,
        skiprows=4,
        index_col=list(range(4)),
        na_values=":",
    )
    sheet.pop("Cover")
    return pd.concat(sheet)


def build_eurostat(
    input_eurostat: str,
    countries: list[str],
    nprocesses: int = 1,
    disable_progressbar: bool = False,
) -> pd.DataFrame:
    """
    Return multi-index for all countries' energy data in TWh/a.

    Parameters
    ----------
    input_eurostat : str
        Path to the Eurostat database.
    countries : list[str]
        List of countries for which energy data is to be retrieved.
    nprocesses : int, optional
        Number of processes to use for parallel execution, by default 1.
    disable_progressbar : bool, optional
        Whether to disable the progress bar, by default False.

    Returns
    -------
    pd.DataFrame
        Multi-index DataFrame containing energy data for all countries in TWh/a.

    Notes
    -----
    - The function first renames the countries in the input list using the `idees_rename` mapping and removes "CH".
    - It then reads country-wise data using :func:`eurostat_per_country` into a single DataFrame.
    - The data is reordered, converted to TWh/a, and missing values are filled.
    """

    countries = {idees_rename.get(country, country) for country in countries} - {"CH"}

    func = partial(eurostat_per_country, input_eurostat)
    tqdm_kwargs = dict(
        ascii=False,
        unit=" country",
        total=len(countries),
        desc="Build from eurostat database",
        disable=disable_progressbar,
    )
    with mute_print():
        with mp.Pool(processes=nprocesses) as pool:
            dfs = list(tqdm(pool.imap(func, countries), **tqdm_kwargs))

    index_names = ["country", "year", "lvl1", "lvl2", "lvl3", "lvl4"]
    df = pd.concat(dfs, keys=countries, names=index_names)
    df.index = df.index.set_levels(df.index.levels[1].astype(int), level=1)

    # drop columns with all NaNs
    unnamed_cols = df.columns[df.columns.astype(str).str.startswith("Unnamed")]
    df.drop(unnamed_cols, axis=1, inplace=True)
    df.drop(list(range(1990, 2022)), axis=1, inplace=True, errors="ignore")

    # make numeric values where possible
    df.replace("Z", 0, inplace=True)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.select_dtypes(include=[np.number])

    # write 'International aviation' to the lower level of the multiindex
    int_avia = df.index.get_level_values(3) == "International aviation"
    temp = df.loc[int_avia]
    temp.index = pd.MultiIndex.from_frame(
        temp.index.to_frame().fillna("International aviation")
    )
    df = pd.concat([temp, df.loc[~int_avia]]).sort_index()

    # Fill in missing data on "Domestic aviation" for each country.
    for country in countries:
        slicer = idx[country, :, :, :, "Domestic aviation"]
        # For the Total and Fossil energy columns, fill in zeros with
        # the closest non-zero value in the year index.
        for col in ["Total", "Fossil energy"]:
            df.loc[slicer, col] = (
                df.loc[slicer, col].replace(0.0, np.nan).ffill().bfill()
            )

    # Renaming some indices
    index_rename = {
        "Households": "Residential",
        "Commercial & public services": "Services",
        "Domestic navigation": "Domestic Navigation",
        "International maritime bunkers": "Bunkers",
        "UK": "GB",
        "EL": "GR",
    }
    columns_rename = {"Total": "Total all products"}
    df.rename(index=index_rename, columns=columns_rename, inplace=True)
    df.sort_index(inplace=True)

    # convert to TWh/a from ktoe/a
    df *= 11.63 / 1e3

    return df


def build_swiss() -> pd.DataFrame:
    """
    Return a pd.DataFrame of Swiss energy data in TWh/a.

    Returns
    -------
    pd.DataFrame
        Swiss energy data in TWh/a.

    Notes
    -----
    - Reads Swiss energy data from `data/switzerland-new_format-all_years.csv`.
    - Reshapes and renames data.
    - Converts energy units from PJ/a to TWh/a.
    """
    fn = snakemake.input.swiss

    df = pd.read_csv(fn, index_col=[0, 1])

    df.columns = df.columns.astype(int)

    df.columns.name = "year"

    df = df.stack().unstack("item")

    df.columns.name = None

    # convert PJ/a to TWh/a
    df /= 3.6

    return df


def idees_per_country(ct: str, base_dir: str) -> pd.DataFrame:
    """
    Calculate energy totals per country using JRC-IDEES data.

    Parameters
    ----------
    ct : str
        The country code.
    base_dir : str
        The base directory where the JRC-IDEES data files are located.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the energy totals per country. Columns are energy uses.

    Notes
    -----
    - Retrieves JRC-IDEES data for the specified country from `base_dir` for residential, tertiary, and transport sectors.
    - Calculates energy totals for each sector, stores them in a dictionary and returns them as data frame.
    - Assertions ensure indices of JRC-IDEES data are as expected.
    """

    ct_idees = idees_rename.get(ct, ct)
    fn_residential = f"{base_dir}/{ct_idees}/JRC-IDEES-2021_Residential_{ct_idees}.xlsx"
    fn_tertiary = f"{base_dir}/{ct_idees}/JRC-IDEES-2021_Tertiary_{ct_idees}.xlsx"
    fn_transport = f"{base_dir}/{ct_idees}/JRC-IDEES-2021_Transport_{ct_idees}.xlsx"

    ct_totals = {}

    # residential

    df = pd.read_excel(fn_residential, "RES_hh_fec", index_col=0)

    rows = ["Advanced electric heating", "Conventional electric heating"]
    ct_totals["electricity residential space"] = df.loc[rows].sum()
    ct_totals["total residential space"] = df.loc["Space heating"]
    ct_totals["total residential water"] = df.loc["Water heating"]

    assert df.index[23] == "Electricity"
    ct_totals["electricity residential water"] = df.iloc[23]

    ct_totals["total residential cooking"] = df.loc["Cooking"]

    assert df.index[30] == "Electricity"
    ct_totals["electricity residential cooking"] = df.iloc[30]

    df = pd.read_excel(fn_residential, "RES_summary", index_col=0)

    row = "Energy consumption by fuel - Eurostat structure (ktoe)"
    ct_totals["total residential"] = df.loc[row]

    assert df.index[40] == "Electricity"
    ct_totals["electricity residential"] = df.iloc[40]

    assert df.index[39] == "Distributed heat"
    ct_totals["distributed heat residential"] = df.iloc[39]

    assert df.index[43] == "Thermal uses"
    ct_totals["thermal uses residential"] = df.iloc[43]

    df = pd.read_excel(fn_residential, "RES_hh_eff", index_col=0)

    ct_totals["total residential space efficiency"] = df.loc["Space heating"]

    assert df.index[5] == "Diesel oil"
    ct_totals["oil residential space efficiency"] = df.iloc[5]

    assert df.index[6] == "Natural gas"
    ct_totals["gas residential space efficiency"] = df.iloc[6]

    ct_totals["total residential water efficiency"] = df.loc["Water heating"]

    assert df.index[18] == "Diesel oil"
    ct_totals["oil residential water efficiency"] = df.iloc[18]

    assert df.index[19] == "Natural gas"
    ct_totals["gas residential water efficiency"] = df.iloc[19]

    # services

    df = pd.read_excel(fn_tertiary, "SER_hh_fec", index_col=0)

    ct_totals["total services space"] = df.loc["Space heating"]

    rows = ["Advanced electric heating", "Conventional electric heating"]
    ct_totals["electricity services space"] = df.loc[rows].sum()

    ct_totals["total services water"] = df.loc["Hot water"]

    assert df.index[24] == "Electricity"
    ct_totals["electricity services water"] = df.iloc[24]

    ct_totals["total services cooking"] = df.loc["Catering"]

    assert df.index[31] == "Electricity"
    ct_totals["electricity services cooking"] = df.iloc[31]

    df = pd.read_excel(fn_tertiary, "SER_summary", index_col=0)

    row = "Energy consumption by fuel - Eurostat structure (ktoe)"
    ct_totals["total services"] = df.loc[row]

    assert df.index[43] == "Electricity"
    ct_totals["electricity services"] = df.iloc[43]

    assert df.index[42] == "Distributed heat"
    ct_totals["distributed heat services"] = df.iloc[42]

    assert df.index[46] == "Thermal uses"
    ct_totals["thermal uses services"] = df.iloc[46]

    df = pd.read_excel(fn_tertiary, "SER_hh_eff", index_col=0)

    ct_totals["total services space efficiency"] = df.loc["Space heating"]

    assert df.index[5] == "Diesel oil"
    ct_totals["oil services space efficiency"] = df.iloc[5]

    assert df.index[7] == "Conventional gas heaters"
    ct_totals["gas services space efficiency"] = df.iloc[7]

    ct_totals["total services water efficiency"] = df.loc["Hot water"]

    assert df.index[20] == "Diesel oil"
    ct_totals["oil services water efficiency"] = df.iloc[20]

    assert df.index[21] == "Natural gas"
    ct_totals["gas services water efficiency"] = df.iloc[21]

    # agriculture, forestry and fishing

    start = "Detailed split of energy consumption (ktoe)"
    end = "Market shares of energy uses (%)"

    df = pd.read_excel(fn_tertiary, "AGR_fec", index_col=0).loc[start:end]

    rows = [
        "Lighting",
        "Ventilation",
        "Specific electricity uses",
        "Pumping devices (electricity)",
    ]
    ct_totals["total agriculture electricity"] = df.loc[rows].sum()

    rows = ["Specific heat uses", "Low enthalpy heat"]
    ct_totals["total agriculture heat"] = df.loc[rows].sum()

    rows = [
        "Motor drives",
        "Farming machine drives (diesel oil and liquid biofuels)",
        "Pumping devices (diesel oil and liquid biofuels)",
    ]
    ct_totals["total agriculture machinery"] = df.loc[rows].sum()

    row = "Agriculture, forestry and fishing"
    ct_totals["total agriculture"] = df.loc[row]

    # transport

    df = pd.read_excel(fn_transport, "TrRoad_ene", index_col=0)

    ct_totals["total road"] = df.loc["by fuel (EUROSTAT DATA)"]

    ct_totals["electricity road"] = df.loc["Electricity"]

    ct_totals["total two-wheel"] = df.loc["Powered two-wheelers (Gasoline)"]

    assert df.index[19] == "Passenger cars"
    ct_totals["total passenger cars"] = df.iloc[19]

    assert df.index[30] == "Battery electric vehicles"
    ct_totals["electricity passenger cars"] = df.iloc[30]

    assert df.index[31] == "Motor coaches, buses and trolley buses"
    ct_totals["total other road passenger"] = df.iloc[31]

    assert df.index[39] == "Battery electric vehicles"
    ct_totals["electricity other road passenger"] = df.iloc[39]

    assert df.index[41] == "Light commercial vehicles"
    ct_totals["total light duty road freight"] = df.iloc[41]

    assert df.index[49] == "Battery electric vehicles"
    ct_totals["electricity light duty road freight"] = df.iloc[49]

    row = "Heavy goods vehicles (Diesel oil incl. biofuels)"
    ct_totals["total heavy duty road freight"] = df.loc[row]

    assert df.index[61] == "Passenger cars"
    ct_totals["passenger car efficiency"] = df.iloc[61]

    df = pd.read_excel(fn_transport, "TrRail_ene", index_col=0)

    ct_totals["total rail"] = df.loc["by fuel"]

    ct_totals["electricity rail"] = df.loc["Electricity"]

    assert df.index[9] == "Passenger transport"
    ct_totals["total rail passenger"] = df.iloc[9]

    assert df.index[10] == "Metro and tram, urban light rail"
    assert df.index[13] == "Electric"
    assert df.index[14] == "High speed passenger trains"
    ct_totals["electricity rail passenger"] = df.iloc[[10, 13, 14]].sum()

    assert df.index[15] == "Freight transport"
    ct_totals["total rail freight"] = df.iloc[15]

    assert df.index[17] == "Electric"
    ct_totals["electricity rail freight"] = df.iloc[17]

    df = pd.read_excel(fn_transport, "TrAvia_ene", index_col=0)

    assert df.index[4] == "Passenger transport"
    ct_totals["total aviation passenger"] = df.iloc[4]

    assert df.index[8] == "Freight transport"
    ct_totals["total aviation freight"] = df.iloc[8]

    assert df.index[2] == "Domestic"
    ct_totals["total domestic aviation passenger"] = df.iloc[2]

    assert df.index[6] == "International - Intra-EEAwUK"
    assert df.index[7] == "International - Extra-EEAwUK"
    ct_totals["total international aviation passenger"] = df.iloc[[6, 7]].sum()

    assert df.index[9] == "Domestic"
    assert df.index[10] == "International - Intra-EEAwUK"
    ct_totals["total domestic aviation freight"] = df.iloc[[9, 10]].sum()

    assert df.index[11] == "International - Extra-EEAwUK"
    ct_totals["total international aviation freight"] = df.iloc[11]

    ct_totals["total domestic aviation"] = (
        ct_totals["total domestic aviation freight"]
        + ct_totals["total domestic aviation passenger"]
    )

    ct_totals["total international aviation"] = (
        ct_totals["total international aviation freight"]
        + ct_totals["total international aviation passenger"]
    )

    df = pd.read_excel(fn_transport, "TrNavi_ene", index_col=0)

    # coastal and inland
    ct_totals["total domestic navigation"] = df.loc["Energy consumption (ktoe)"]

    df = pd.read_excel(fn_transport, "TrRoad_act", index_col=0)

    assert df.index[85] == "Passenger cars"
    ct_totals["passenger cars"] = df.iloc[85]

    return pd.DataFrame(ct_totals)


def build_idees(countries: list[str]) -> pd.DataFrame:
    """
    Build energy totals from IDEES database for the given list of countries
    using :func:`idees_per_country`.

    Parameters
    ----------
    countries : list[str]
        List of country names for which energy totals need to be built.

    Returns
    -------
    pd.DataFrame
        Energy totals for the given countries.

    Notes
    -----
    - Retrieves energy totals per country and year using :func:`idees_per_country`.
    - Returns a DataFrame with columns: country, year, and energy totals for different categories.
    """

    nprocesses = snakemake.threads
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    func = partial(idees_per_country, base_dir=snakemake.input.idees)
    tqdm_kwargs = dict(
        ascii=False,
        unit=" country",
        total=len(countries),
        desc="Build from IDEES database",
        disable=disable_progress,
    )
    with mute_print():
        with mp.Pool(processes=nprocesses) as pool:
            totals_list = list(tqdm(pool.imap(func, countries), **tqdm_kwargs))

    totals = pd.concat(
        totals_list,
        keys=countries,
        names=["country", "year"],
    )

    # clean up dataframe
    years = np.arange(2000, 2022)
    totals = totals[totals.index.get_level_values(1).isin(years)]

    # efficiency kgoe/100km -> ktoe/100km so that after conversion TWh/100km
    totals.loc[:, "passenger car efficiency"] /= 1e6
    # convert ktoe to TWh
    patterns = ["passenger cars", ".*space efficiency", ".*water efficiency"]
    exclude = totals.columns.str.fullmatch("|".join(patterns))
    totals = totals.copy()
    totals.loc[:, ~exclude] *= 11.63 / 1e3

    return totals


def fill_missing_years(fill_values: pd.Series) -> pd.Series:
    """
    Fill missing years for some countries by first using forward fill (ffill)
    and then backward fill (bfill).

    Parameters
    ----------
    fill_values : pd.Series
        A pandas Series with a MultiIndex (levels: country and year) representing
        energy values, where some values may be zero and need to be filled.

    Returns
    -------
    pd.Series
        A pandas Series with zero values replaced by the forward-filled and
        backward-filled values of the corresponding country.

    Notes
    -----
    - The function groups the data by the 'country' level and performs forward fill
      and backward fill to fill zero values.
    - Zero values in the original Series are replaced by the ffilled and bfilled
      value of their respective country group.
    """

    # Forward fill and then backward fill within each country group
    fill_values = fill_values.groupby(level="country").ffill().bfill()

    return fill_values


def build_energy_totals(
    countries: list[str],
    eurostat: pd.DataFrame,
    swiss: pd.DataFrame,
    idees: pd.DataFrame,
) -> pd.DataFrame:
    """
    Combine energy totals for the specified countries from Eurostat, Swiss, and
    IDEES data.

    Parameters
    ----------
    countries : list[str]
        List of country codes for which energy totals are to be calculated.
    eurostat : pd.DataFrame
        Eurostat energy balances dataframe.
    swiss : pd.DataFrame
        Swiss energy data dataframe.
    idees : pd.DataFrame
        IDEES energy data dataframe.

    Returns
    -------
    pd.DataFrame
        Energy totals dataframe for the given countries.

    Notes
    -----
    - Missing values are filled based on Eurostat energy balances and average values in EU28.
    - The function also performs specific calculations for Norway and splits road, rail, and aviation traffic for non-IDEES data.

    References
    ----------
    - `Norway heating data <http://www.ssb.no/en/energi-og-industri/statistikker/husenergi/hvert-3-aar/2014-07-14>`_
    """

    eurostat_fuels = {"electricity": "Electricity", "total": "Total all products"}
    eurostat_countries = eurostat.index.unique(0)
    eurostat_years = eurostat.index.unique(1)

    new_index = pd.MultiIndex.from_product(
        [countries, eurostat_years], names=["country", "year"]
    )

    efficiency_keywords = ["space efficiency", "water efficiency"]
    to_drop = idees.columns[idees.columns.str.contains("|".join(efficiency_keywords))]
    to_drop = to_drop.append(pd.Index(["passenger cars", "passenger car efficiency"]))

    df = idees.reindex(new_index).drop(to_drop, axis=1)

    in_eurostat = df.index.levels[0].intersection(eurostat_countries)

    # add international navigation

    slicer = idx[in_eurostat, :, :, "Bunkers", :]
    fill_values = eurostat.loc[slicer, "Total all products"].groupby(level=[0, 1]).sum()
    # fill missing years for some countries by mean over the other years
    fill_values = fill_missing_years(fill_values)
    df.loc[in_eurostat, "total international navigation"] = fill_values

    # add swiss energy data

    df = pd.concat([df.drop("CH", errors="ignore"), swiss]).sort_index()

    # get values for missing countries based on Eurostat EnergyBalances

    # agriculture

    to_fill = df.index[
        df["total agriculture"].isna()
        & df.index.get_level_values("country").isin(eurostat_countries)
    ]
    c = to_fill.get_level_values("country")
    y = to_fill.get_level_values("year")

    # take total final energy consumption from Eurostat
    eurostat_sector = "Agriculture & forestry"
    slicer = idx[c, y, :, :, eurostat_sector]

    fill_values = eurostat.loc[slicer]["Total all products"].groupby(level=[0, 1]).sum()
    # fill missing years for some countries by mean over the other years
    fill_values = fill_missing_years(fill_values)
    df.loc[to_fill, "total agriculture"] = fill_values

    # split into end uses by average EU data from IDEES
    uses = ["electricity", "heat", "machinery"]

    for use in uses:
        avg = (
            idees["total agriculture electricity"] / idees["total agriculture"]
        ).mean()
        df.loc[to_fill, f"total agriculture {use}"] = (
            df.loc[to_fill, "total agriculture"] * avg
        )

    # divide cooking/space/water according to averages in EU28

    uses = ["space", "cooking", "water"]

    to_fill = df.index[
        df["total residential"].isna()
        & df.index.get_level_values("country").isin(eurostat_countries)
    ]
    c = to_fill.get_level_values("country")
    y = to_fill.get_level_values("year")

    for sector in ["residential", "services", "road", "rail"]:
        eurostat_sector = sector.capitalize()

        # fuel use

        for fuel in ["electricity", "total"]:
            slicer = idx[c, y, :, :, eurostat_sector]
            fill_values = (
                eurostat.loc[slicer, eurostat_fuels[fuel]].groupby(level=[0, 1]).sum()
            )
            # fill missing years for some countries by mean over the other years
            fill_values = fill_missing_years(fill_values)
            df.loc[to_fill, f"{fuel} {sector}"] = fill_values

    for sector in ["residential", "services"]:
        # electric use

        for use in uses:
            fuel_use = df[f"electricity {sector} {use}"]
            fuel = (
                df[f"electricity {sector}"].replace(0, np.nan).infer_objects(copy=False)
            )
            avg = fuel_use.div(fuel).mean()
            logger.debug(
                f"{sector}: average fraction of electricity for {use} is {avg:.3f}"
            )
            df.loc[to_fill, f"electricity {sector} {use}"] = (
                avg * df.loc[to_fill, f"electricity {sector}"]
            )

        # non-electric use

        for use in uses:
            nonelectric_use = (
                df[f"total {sector} {use}"] - df[f"electricity {sector} {use}"]
            )
            nonelectric = df[f"total {sector}"] - df[f"electricity {sector}"]
            nonelectric = nonelectric.copy().replace(0, np.nan)
            avg = nonelectric_use.div(nonelectric).mean()
            logger.debug(
                f"{sector}: average fraction of non-electric for {use} is {avg:.3f}"
            )
            electric_use = df.loc[to_fill, f"electricity {sector} {use}"]
            nonelectric = (
                df.loc[to_fill, f"total {sector}"]
                - df.loc[to_fill, f"electricity {sector}"]
            )
            df.loc[to_fill, f"total {sector} {use}"] = electric_use + avg * nonelectric

    # Fix Norway space and water heating fractions
    # http://www.ssb.no/en/energi-og-industri/statistikker/husenergi/hvert-3-aar/2014-07-14
    # The main heating source for about 73 per cent of the households is based on electricity
    # => 26% is non-electric

    if "NO" in df.index:
        elec_fraction = 0.73

        no_norway = df.drop("NO")

        for sector in ["residential", "services"]:
            # assume non-electric is heating
            nonelectric = (
                df.loc["NO", f"total {sector}"] - df.loc["NO", f"electricity {sector}"]
            )
            total_heating = nonelectric / (1 - elec_fraction)

            for use in uses:
                nonelectric_use = (
                    no_norway[f"total {sector} {use}"]
                    - no_norway[f"electricity {sector} {use}"]
                )
                nonelectric = (
                    no_norway[f"total {sector}"] - no_norway[f"electricity {sector}"]
                )
                nonelectric = nonelectric.copy().replace(0, np.nan)
                fraction = nonelectric_use.div(nonelectric).mean()
                df.loc["NO", f"total {sector} {use}"] = (
                    total_heating * fraction
                ).values
                df.loc["NO", f"electricity {sector} {use}"] = (
                    total_heating * fraction * elec_fraction
                ).values

    # Missing aviation

    slicer = idx[c, y, :, :, "Domestic aviation"]
    fill_values = eurostat.loc[slicer, "Total all products"].groupby(level=[0, 1]).sum()
    # fill missing years for some countries by mean over the other years
    fill_values = fill_missing_years(fill_values)
    df.loc[to_fill, "total domestic aviation"] = fill_values

    slicer = idx[c, y, :, :, "International aviation"]
    fill_values = eurostat.loc[slicer, "Total all products"].groupby(level=[0, 1]).sum()
    # fill missing years for some countries by mean over the other years
    fill_values = fill_missing_years(fill_values)
    df.loc[to_fill, "total international aviation"] = fill_values

    # missing domestic navigation

    slicer = idx[c, y, :, :, "Domestic Navigation"]
    fill_values = eurostat.loc[slicer, "Total all products"].groupby(level=[0, 1]).sum()
    # fill missing years for some countries by mean over the other years
    fill_values = fill_missing_years(fill_values)
    df.loc[to_fill, "total domestic navigation"] = fill_values

    # split road traffic for non-IDEES
    missing = df.index[df["total passenger cars"].isna()]
    for fuel in ["total", "electricity"]:
        selection = [
            f"{fuel} passenger cars",
            f"{fuel} other road passenger",
            f"{fuel} light duty road freight",
        ]
        if fuel == "total":
            selection.extend([f"{fuel} two-wheel", f"{fuel} heavy duty road freight"])
        road = df[selection].sum()
        road_fraction = road / road.sum()
        fill_values = cartesian(df.loc[missing, f"{fuel} road"], road_fraction)
        df.loc[missing, road_fraction.index] = fill_values

    # split rail traffic for non-IDEES
    missing = df.index[df["total rail passenger"].isna()]
    for fuel in ["total", "electricity"]:
        selection = [f"{fuel} rail passenger", f"{fuel} rail freight"]
        rail = df[selection].sum()
        rail_fraction = rail / rail.sum()
        fill_values = cartesian(df.loc[missing, f"{fuel} rail"], rail_fraction)
        df.loc[missing, rail_fraction.index] = fill_values

    # split aviation traffic for non-IDEES
    missing = df.index[df["total domestic aviation passenger"].isna()]
    for destination in ["domestic", "international"]:
        selection = [
            f"total {destination} aviation passenger",
            f"total {destination} aviation freight",
        ]
        aviation = df[selection].sum()
        aviation_fraction = aviation / aviation.sum()
        fill_values = cartesian(
            df.loc[missing, f"total {destination} aviation"], aviation_fraction
        )
        df.loc[missing, aviation_fraction.index] = fill_values

    for purpose in ["passenger", "freight"]:
        attrs = [
            f"total domestic aviation {purpose}",
            f"total international aviation {purpose}",
        ]
        df.loc[missing, f"total aviation {purpose}"] = df.loc[missing, attrs].sum(
            axis=1
        )

    if "BA" in df.index:
        # fill missing data for BA (services and road energy data)
        # proportional to RS with ratio of total residential demand
        mean_BA = df.loc["BA"].loc[2014:2021, "total residential"].mean()
        mean_RS = df.loc["RS"].loc[2014:2021, "total residential"].mean()
        ratio = mean_BA / mean_RS
        df.loc["BA"] = (
            df.loc["BA"].replace(0.0, np.nan).infer_objects(copy=False).values
        )
        df.loc["BA"] = df.loc["BA"].combine_first(ratio * df.loc["RS"]).values

    return df


def build_district_heat_share(countries: list[str], idees: pd.DataFrame) -> pd.Series:
    """
    Calculate the share of district heating for each country.

    Parameters
    ----------
    countries : list[str]
        List of country codes for which to calculate district heating share.
    idees : pd.DataFrame
        IDEES energy data dataframe.

    Returns
    -------
    pd.Series
        Series with the district heating share for each country.

    Notes
    -----
    - The function calculates the district heating share as the sum of residential and services distributed heat, divided by the sum of residential and services thermal uses.
    - The district heating share is then reindexed to match the provided list of countries.
    - Missing district heating shares are filled from `data/district_heat_share.csv`.
    - The function makes a conservative assumption and takes the minimum district heating share from both the IDEES data and `data/district_heat_share.csv`.
    """

    # district heating share
    district_heat = idees[
        ["distributed heat residential", "distributed heat services"]
    ].sum(axis=1)
    total_heat = (
        idees[["thermal uses residential", "thermal uses services"]]
        .sum(axis=1)
        .replace(0, np.nan)
    )

    district_heat_share = district_heat / total_heat

    district_heat_share = district_heat_share.reindex(countries, level="country")

    # Missing district heating share
    dh_share = (
        pd.read_csv(snakemake.input.district_heat_share, index_col=0, usecols=[0, 1])
        .div(100)
        .squeeze()
    )
    # make conservative assumption and take minimum from both data sets
    new_index = pd.MultiIndex.from_product(
        [dh_share.index, district_heat_share.index.get_level_values(1).unique()]
    )
    district_heat_share = pd.concat(
        [district_heat_share, dh_share.reindex(new_index, level=0)], axis=1
    ).min(axis=1)

    district_heat_share = district_heat_share.reindex(countries, level=0)

    district_heat_share.name = "district heat share"

    # restrict to available years
    district_heat_share = (
        district_heat_share.unstack().dropna(how="all", axis=1).ffill(axis=1)
    )

    return district_heat_share


def build_eea_co2(
    input_co2: str, year: int = 1990, emissions_scope: str = "CO2"
) -> pd.DataFrame:
    """
    Calculate CO2 emissions for a given year based on EEA data in Mt.

    Parameters
    ----------
    input_co2 : str
        Path to the input CSV file with CO2 data.
    year : int, optional
        Year for which to calculate emissions, by default 1990.
    emissions_scope : str, optional
        Scope of the emissions to consider, by default "CO2".

    Returns
    -------
    pd.DataFrame
        DataFrame with CO2 emissions for the given year.

    Notes
    -----
    - The function reads the `input_co2` data and for a specific `year` and `emission scope`
    - It calculates "industrial non-elec" and "agriculture" emissions from that data
    - It drops unneeded columns and converts the emissions to Mt.

    References
    ----------
    - `EEA CO2 data <https://www.eea.europa.eu/data-and-maps/data/national-emissions-reported-to-the-unfccc-and-to-the-eu-greenhouse-gas-monitoring-mechanism-16>`_ (downloaded 201228, modified by EEA last on 201221)
    """

    df = pd.read_csv(input_co2, encoding="latin-1", low_memory=False)

    df.replace(dict(Year="1985-1987"), 1986, inplace=True)
    df.Year = df.Year.astype(int)
    index_col = ["Country_code", "Pollutant_name", "Year", "Sector_name"]
    df = df.set_index(index_col).sort_index()

    cts = ["CH", "EUA", "NO"] + eu28_eea

    slicer = idx[cts, emissions_scope, year, to_ipcc.values()]
    emissions = (
        df.loc[slicer, "emissions"]
        .unstack("Sector_name")
        .rename(columns=reverse(to_ipcc))
        .droplevel([1, 2])
    )

    emissions.rename(index={"EUA": "EU28", "UK": "GB"}, inplace=True)

    to_subtract = [
        "electricity",
        "services non-elec",
        "residential non-elec",
        "road non-elec",
        "rail non-elec",
        "domestic aviation",
        "international aviation",
        "domestic navigation",
        "international navigation",
        "agriculture, forestry and fishing",
    ]
    emissions["industrial non-elec"] = emissions["total energy"] - emissions[
        to_subtract
    ].sum(axis=1)

    emissions["agriculture"] += emissions["agriculture, forestry and fishing"]

    to_drop = [
        "total energy",
        "total wL",
        "total woL",
        "agriculture, forestry and fishing",
    ]
    emissions.drop(columns=to_drop, inplace=True)

    # convert from Gt to Mt
    return emissions / 1e3


def build_eurostat_co2(eurostat: pd.DataFrame, year: int = 1990) -> pd.Series:
    """
    Calculate CO2 emissions for a given year based on Eurostat fuel consumption
    data and fuel-specific emissions.

    Parameters
    ----------
    eurostat : pd.DataFrame
        DataFrame with Eurostat data.
    year : int, optional
        Year for which to calculate emissions, by default 1990.

    Returns
    -------
    pd.Series
        Series with CO2 emissions for the given year.

    Notes
    -----
    - The function hard-sets fuel-specific emissions:
        - solid fuels: 0.36 tCO2_equi/MW_th (approximates coal)
        - oil: 0.285 tCO2_equi/MW_th (average of distillate and residue)
        - natural gas: 0.2 tCO2_equi/MW_th
    - It then multiplies the Eurostat fuel consumption data for `year` by the specific emissions and sums the result.

    References
    ----------
    - Oil values from `EIA <https://www.eia.gov/tools/faqs/faq.cfm?id=74&t=11>`_
    - Distillate oil (No. 2)  0.276
    - Residual oil (No. 6)  0.298
    - `EIA Electricity Annual <https://www.eia.gov/electricity/annual/html/epa_a_03.html>`_
    """

    eurostat_year = eurostat.xs(year, level="year")

    specific_emissions = pd.Series(index=eurostat.columns, dtype=float)

    # emissions in tCO2_equiv per MWh_th
    specific_emissions["Solid fossil fuels"] = 0.36  # Approximates coal
    specific_emissions["Oil and petroleum products"] = (
        0.285  # Average of distillate and residue
    )
    specific_emissions["Natural gas"] = 0.2  # For natural gas

    return eurostat_year.multiply(specific_emissions).sum(axis=1)


def build_co2_totals(
    countries: list[str], eea_co2: pd.DataFrame, eurostat_co2: pd.DataFrame
) -> pd.DataFrame:
    """
    Combine CO2 emissions data from EEA and Eurostat for a list of countries.

    Parameters
    ----------
    countries : list[str]
        List of country codes for which CO2 totals need to be built.
    eea_co2 : pd.DataFrame
        DataFrame with EEA CO2 emissions data.
    eurostat_co2 : pd.DataFrame
        DataFrame with Eurostat CO2 emissions data.

    Returns
    -------
    pd.DataFrame
        Combined CO2 emissions data for the given countries.

    Notes
    -----
    - The function combines the CO2 emissions from EEA and Eurostat into a single DataFrame for the given countries.
    """

    co2 = eea_co2.reindex(countries)

    for ct in pd.Index(countries).intersection(
        ["BA", "RS", "XK", "AL", "ME", "MK", "UA", "MD"]
    ):
        mappings = {
            "electricity": (ct, "+", "Electricity & heat generation", np.nan),
            "residential non-elec": (ct, "+", "+", "Residential"),
            "services non-elec": (ct, "+", "+", "Services"),
            "road non-elec": (ct, "+", "+", "Road"),
            "rail non-elec": (ct, "+", "+", "Rail"),
            "domestic navigation": (ct, "+", "+", "Domestic Navigation"),
            "international navigation": (ct, "-", "Bunkers"),
            "domestic aviation": (ct, "+", "+", "Domestic aviation"),
            "international aviation": (ct, "-", "International aviation"),
            # does not include industrial process emissions or fuel processing/refining
            "industrial non-elec": (ct, "+", "Industry sector"),
            # does not include non-energy emissions
            "agriculture": (eurostat_co2.index.get_level_values(0) == ct)
            & eurostat_co2.index.isin(["Agriculture & forestry", "Fishing"], level=3),
        }

        for i, mi in mappings.items():
            co2.at[ct, i] = eurostat_co2.loc[mi].sum()

    return co2


def build_transport_data(
    countries: list[str], population: pd.DataFrame, idees: pd.DataFrame
) -> pd.DataFrame:
    """
    Build transport data for a set of countries based on IDEES data.

    Parameters
    ----------
    countries : list[str]
        List of country codes.
    population : pd.DataFrame
        DataFrame with population data.
    idees : pd.DataFrame
        DataFrame with IDEES data.

    Returns
    -------
    pd.DataFrame
        DataFrame with transport data.

    Notes
    -----
    - The function first collects the number of passenger cars.
    - For Switzerland, it reads the data from `data/gr-e-11.03.02.01.01-cc.csv`.
    - It fills missing data on the number of cars and fuel efficiency with average data.

    References
    ----------
    - Swiss transport data: `BFS <https://www.bfs.admin.ch/bfs/en/home/statistics/mobility-transport/transport-infrastructure-vehicles/vehicles/road-vehicles-stock-level-motorisation.html>`_
    """
    years = np.arange(2000, 2022)

    # first collect number of cars
    transport_data = pd.DataFrame(idees["passenger cars"])

    countries_without_ch = set(countries) - {"CH"}
    new_index = pd.MultiIndex.from_product(
        [countries_without_ch, transport_data.index.unique(1)],
        names=["country", "year"],
    )

    transport_data = transport_data.reindex(index=new_index)

    if "CH" in countries:
        fn = snakemake.input.swiss_transport
        swiss_cars = pd.read_csv(fn, index_col=0).loc[years, ["passenger cars"]]

        swiss_cars.index = pd.MultiIndex.from_product(
            [["CH"], swiss_cars.index], names=["country", "year"]
        )

        transport_data = pd.concat([transport_data, swiss_cars]).sort_index()

    transport_data = transport_data.rename(columns={"passenger cars": "number cars"})
    # clean up dataframe
    transport_data = transport_data[
        transport_data.index.get_level_values(1).isin(years)
    ]

    missing = transport_data.index[transport_data["number cars"].isna()]
    if not missing.empty:
        logger.info(
            f"Missing data on cars from:\n{list(missing)}\nFilling gaps with averaged data."
        )
        cars_pp = transport_data["number cars"] / population

        fill_values = {
            year: cars_pp.mean() * population for year in transport_data.index.unique(1)
        }
        fill_values = pd.DataFrame(fill_values).stack()
        fill_values = pd.DataFrame(fill_values, columns=["number cars"])
        fill_values.index.names = ["country", "year"]
        fill_values = fill_values.reindex(transport_data.index)

        transport_data = transport_data.combine_first(fill_values)

    # collect average fuel efficiency in MWh/100km, taking passengar car efficiency in TWh/100km
    transport_data["average fuel efficiency"] = idees["passenger car efficiency"] * 1e6

    missing = transport_data.index[transport_data["average fuel efficiency"].isna()]
    if not missing.empty:
        logger.info(
            f"Missing data on fuel efficiency from:\n{list(missing)}\nFilling gaps with averaged data."
        )

        fill_values = transport_data["average fuel efficiency"].mean()
        transport_data.loc[missing, "average fuel efficiency"] = fill_values

    return transport_data


def rescale_idees_from_eurostat(
    idees_countries: list[str], energy: pd.DataFrame, eurostat: pd.DataFrame
) -> pd.DataFrame:
    """
    Takes JRC IDEES data from 2021 and rescales it by the ratio of the Eurostat
    data and the 2021 Eurostat data.
    Missing data: ['passenger car efficiency', 'passenger cars']

    Parameters
    ----------
    idees_countries : list[str]
        List of IDEES country codes.
    energy : pd.DataFrame
        DataFrame with JRC IDEES data.
    eurostat : pd.DataFrame
        DataFrame with Eurostat data.

    Returns
    -------
    pd.DataFrame
        DataFrame with rescaled IDEES data.

    Notes
    -----
    - The function first reads in the Eurostat data for 2015 and calculates the ratio of that data with other Eurostat data.
    - This ratio is mapped to the IDEES data.

    References
    ----------
    - JRC IDEES data: `JRC IDEES <https://ec.europa.eu/jrc/en/publication/eur-scientific-and-technical-research-reports/jrc-idees>`_
    - Eurostat data: `Eurostat <https://ec.europa.eu/eurostat/data/database>`_
    """

    main_cols = ["Total all products", "Electricity"]
    # read in the eurostat data for 2015
    eurostat_2021 = eurostat.xs(2021, level="year")[main_cols]
    # calculate the ratio of the two data sets
    ratio = eurostat[main_cols] / eurostat_2021
    ratio = ratio.droplevel([2, 5])
    cols_rename = {"Total all products": "total", "Electricity": "ele"}
    index_rename = {v: k for k, v in idees_rename.items()}
    ratio.rename(columns=cols_rename, index=index_rename, inplace=True)

    mappings = {
        "Residential": {
            "total": [
                "total residential space",
                "total residential water",
                "total residential cooking",
                "total residential",
                "distributed heat residential",
                "thermal uses residential",
            ],
            "elec": [
                "electricity residential space",
                "electricity residential water",
                "electricity residential cooking",
                "electricity residential",
            ],
        },
        "Services": {
            "total": [
                "total services space",
                "total services water",
                "total services cooking",
                "total services",
                "distributed heat services",
                "thermal uses services",
            ],
            "elec": [
                "electricity services space",
                "electricity services water",
                "electricity services cooking",
                "electricity services",
            ],
        },
        "Agriculture & forestry": {
            "total": [
                "total agriculture heat",
                "total agriculture machinery",
                "total agriculture",
            ],
            "elec": [
                "total agriculture electricity",
            ],
        },
        "Road": {
            "total": [
                "total road",
                "total passenger cars",
                "total other road passenger",
                "total light duty road freight",
                "total heavy duty road freight",
            ],
            "elec": [
                "electricity road",
                "electricity passenger cars",
                "electricity other road passenger",
                "electricity light duty road freight",
            ],
        },
        "Rail": {
            "total": [
                "total rail",
                "total rail passenger",
                "total rail freight",
            ],
            "elec": [
                "electricity rail",
                "electricity rail passenger",
                "electricity rail freight",
            ],
        },
    }

    avia_inter = [
        "total aviation passenger",
        "total aviation freight",
        "total international aviation passenger",
        "total international aviation freight",
        "total international aviation",
    ]
    avia_domestic = [
        "total domestic aviation passenger",
        "total domestic aviation freight",
        "total domestic aviation",
    ]
    navigation = [
        "total domestic navigation",
    ]
    # international navigation is already read in from the eurostat data directly

    for country in idees_countries:
        filling_years = [(2015, slice(2016, 2021)), (2000, slice(1990, 1999))]

        for source_year, target_years in filling_years:
            slicer_source = idx[country, source_year, :, :]
            slicer_target = idx[country, target_years, :, :]

            for sector, mapping in mappings.items():
                sector_ratio = ratio.loc[
                    (country, slice(None), slice(None), sector)
                ].droplevel("lvl2")

                energy.loc[slicer_target, mapping["total"]] = cartesian(
                    sector_ratio.loc[target_years, "total"],
                    energy.loc[slicer_source, mapping["total"]].squeeze(axis=0),
                ).values
                energy.loc[slicer_target, mapping["elec"]] = cartesian(
                    sector_ratio.loc[target_years, "ele"],
                    energy.loc[slicer_source, mapping["elec"]].squeeze(axis=0),
                ).values

            level_drops = ["country", "lvl2", "lvl3"]

            slicer = idx[country, :, :, "Domestic aviation"]
            avi_d = ratio.loc[slicer, "total"].droplevel(level_drops)

            slicer = idx[country, :, :, "International aviation"]
            avi_i = ratio.loc[slicer, "total"].droplevel(level_drops)

            slicer = idx[country, :, :, "Domestic Navigation"]
            nav = ratio.loc[slicer, "total"].droplevel(level_drops)

            energy.loc[slicer_target, avia_inter] = cartesian(
                avi_i.loc[target_years],
                energy.loc[slicer_source, avia_inter].squeeze(axis=0),
            ).values

            energy.loc[slicer_target, avia_domestic] = cartesian(
                avi_d.loc[target_years],
                energy.loc[slicer_source, avia_domestic].squeeze(axis=0),
            ).values

            energy.loc[slicer_target, navigation] = cartesian(
                nav.loc[target_years],
                energy.loc[slicer_source, navigation].squeeze(axis=0),
            ).values

        # set the total of agriculture/road to the sum of all agriculture/road categories (corresponding to the IDEES data)
        rows = idx[country, :]
        cols = [
            "total agriculture electricity",
            "total agriculture heat",
            "total agriculture machinery",
        ]
        energy.loc[rows, "total agriculture"] = energy.loc[rows, cols].sum(axis=1)

        cols = [
            "total passenger cars",
            "total other road passenger",
            "total light duty road freight",
            "total heavy duty road freight",
        ]
        energy.loc[rows, "total road"] = energy.loc[rows, cols].sum(axis=1)

    return energy


def update_residential_from_eurostat(energy: pd.DataFrame) -> pd.DataFrame:
    """
    Updates energy balances for residential from disaggregated data from
    Eurostat by mutating input data DataFrame.

    Parameters
    ----------
    energy : pd.DataFrame
        DataFrame with energy data.

    Returns
    -------
    pd.DataFrame
        DataFrame with updated energy balances.

    Notes
    -----
    - The function first reads in the Eurostat data for households and maps the energy types to the corresponding Eurostat codes.
    - For each energy type, it selects the corresponding data, converts units, and drops unnecessary data.
    """
    eurostat_households = pd.read_csv(snakemake.input.eurostat_households)

    # Column mapping for energy type
    nrg_type = {
        "total residential": ("FC_OTH_HH_E", "TOTAL"),
        "total residential space": ("FC_OTH_HH_E_SH", "TOTAL"),
        "total residential water": ("FC_OTH_HH_E_WH", "TOTAL"),
        "total residential cooking": ("FC_OTH_HH_E_CK", "TOTAL"),
        "electricity residential": ("FC_OTH_HH_E", "E7000"),
        "electricity residential space": ("FC_OTH_HH_E_SH", "E7000"),
        "electricity residential water": ("FC_OTH_HH_E_WH", "E7000"),
        "electricity residential cooking": ("FC_OTH_HH_E_CK", "E7000"),
    }

    for nrg_name, (code, siec) in nrg_type.items():
        # Select energy balance type, rename columns and countries to match IDEES data,
        # convert TJ to TWh
        col_to_rename = {"geo": "country", "TIME_PERIOD": "year", "OBS_VALUE": nrg_name}
        idx_to_rename = {v: k for k, v in idees_rename.items()}
        drop_geo = ["EU27_2020", "EA20"]  # noqa: F841
        nrg_data = eurostat_households.query(
            "nrg_bal == @code and siec == @siec and geo not in @drop_geo and OBS_VALUE > 0"
        ).copy()
        nrg_data.rename(columns=col_to_rename, inplace=True)
        nrg_data = nrg_data.set_index(["country", "year"])[nrg_name] / 3.6e3
        nrg_data.rename(index=idx_to_rename, inplace=True)

        # update energy balance from household-specific eurostat data
        idx = nrg_data.index.intersection(energy.index)
        energy.loc[idx, nrg_name] = nrg_data[idx]

    logger.info(
        "Updated energy balances for residential using disaggregate final energy consumption data in Households from Eurostat"
    )


def build_transformation_output_coke(eurostat, fn):
    """
    Extracts and builds the transformation output data for coke ovens from the
    Eurostat dataset.

    This function specifically filters the Eurostat data to extract
    transformation output related to coke ovens.
    Since the transformation output for coke ovens
    is not included in the final energy consumption of the iron and steel sector,
    it needs to be processed and added separately. The filtered data is saved
    as a CSV file.

    Parameters
    ----------
    eurostat (pd.DataFrame): A pandas DataFrame containing Eurostat data with
                             a multi-level index
    fn (str): The file path where the resulting CSV file should be saved.

    Output:
    The resulting transformation output data for coke ovens is saved as a CSV
    file at the path specified in fn.
    """
    slicer = pd.IndexSlice[:, :, :, "Coke ovens", "Other sources", :]
    df = eurostat.loc[slicer, :].droplevel(level=[2, 3, 4, 5])
    df.to_csv(fn)


def build_heating_efficiencies(
    countries: list[str], idees: pd.DataFrame
) -> pd.DataFrame:
    """
    Build heating efficiencies for a set of countries based on IDEES data.

    Parameters
    ----------
    countries : list[str]
        List of country codes.
    idees : pd.DataFrame
        DataFrame with IDEES data.

    Returns
    -------
    pd.DataFrame
        DataFrame with heating efficiencies.


    Notes
    -----
    - It fills missing data with average data.
    """

    cols = idees.columns[
        idees.columns.str.contains("space efficiency")
        ^ idees.columns.str.contains("water efficiency")
    ]

    heating_efficiencies = pd.DataFrame(idees[cols])

    new_index = pd.MultiIndex.from_product(
        [countries, heating_efficiencies.index.unique(1)],
        names=["country", "year"],
    )

    heating_efficiencies = heating_efficiencies.reindex(index=new_index)

    for col in cols:
        unstacked = heating_efficiencies[col].unstack()

        fillvalue = unstacked.mean()

        for ct in unstacked.index:
            mask = unstacked.loc[ct].isna()
            unstacked.loc[ct, mask] = fillvalue[mask]
        heating_efficiencies[col] = unstacked.stack()

    return heating_efficiencies


# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_energy_totals")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params.energy

    nuts3 = gpd.read_file(snakemake.input.nuts3_shapes).set_index("index")
    population = nuts3["pop"].groupby(nuts3.country).sum()

    countries = snakemake.params.countries
    idees_countries = pd.Index(countries).intersection(eu27)

    input_eurostat = snakemake.input.eurostat
    eurostat = build_eurostat(
        input_eurostat,
        countries,
        nprocesses=snakemake.threads,
        disable_progressbar=snakemake.config["run"].get("disable_progressbar", False),
    )

    build_transformation_output_coke(
        eurostat, snakemake.output.transformation_output_coke
    )

    swiss = build_swiss()
    idees = build_idees(idees_countries)

    energy = build_energy_totals(countries, eurostat, swiss, idees)

    update_residential_from_eurostat(energy)

    energy.to_csv(snakemake.output.energy_name)

    # use rescaled idees data to calculate district heat share
    district_heat_share = build_district_heat_share(
        countries, energy.loc[idees_countries]
    )
    district_heat_share.to_csv(snakemake.output.district_heat_share)

    base_year_emissions = params["base_emissions_year"]
    emissions_scope = snakemake.params.energy["emissions"]
    eea_co2 = build_eea_co2(snakemake.input.co2, base_year_emissions, emissions_scope)
    eurostat_co2 = build_eurostat_co2(eurostat, base_year_emissions)

    co2 = build_co2_totals(countries, eea_co2, eurostat_co2)
    co2.to_csv(snakemake.output.co2_name)

    transport = build_transport_data(countries, population, idees)
    transport.to_csv(snakemake.output.transport_name)

    heating_efficiencies = build_heating_efficiencies(countries, idees)
    heating_efficiencies.to_csv(snakemake.output.heating_efficiencies)
