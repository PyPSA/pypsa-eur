#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script calculates the space heating savings through better insulation of
the thermal envelope of a building and corresponding costs for different
building types in different countries.

Methodology
-----------

The energy savings calculations are based on the

  EN ISO 13790 / seasonal method https://www.iso.org/obp/ui/#iso:std:iso:13790:ed-2:v1:en:

  - calculations heavily oriented on the TABULAWebTool
  http://webtool.building-typology.eu/
  http://www.episcope.eu/fileadmin/tabula/public/docs/report/TABULA_CommonCalculationMethod.pdf
  which is following the EN ISO 13790 / seasonal method

  - building stock data:
      mainly: hotmaps project https://gitlab.com/hotmaps/building-stock
      missing: EU building observatory https://ec.europa.eu/energy/en/eu-buildings-database

  - building types with typical surfaces/ standard values:
      - tabula https://episcope.eu/fileadmin/tabula/public/calc/tabula-calculator.xlsx


Basic Equations
---------------

The basic equations:

    The Energy needed for space heating E_space [W/m²] are calculated as the
    sum of heat losses and heat gains:

        E_space = H_losses - H_gains

    Heat losses constitute from the losses through heat transmission (H_tr [W/m²K])
    (this includes heat transfer through building elements and thermal bridges)
    and losses by ventilation (H_ve [W/m²K]):

        H_losses = (H_tr + H_ve) * F_red * (T_threshold - T_averaged_d_heat) * d_heat * 1/365

        F_red : reduction factor, considering non-uniform heating [°C], p.16 chapter 2.6 [-]
        T_threshold : heating temperature threshold, assumed 15 C
        d_heat : Length of heating season, number of days with daily averaged temperature below T_threshold
        T_averaged_d_heat : mean daily averaged temperature of the days within heating season d_heat

    Heat gains constitute from the gains by solar radiation (H_solar) and
    internal heat gains (H_int) weighted by a gain utilisation factor nu:

        H_gains = nu * (H_solar + H_int)

Structure
---------

The script has the following structure:

    (0) fixed parameters are set
    (1) prepare data, bring to same format
    (2) calculate space heat demand depending on additional insulation material
    (3) calculate costs for corresponding additional insulation material
    (4) get cost savings per retrofitting measures for each sector by weighting
        with heated floor area
"""
import pandas as pd
import xarray as xr

# (i) --- FIXED PARAMETER / STANDARD VALUES -----------------------------------

# thermal conductivity standard value
k = 0.035
# strength of relative retrofitting depending on the component
# determined by historical data of insulation thickness for retrofitting
l_weight = pd.DataFrame({"weight": [1.95, 1.48, 1.0]}, index=["Roof", "Wall", "Floor"])

# standard room height [m], used to calculate heat transfer by ventilation
h_room = 2.5
# volume specific heat capacity air [Wh/m^3K]
c_p_air = 0.34
# internal heat capacity per m² A_c_ref [Wh/(m^2K)]
c_m = 45
# average thermal output of the internal heat sources per m^2 reference area [W/m^2]
phi_int = 3
# constant parameter tau_H_0 [h] according to EN 13790 seasonal method
tau_H_0 = 30
# constant parameter alpha_H_0 [-] according to EN 13790 seasonal method
alpha_H_0 = 0.8

# parameter for solar heat load during heating season -------------------------
# tabular standard values table p.8 in documentation
external_shading = 0.6  # vertical orientation: fraction of window area shaded [-]
frame_area_fraction = 0.3  # fraction of frame area of window [-]
non_perpendicular = (
    0.9  # reduction factor, considering radiation non perpendicular to the glazing[-]
)
solar_energy_transmittance = (
    0.5  # solar energy transmiitance for radiation perpecidular to the glazing [-]
)
# solar global radiation [kWh/(m^2a)]
solar_global_radiation = pd.Series(
    [271, 392, 271, 160],
    index=["east", "south", "west", "north"],
    name="solar_global_radiation [kWh/(m^2a)]",
)

# threshold temperature for heating [Celsius] --------------------------------
t_threshold = 15

# rename sectors
# rename residential sub sectors
rename_sectors = {
    "Single family- Terraced houses": "SFH",
    "Multifamily houses": "MFH",
    "Appartment blocks": "AB",
}


# additional insulation thickness, determines maximum possible savings [m]
l_strength = ["0.07", "0.075", "0.08", "0.1", "0.15", "0.22", "0.24", "0.26"]


# (ii) --- FUNCTIONS ----------------------------------------------------------


def get_average_temperature_during_heating_season(temperature, t_threshold=15):
    """
    returns average temperature during heating season
    input:
        temperature : pd.Series(Index=time, values=temperature)
        t_threshold : threshold temperature for heating degree days (HDD)
    returns:
        average temperature
    """
    t_average_daily = temperature.resample("1D").mean()
    return t_average_daily.loc[t_average_daily < t_threshold].mean()


def prepare_building_stock_data():
    """
    reads building stock data and cleans up the format, returns
    --------
    u_values:          pd.DataFrame current U-values
    area_tot:          heated floor area per country and sector [Mm²]
    area:              heated floor area [Mm²] for country, sector, building
                       type and period

    """
    building_data = pd.read_csv(snakemake.input.building_stock, usecols=list(range(13)))

    # standardize data
    building_data["type"].replace(
        {
            "Covered area: heated  [Mm²]": "Heated area [Mm²]",
            "Windows ": "Window",
            "Windows": "Window",
            "Walls ": "Wall",
            "Walls": "Wall",
            "Roof ": "Roof",
            "Floor ": "Floor",
        },
        inplace=True,
    )
    building_data["feature"].replace(
        {
            "Construction features (U-value)": "Construction features (U-values)",
        },
        inplace=True,
    )

    building_data.country_code = building_data.country_code.str.upper()
    building_data["subsector"].replace(
        {"Hotels and Restaurants": "Hotels and restaurants"}, inplace=True
    )
    building_data["sector"].replace(
        {"Residential sector": "residential", "Service sector": "services"},
        inplace=True,
    )

    # extract u-values
    u_values = building_data[
        (building_data.feature.str.contains("U-values"))
        & (building_data.subsector != "Total")
    ]

    components = list(u_values.type.unique())

    country_iso_dic = building_data.set_index("country")["country_code"].to_dict()

    # add missing /rename countries
    country_iso_dic.update(
        {
            "Norway": "NO",
            "Iceland": "IS",
            "Montenegro": "ME",
            "Serbia": "RS",
            "Albania": "AL",
            "United Kingdom": "GB",
            "Bosnia and Herzegovina": "BA",
            "Switzerland": "CH",
        }
    )

    building_data["country_code"] = building_data["country"].map(country_iso_dic)

    # heated floor area ----------------------------------------------------------
    area = building_data[
        (building_data.type == "Heated area [Mm²]")
        & (building_data.subsector != "Total")
    ]
    area_tot = area[["country", "sector", "value"]].groupby(["country", "sector"]).sum()
    area = pd.concat(
        [
            area,
            area.apply(
                lambda x: x.value / area_tot.value.loc[(x.country, x.sector)], axis=1
            ).rename("weight"),
        ],
        axis=1,
    )
    area = area.groupby(["country", "sector", "subsector", "bage"]).sum()
    area_tot.rename(index=country_iso_dic, inplace=True)

    # add for some missing countries floor area from other data sources
    area_missing = pd.read_csv(
        snakemake.input.floor_area_missing,
        index_col=[0, 1],
        usecols=[0, 1, 2, 3],
        encoding="ISO-8859-1",
    )
    area_tot = pd.concat([area_tot, area_missing.unstack(level=-1).dropna().stack()])
    area_tot = area_tot.loc[~area_tot.index.duplicated(keep="last")]

    # for still missing countries calculate floor area by population size
    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    pop_layout["ct"] = pop_layout.index.str[:2]
    ct_total = pop_layout.total.groupby(pop_layout["ct"]).sum()

    area_per_pop = (
        area_tot.unstack()
        .reindex(index=ct_total.index)
        .apply(lambda x: x / ct_total[x.index])
    )
    missing_area_ct = ct_total.index.difference(area_tot.index.levels[0])
    for ct in missing_area_ct.intersection(ct_total.index):
        averaged_data = pd.DataFrame(
            area_per_pop.value.reindex(map_for_missings[ct]).mean() * ct_total[ct],
            columns=["value"],
        )
        index = pd.MultiIndex.from_product([[ct], averaged_data.index.to_list()])
        averaged_data.index = index
        averaged_data["estimated"] = 1
        if ct not in area_tot.index.levels[0]:
            area_tot = pd.concat([area_tot, averaged_data], sort=True)
        else:
            area_tot.loc[averaged_data.index] = averaged_data

    # u_values for Poland are missing -> take them from eurostat -----------
    u_values_PL = pd.read_csv(snakemake.input.u_values_PL)
    u_values_PL.component.replace({"Walls": "Wall", "Windows": "Window"}, inplace=True)
    area_PL = area.loc["Poland"].reset_index()
    data_PL = pd.DataFrame(columns=u_values.columns, index=area_PL.index)
    data_PL["country"] = "Poland"
    data_PL["country_code"] = "PL"
    # data from area
    for col in ["sector", "subsector", "bage"]:
        data_PL[col] = area_PL[col]
    data_PL["btype"] = area_PL["subsector"]

    data_PL_final = pd.DataFrame()
    for component in components:
        data_PL["type"] = component
        data_PL["value"] = data_PL.apply(
            lambda x: u_values_PL[
                (u_values_PL.component == component)
                & (u_values_PL.sector == x["sector"])
            ][x["bage"]].iloc[0],
            axis=1,
        )
        data_PL_final = pd.concat([data_PL_final, data_PL])

    u_values = pd.concat([u_values, data_PL_final]).reset_index(drop=True)

    # clean data ---------------------------------------------------------------
    # smallest possible today u values for windows 0.8 (passive house standard)
    # maybe the u values for the glass and not the whole window including frame
    # for those types assumed in the dataset
    u_values.loc[(u_values.type == "Window") & (u_values.value < 0.8), "value"] = 0.8
    # drop unnecessary columns
    u_values.drop(
        ["topic", "feature", "detail", "estimated", "unit"],
        axis=1,
        inplace=True,
        errors="ignore",
    )

    u_values.subsector.replace(rename_sectors, inplace=True)
    u_values.btype.replace(rename_sectors, inplace=True)

    # for missing weighting of surfaces of building types assume MFH
    u_values["assumed_subsector"] = u_values.subsector
    u_values.loc[
        ~u_values.subsector.isin(rename_sectors.values()), "assumed_subsector"
    ] = "MFH"

    u_values.country_code.replace({"UK": "GB"}, inplace=True)
    u_values.bage.replace({"Berfore 1945": "Before 1945"}, inplace=True)
    u_values = u_values[~u_values.bage.isna()]

    u_values.set_index(["country_code", "subsector", "bage", "type"], inplace=True)

    #  only take in config.yaml specified countries into account
    countries = snakemake.params.countries
    area_tot = area_tot.loc[countries]

    return u_values, country_iso_dic, countries, area_tot, area


def prepare_building_topology(u_values, same_building_topology=True):
    """
    Reads in typical building topologies (e.g. average surface of building
    elements) and typical losses through thermal bridging and air ventilation.
    """
    data_tabula = pd.read_csv(
        snakemake.input.data_tabula,
        skiprows=lambda x: x in range(1, 11),
        low_memory=False,
    ).iloc[:2974]

    parameters = [
        "Code_Country",
        # building type (SFH/MFH/AB)
        "Code_BuildingSizeClass",
        # time period of build year
        "Year1_Building",
        "Year2_Building",
        # areas [m^2]
        "A_C_Ref",  # conditioned area, internal
        "A_Roof_1",
        "A_Roof_2",
        "A_Wall_1",
        "A_Wall_2",
        "A_Floor_1",
        "A_Floor_2",
        "A_Window_1",
        "A_Window_2",
        # for air ventilation loses [1/h]
        "n_air_use",
        "n_air_infiltration",
        # for losses due to thermal bridges, standard values [W/(m^2K)]
        "delta_U_ThermalBridging",
        # floor area related heat transfer coefficient by transmission [-]
        "F_red_temp",
        # refurbishment state [1: not refurbished, 2: moderate ,3: strong refurbishment]
        "Number_BuildingVariant",
    ]

    data_tabula = data_tabula[parameters]

    building_elements = ["Roof", "Wall", "Floor", "Window"]

    # get total area of building components
    for element in building_elements:
        elements = ["A_{}_1".format(element), "A_{}_2".format(element)]
        data_tabula = pd.concat(
            [
                data_tabula.drop(elements, axis=1),
                data_tabula[elements].sum(axis=1).rename("A_{}".format(element)),
            ],
            axis=1,
        )

    # clean data
    data_tabula = data_tabula.loc[
        pd.concat(
            [
                data_tabula[col] != 0
                for col in ["A_Wall", "A_Floor", "A_Window", "A_Roof", "A_C_Ref"]
            ],
            axis=1,
        ).all(axis=1)
    ]
    data_tabula = data_tabula[data_tabula.Number_BuildingVariant.isin([1, 2, 3])]
    data_tabula = data_tabula[
        data_tabula.Code_BuildingSizeClass.isin(["AB", "SFH", "MFH", "TH"])
    ]

    # map tabula building periods to hotmaps building periods
    def map_periods(build_year1, build_year2):
        periods = {
            (0, 1945): "Before 1945",
            (1945, 1969): "1945 - 1969",
            (1970, 1979): "1970 - 1979",
            (1980, 1989): "1980 - 1989",
            (1990, 1999): "1990 - 1999",
            (2000, 2010): "2000 - 2010",
            (2010, 10000): "Post 2010",
        }
        minimum = 1e5
        for key in periods:
            diff = abs(build_year1 - key[0]) + abs(build_year2 - key[1])
            if diff < minimum:
                minimum = diff
                searched_period = periods[key]
        return searched_period

    data_tabula["bage"] = data_tabula.apply(
        lambda x: map_periods(x.Year1_Building, x.Year2_Building), axis=1
    )

    # set new index
    data_tabula = data_tabula.set_index(
        ["Code_Country", "Code_BuildingSizeClass", "bage", "Number_BuildingVariant"]
    )

    # get typical building topology
    area_cols = ["A_C_Ref", "A_Floor", "A_Roof", "A_Wall", "A_Window"]
    typical_building = (
        data_tabula.groupby(level=[1, 2])
        .mean()
        .rename(index={"TH": "SFH"})
        .groupby(level=[0, 1])
        .mean()
    )

    # drop duplicates
    data_tabula = data_tabula[~data_tabula.index.duplicated(keep="first")]

    # fill missing values
    hotmaps_data_i = (
        u_values.reset_index()
        .set_index(["country_code", "assumed_subsector", "bage"])
        .index
    )
    # missing countries in tabular
    missing_ct = data_tabula.unstack().reindex(hotmaps_data_i.unique())
    # areas should stay constant for different retrofitting measures
    cols_constant = [
        "Year1_Building",
        "Year2_Building",
        "A_C_Ref",
        "A_Roof",
        "A_Wall",
        "A_Floor",
        "A_Window",
    ]
    for col in cols_constant:
        missing_ct[col] = missing_ct[col].combine_first(
            missing_ct[col].groupby(level=[0, 1, 2]).mean()
        )
    missing_ct = (
        missing_ct.unstack().unstack().fillna(missing_ct.unstack().unstack().mean())
    )
    data_tabula = missing_ct.stack(level=[-1, -2, -3], dropna=False)

    # sets for different countries same building topology which only depends on
    # build year and subsector (MFH, SFH, AB)
    if same_building_topology:
        typical_building = (
            typical_building.reindex(data_tabula.droplevel(0).index)
        ).set_index(data_tabula.index)
        data_tabula.update(typical_building[area_cols])

    # total buildings envelope surface [m^2]
    data_tabula["A_envelope"] = data_tabula[
        ["A_{}".format(element) for element in building_elements]
    ].sum(axis=1)

    return data_tabula


def prepare_cost_retro(country_iso_dic):
    """
    Read and prepare retro costs, annualises them if annualise_cost=True.
    """
    cost_retro = pd.read_csv(
        snakemake.input.cost_germany, nrows=4, index_col=0, usecols=[0, 1, 2, 3]
    )
    cost_retro.rename(lambda x: x.capitalize(), inplace=True)

    window_assumptions = pd.read_csv(
        snakemake.input.window_assumptions, skiprows=[1], usecols=[0, 1, 2, 3], nrows=2
    )

    if annualise_cost:
        cost_retro[["cost_fix", "cost_var"]] = cost_retro[
            ["cost_fix", "cost_var"]
        ].apply(
            lambda x: x
            * interest_rate
            / (1 - (1 + interest_rate) ** -cost_retro.loc[x.index, "life_time"])
        )

    # weightings of costs ---------------------------------------------
    if construction_index:
        cost_w = pd.read_csv(
            snakemake.input.construction_index, skiprows=3, nrows=32, index_col=0
        )
        # since German retrofitting costs are assumed
        cost_w = (cost_w["2018"] / cost_w.loc["Germany", "2018"]).rename(
            index=country_iso_dic
        )
    else:
        cost_w = None

    if tax_weighting:
        tax_w = pd.read_csv(
            snakemake.input.tax_w, header=12, nrows=39, index_col=0, usecols=[0, 4]
        )
        tax_w.rename(index=country_iso_dic, inplace=True)
        tax_w = tax_w.apply(pd.to_numeric, errors="coerce").iloc[:, 0]
        tax_w.dropna(inplace=True)
    else:
        tax_w = None

    return cost_retro, window_assumptions, cost_w, tax_w


def prepare_temperature_data():
    """
    Returns the temperature dependent data for each country:

    d_heat             : length of heating season pd.Series(index=countries) [days/year]
                         on those days, daily average temperature is below
                         threshold temperature t_threshold
    temperature_factor : accumulated difference between internal and
                         external temperature pd.Series(index=countries) ([K]) * [days/year]

    temperature_factor = (t_threshold - temperature_average_d_heat) * d_heat * 1/365
    """
    temperature = xr.open_dataarray(snakemake.input.air_temperature).to_pandas()
    d_heat = (
        temperature.T.groupby(temperature.columns.str[:2])
        .mean()
        .T.resample("1D")
        .mean()
        < t_threshold
    ).sum()
    temperature_average_d_heat = (
        temperature.T.groupby(temperature.columns.str[:2])
        .mean()
        .T.apply(
            lambda x: get_average_temperature_during_heating_season(x, t_threshold=15)
        )
    )
    # accumulated difference between internal and external temperature
    # units ([K]-[K]) * [days/year]
    temperature_factor = (t_threshold - temperature_average_d_heat) * d_heat * 1 / 365

    return d_heat, temperature_factor


# windows ---------------------------------------------------------------
def window_limit(l, window_assumptions):  # noqa: E741
    """
    Define limit u value from which on window is retrofitted.
    """
    m = (
        (window_assumptions.diff()["u_limit"] / window_assumptions.diff()["strength"])
        .dropna()
        .iloc[0]
    )
    a = window_assumptions["u_limit"][0] - m * window_assumptions["strength"][0]
    return m * l + a


def u_retro_window(l, window_assumptions):  # noqa: E741
    """
    Define retrofitting value depending on renovation strength.
    """
    m = (
        (window_assumptions.diff()["u_value"] / window_assumptions.diff()["strength"])
        .dropna()
        .iloc[0]
    )
    a = window_assumptions["u_value"][0] - m * window_assumptions["strength"][0]
    return max(m * l + a, 0.8)


def window_cost(u, cost_retro, window_assumptions):  # noqa: E741
    """
    Get costs for new windows depending on u value.
    """
    m = (
        (window_assumptions.diff()["cost"] / window_assumptions.diff()["u_value"])
        .dropna()
        .iloc[0]
    )
    a = window_assumptions["cost"][0] - m * window_assumptions["u_value"][0]
    window_cost = m * u + a
    if annualise_cost:
        window_cost = (
            window_cost
            * interest_rate
            / (1 - (1 + interest_rate) ** -cost_retro.loc["Window", "life_time"])
        )
    return window_cost


def calculate_costs(u_values, l, cost_retro, window_assumptions):  # noqa: E741
    """
    Returns costs for a given retrofitting strength weighted by the average
    surface/volume ratio of the component for each building type.
    """
    return u_values.apply(
        lambda x: (
            (
                cost_retro.loc[x.name[3], "cost_var"]
                * 100
                * float(l)
                * l_weight.loc[x.name[3]].iloc[0]
                + cost_retro.loc[x.name[3], "cost_fix"]
            )
            * x.A_element
            / x.A_C_Ref
            if x.name[3] != "Window"
            else (
                (
                    (
                        window_cost(x[f"new_U_{l}"], cost_retro, window_assumptions)
                        * x.A_element
                    )
                    / x.A_C_Ref
                )
                if x.value > window_limit(float(l), window_assumptions)
                else 0
            )
        ),
        axis=1,
    )


def calculate_new_u(u_values, l, l_weight, window_assumptions, k=0.035):  # noqa: E741
    """
    Calculate U-values after building retrofitting, depending on the old
    U-values (u_values). This is for simple insulation measuers, adding an
    additional layer of insulation.

    They depend for the components Roof, Wall, Floor on the additional
    insulation thickness (l), and the weighting for the corresponding
    component (l_weight).

    Windows are renovated to new ones with U-value (function: u_retro_window(l))
    only if the are worse insulated than a certain limit value
    (function: window_limit).

    Parameters
    ----------
    u_values: pd.DataFrame
    l: string
    l_weight: pd.DataFrame (component, weight)
    k: thermal conductivity
    """
    return u_values.apply(
        lambda x: (
            k / ((k / x.value) + (float(l) * l_weight.loc[x.name[3]]))
            if x.name[3] != "Window"
            else (
                min(x.value, u_retro_window(float(l), window_assumptions))
                if x.value > window_limit(float(l), window_assumptions)
                else x.value
            )
        ),
        axis=1,
    )


def map_tabula_to_hotmaps(df_tabula, df_hotmaps, column_prefix):
    """
    Maps tabula data to hotmaps data with wished column name prefix.

    Parameters
    ----------
    df_tabula : pd.Series
        tabula data with pd.MultiIndex
    df_hotmaps : pd.DataFrame
        dataframe with hotmaps pd.MultiIndex
    column_prefix : string
        column prefix to rename column names of df_tabula

    Returns
    -------
    pd.DataFrame (index=df_hotmaps.index)
        returns df_tabula with hotmaps index
    """
    values = df_tabula.unstack().reindex(
        df_hotmaps.rename(
            index=lambda x: "MFH" if x not in rename_sectors.values() else x, level=1
        ).index
    )
    values.columns = pd.MultiIndex.from_product([[column_prefix], values.columns])
    values.index = df_hotmaps.index
    return values


def get_solar_gains_per_year(window_area):
    """
    Returns solar heat gains during heating season in [kWh/a] depending on the
    window area [m^2] of the building, assuming a equal distributed window
    orientation (east, south, north, west)
    """
    return sum(
        external_shading
        * frame_area_fraction
        * non_perpendicular
        * 0.25
        * window_area
        * solar_global_radiation
    )


def map_to_lstrength(l_strength, df):
    """
    Renames column names from a pandas dataframe to map tabula retrofitting
    strengths [2 = moderate, 3 = ambitious] to l_strength.
    """
    middle = len(l_strength) // 2
    map_to_l = pd.MultiIndex.from_arrays(
        [middle * [2] + len(l_strength[middle:]) * [3], l_strength]
    )
    l_strength_df = (
        df.stack(-2)
        .reindex(map_to_l, axis=1, level=0)
        .droplevel(0, axis=1)
        .unstack()
        .swaplevel(axis=1)
        .dropna(axis=1)
    )

    return pd.concat([df.drop([2, 3], axis=1, level=1), l_strength_df], axis=1)


def calculate_heat_losses(u_values, data_tabula, l_strength, temperature_factor):
    """
    Calculates total annual heat losses Q_ht for different insulation
    thicknesses (l_strength), depending on current insulation state (u_values),
    standard building topologies and air ventilation from TABULA (data_tabula)
    and the accumulated difference between internal and external temperature
    during the heating season (temperature_factor).

    Total annual heat losses Q_ht constitute from losses by:
        (1) transmission (H_tr_e)
        (2) thermal bridges (H_tb)
        (3) ventilation (H_ve)
    weighted by  a factor (F_red_temp) which is taken account for non-uniform heating
    and the temperature factor of the heating season

     Q_ht [W/m^2] = (H_tr_e + H_tb + H_ve) [W/m^2K] * F_red_temp * temperature_factor [K]

     returns Q_ht as pd.DataFrame(index=['country_code', 'subsector', 'bage'],
                          columns=[current (1.) + retrofitted (l_strength)])
    """
    #  (1) by transmission
    # calculate new U values of building elements due to additional insulation
    for l in l_strength:  # noqa: E741
        u_values[f"new_U_{l}"] = calculate_new_u(
            u_values, l, l_weight, window_assumptions
        )
    # surface area of building components [m^2]
    area_element = (
        data_tabula[[f"A_{e}" for e in u_values.index.levels[3]]]
        .rename(columns=lambda x: x[2:])
        .stack()
        .unstack(-2)
        .stack()
    )
    u_values["A_element"] = map_tabula_to_hotmaps(
        area_element, u_values, "A_element"
    ).xs(1, level=1, axis=1)

    # heat transfer H_tr_e [W/m^2K] through building element
    # U_e * A_e / A_C_Ref
    columns = ["value"] + [f"new_U_{l}" for l in l_strength]
    heat_transfer = pd.concat(
        [u_values[columns].mul(u_values.A_element, axis=0), u_values.A_element], axis=1
    )
    # get real subsector back in index
    heat_transfer.index = u_values.index
    heat_transfer = heat_transfer.groupby(level=[0, 1, 2]).sum()

    # rename columns of heat transfer H_tr_e [W/K] and envelope surface A_envelope [m^2]
    heat_transfer.rename(
        columns={
            "A_element": "A_envelope",
        },
        inplace=True,
    )

    # map reference area
    heat_transfer["A_C_Ref"] = map_tabula_to_hotmaps(
        data_tabula.A_C_Ref, heat_transfer, "A_C_Ref"
    ).xs(1.0, level=1, axis=1)
    u_values["A_C_Ref"] = map_tabula_to_hotmaps(
        data_tabula.A_C_Ref, u_values, "A_C_Ref"
    ).xs(1.0, level=1, axis=1)

    # get heat transfer  by transmission through building element [W/(m^2K)]
    heat_transfer_perm2 = heat_transfer[columns].div(heat_transfer.A_C_Ref, axis=0)
    heat_transfer_perm2.columns = pd.MultiIndex.from_product(
        [["H_tr_e"], [1.0] + l_strength]
    )

    # (2) heat transfer by thermal bridges H_tb [W/(m^2K)]
    # H_tb = delta_U [W/(m^2K)]* A_envelope [m^2] / A_C_Ref [m^2]
    H_tb_tabula = (
        data_tabula.delta_U_ThermalBridging
        * data_tabula.A_envelope
        / data_tabula.A_C_Ref
    )

    heat_transfer_perm2 = pd.concat(
        [
            heat_transfer_perm2,
            map_tabula_to_hotmaps(H_tb_tabula, heat_transfer_perm2, "H_tb"),
        ],
        axis=1,
    )

    # (3) by ventilation H_ve [W/(m²K)]
    # = c_p_air [Wh/(m^3K)] * (n_air_use + n_air_infilitraion) [1/h] * h_room [m]
    H_ve_tabula = (
        (data_tabula.n_air_infiltration + data_tabula.n_air_use) * c_p_air * h_room
    )
    heat_transfer_perm2 = pd.concat(
        [
            heat_transfer_perm2,
            map_tabula_to_hotmaps(H_ve_tabula, heat_transfer_perm2, "H_ve"),
        ],
        axis=1,
    )

    # F_red_temp factor which is taken account for non-uniform heating e.g.
    # lower heating/switch point during night times/weekends
    # effect is significant for buildings with poor insulation
    # for well insulated buildings/passive houses it has nearly no effect
    # based on tabula values depending on the building type
    F_red_temp = map_tabula_to_hotmaps(
        data_tabula.F_red_temp, heat_transfer_perm2, "F_red_temp"
    )
    # total heat transfer Q_ht [W/m^2] =
    # (H_tr_e + H_tb + H_ve) [W/m^2K] * F_red_temp * temperature_factor [K]
    # temperature_factor = (t_threshold - temperature_average_d_heat) * d_heat * 1/365
    heat_transfer_perm2 = map_to_lstrength(l_strength, heat_transfer_perm2)
    F_red_temp = map_to_lstrength(l_strength, F_red_temp)

    Q_ht = (
        heat_transfer_perm2.T.groupby(level=1)
        .sum()
        .T.mul(F_red_temp.droplevel(0, axis=1))
        .mul(temperature_factor.reindex(heat_transfer_perm2.index, level=0), axis=0)
    )

    return Q_ht, heat_transfer_perm2


def calculate_heat_gains(data_tabula, heat_transfer_perm2, d_heat):
    """
    Calculates heat gains Q_gain [W/m^2], which consititure from gains by:

    (1) solar radiation (2) internal heat gains
    """
    # (1) by solar radiation H_solar [W/m^2]
    # solar radiation [kWhm^2/a] / A_C_Ref [m^2] *1e3[1/k] / 8760 [a/h]
    H_solar = (
        data_tabula.A_Window.apply(lambda x: get_solar_gains_per_year(x))
        / data_tabula.A_C_Ref
        * 1e3
        / 8760
    )

    Q_gain = map_tabula_to_hotmaps(H_solar, heat_transfer_perm2, "H_solar").xs(
        1.0, level=1, axis=1
    )

    # (2) by  internal H_int
    # phi [W/m^2] * d_heat [d/a] * 1/365 [a/d] -> W/m^2
    Q_gain["H_int"] = (phi_int * d_heat * 1 / 365).reindex(
        index=heat_transfer_perm2.index, level=0
    )

    return Q_gain


def calculate_gain_utilisation_factor(heat_transfer_perm2, Q_ht, Q_gain):
    """
    Calculates gain utilisation factor nu.
    """
    # time constant of the building tau [h] = c_m [Wh/(m^2K)] * 1 /(H_tr_e+H_tb*H_ve) [m^2 K /W]
    tau = c_m / heat_transfer_perm2.T.groupby(axis=1).sum().T
    alpha = alpha_H_0 + (tau / tau_H_0)
    # heat balance ratio
    gamma = (1 / Q_ht).mul(Q_gain.sum(axis=1), axis=0)
    return (1 - gamma**alpha) / (1 - gamma ** (alpha + 1))


def calculate_space_heat_savings(
    u_values, data_tabula, l_strength, temperature_factor, d_heat
):
    """
    Calculates space heat savings (dE_space [per unit of unrefurbished state])
    through retrofitting of the thermal envelope by additional insulation
    material (l_strength[m])
    """
    # heat losses Q_ht [W/m^2]
    Q_ht, heat_transfer_perm2 = calculate_heat_losses(
        u_values, data_tabula, l_strength, temperature_factor
    )
    # heat gains Q_gain [W/m^2]
    Q_gain = calculate_heat_gains(data_tabula, heat_transfer_perm2, d_heat)

    # calculate gain utilisation factor nu [dimensionless]
    nu = calculate_gain_utilisation_factor(heat_transfer_perm2, Q_ht, Q_gain)

    # total space heating demand E_space
    E_space = Q_ht - nu.mul(Q_gain.sum(axis=1), axis=0)
    dE_space = E_space.div(E_space[1.0], axis=0).iloc[:, 1:]
    dE_space.columns = pd.MultiIndex.from_product([["dE"], l_strength])

    return dE_space


def calculate_retro_costs(u_values, l_strength, cost_retro):
    """
    Returns costs of different retrofitting measures.
    """
    costs = pd.concat(
        [
            calculate_costs(u_values, l, cost_retro, window_assumptions).rename(l)
            for l in l_strength
        ],
        axis=1,
    )

    # energy and costs per country, sector, subsector and year
    cost_tot = costs.groupby(level=["country_code", "subsector", "bage"]).sum()
    cost_tot.columns = pd.MultiIndex.from_product([["cost"], cost_tot.columns])

    return cost_tot


def sample_dE_costs_area(
    area, area_tot, costs, dE_space, countries, construction_index, tax_weighting
):
    """
    Bring costs and energy savings together, fill area and costs per energy
    savings for missing countries, weight costs, determine "moderate" and
    "ambitious" retrofitting.
    """
    sub_to_sector_dict = (
        area.reset_index()
        .replace(rename_sectors)
        .set_index("subsector")["sector"]
        .to_dict()
    )

    area_reordered = (
        (
            area.rename(index=country_iso_dic, level=0)
            .rename(index=rename_sectors, level=2)
            .reset_index()
        )
        # if uncommented, leads to the second `country_code` column
        # .rename(columns={"country": "country_code"})
        .set_index(["country_code", "subsector", "bage"])
    )

    cost_dE = (
        pd.concat([costs, dE_space], axis=1)
        .mul(area_reordered.weight, axis=0)
        .rename(sub_to_sector_dict, level=1)
        .groupby(level=[0, 1])
        .sum()
    )

    # map missing countries
    for ct in set(countries).difference(cost_dE.index.levels[0]):
        averaged_data = (
            cost_dE.reindex(index=map_for_missings[ct], level=0)
            .groupby(level=1)
            .mean()
            .set_index(pd.MultiIndex.from_product([[ct], cost_dE.index.levels[1]]))
        )
        cost_dE = pd.concat([cost_dE, averaged_data])

    # weights costs after construction index
    if construction_index:
        for ct in list(map_for_missings.keys() - cost_w.index):
            cost_w.loc[ct] = cost_w.reindex(index=map_for_missings[ct]).mean()
        cost_dE.cost = cost_dE.cost.mul(cost_w, level=0, axis=0)

    # weights cost depending on country taxes
    if tax_weighting:
        for ct in list(map_for_missings.keys() - tax_w.index):
            tax_w[ct] = tax_w.reindex(index=map_for_missings[ct]).mean()
        cost_dE.cost = cost_dE.cost.mul(tax_w, level=0, axis=0)

    # drop not considered countries
    cost_dE = cost_dE.reindex(countries, level=0)
    # get share of residential and service floor area
    sec_w = area_tot.div(area_tot.groupby(level=0).transform("sum"))
    # get the total cost-energy-savings weight by sector area
    tot = (
        # sec_w has columns "estimated" and "value"
        cost_dE.mul(sec_w.value, axis=0)
        # for some reasons names of the levels were lost somewhere
        # .groupby(level="country_code")
        .groupby(level=0)
        .sum()
        .set_index(pd.MultiIndex.from_product([cost_dE.index.unique(level=0), ["tot"]]))
    )
    cost_dE = pd.concat([cost_dE, tot]).unstack().stack()

    summed_area = pd.DataFrame(area_tot.groupby(level=0).sum()).set_index(
        pd.MultiIndex.from_product([area_tot.index.unique(level=0), ["tot"]])
    )
    area_tot = pd.concat([area_tot, summed_area]).unstack().stack()

    cost_per_saving = cost_dE["cost"] / (
        1 - cost_dE["dE"]
    )  # .diff(axis=1).dropna(axis=1)

    moderate_min = cost_per_saving.idxmin(axis=1)
    moderate_dE_cost = pd.concat(
        [cost_dE.loc[i].xs(moderate_min.loc[i], level=1) for i in moderate_min.index],
        axis=1,
    ).T
    moderate_dE_cost.columns = pd.MultiIndex.from_product(
        [moderate_dE_cost.columns, ["moderate"]]
    )

    ambitious_dE_cost = cost_dE.xs("0.26", level=1, axis=1)
    ambitious_dE_cost.columns = pd.MultiIndex.from_product(
        [ambitious_dE_cost.columns, ["ambitious"]]
    )

    cost_dE_new = pd.concat([moderate_dE_cost, ambitious_dE_cost], axis=1)

    return cost_dE_new, area_tot


# %% --- MAIN --------------------------------------------------------------
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_retro_cost",
            simpl="",
            clusters=48,
            ll="v1.0",
            sector_opts="Co2L0-168H-T-H-B-I-solar3-dist1",
        )

    #  ********  config  *********************************************************

    retro_opts = snakemake.params.retrofitting
    interest_rate = retro_opts["interest_rate"]
    annualise_cost = retro_opts["annualise_cost"]  # annualise the investment costs
    tax_weighting = retro_opts[
        "tax_weighting"
    ]  # weight costs depending on taxes in countries
    construction_index = retro_opts[
        "construction_index"
    ]  # weight costs depending on labour/material costs per ct

    # mapping missing countries by neighbours
    map_for_missings = {
        "AL": ["BG", "RO", "GR"],
        "BA": ["HR"],
        "RS": ["BG", "RO", "HR", "HU"],
        "MK": ["BG", "GR"],
        "ME": ["BA", "AL", "RS", "HR"],
        "CH": ["SE", "DE"],
        "NO": ["SE"],
    }

    #   (1) prepare data **********************************************************

    # building stock data -----------------------------------------------------
    # hotmaps u_values, heated floor areas per sector
    u_values, country_iso_dic, countries, area_tot, area = prepare_building_stock_data()
    # building topology, thermal bridges, ventilation losses
    data_tabula = prepare_building_topology(u_values)
    # costs for retrofitting -------------------------------------------------
    cost_retro, window_assumptions, cost_w, tax_w = prepare_cost_retro(country_iso_dic)
    # temperature dependent parameters
    d_heat, temperature_factor = prepare_temperature_data()

    #  (2) space heat savings ****************************************************
    dE_space = calculate_space_heat_savings(
        u_values, data_tabula, l_strength, temperature_factor, d_heat
    )

    #  (3) costs *****************************************************************
    costs = calculate_retro_costs(u_values, l_strength, cost_retro)

    #  (4) cost-dE and area per sector *******************************************
    cost_dE, area_tot = sample_dE_costs_area(
        area, area_tot, costs, dE_space, countries, construction_index, tax_weighting
    )

    #   save *********************************************************************
    cost_dE.to_csv(snakemake.output.retro_cost)
    area_tot.to_csv(snakemake.output.floor_area)
