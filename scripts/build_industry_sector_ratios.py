# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build best case specific energy consumption by carrier and category.

Description
-------

This script uses the `JRC-IDEES <https://joint-research-centre.ec.europa.eu/potencia-policy-oriented-tool-energy-and-climate-change-impact-assessment/jrc-idees_en>` data to calculate an EU28 average specific energy consumption by carrier and industries.
The industries are according to the rule `industrial_production_per_country <https://pypsa-eur.readthedocs.io/en/latest/sector.html#module-build_industrial_production_per_country>`.

The following carriers are considered:
- elec
- coal
- coke
- biomass
- methane
- hydrogen
- heat
- naphtha
- process emission
- process emission from feedstock
- (ammonia)

If the `config["industry"]["ammonia"] <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#industry>` is set to true the ammonia demand is not converted to hydrogen and electricity but is considered as a separate carrier.

The unit of the specific energy consumption is MWh/t material and tCO2/t material for process emissions.
"""

import logging

import country_converter as coco
import pandas as pd
from _helpers import configure_logging, mute_print, set_scenario_config

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()

# GWh/ktoe OR MWh/toe
toe_to_MWh = 11.630

eu27 = cc.EU27as("ISO2").ISO2.tolist()


sheet_names = {
    "Iron and steel": "ISI",
    "Chemicals Industry": "CHI",
    "Non-metallic mineral products": "NMM",
    "Pulp, paper and printing": "PPA",
    "Food, beverages and tobacco": "FBT",
    "Non Ferrous Metals": "NFM",
    "Transport equipment": "TRE",
    "Machinery equipment": "MAE",
    "Textiles and leather": "TEL",
    "Wood and wood products": "WWP",
    "Other industrial sectors": "OIS",
}


index = [
    "elec",
    "coal",
    "coke",
    "biomass",
    "methane",
    "hydrogen",
    "heat",
    "naphtha",
    "ammonia",
    "methanol",
    "process emission",
    "process emission from feedstock",
]


def load_idees_data(sector, country="EU27"):
    suffixes = {"out": "", "fec": "_fec", "ued": "_ued", "emi": "_emi"}
    sheets = {k: sheet_names[sector] + v for k, v in suffixes.items()}

    def usecols(x):
        return isinstance(x, str) or x == year

    with mute_print():
        idees = pd.read_excel(
            f"{snakemake.input.idees}/{country}/JRC-IDEES-2021_Industry_{country}.xlsx",
            sheet_name=list(sheets.values()),
            index_col=0,
            header=0,
            usecols=usecols,
        )

    for k, v in sheets.items():
        idees[k] = idees.pop(v).squeeze()
        idees[k] = idees[k][year]

    return idees


def iron_and_steel():
    """
    This function calculates the energy consumption and emissions for different
    approaches to producing iron and steel. The two primary approaches are
    integrated steelworks and electric arc furnaces (EAF). The function assumes
    that integrated steelworks will be replaced entirely by electric arc
    furnaces due to their higher efficiency and greater reliance on
    electricity.

    Returns:
        pd.DataFrame: A DataFrame containing the energy consumption (in MWh/t material)
                      and process emissions (in tCO2/t material) for different steel
                      production approaches, including electric arc, DRI + electric arc,
                      and integrated steelworks.
    """

    sector = "Iron and steel"
    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    ## Electric arc

    sector = "Electric arc"

    df[sector] = 0.0

    s_fec = idees["fec"][52:68]
    assert s_fec.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.at["elec", sector] += s_fec.loc[sel].sum()

    df.at["heat", sector] += s_fec.loc["Low-enthalpy heat"]

    subsector = "Steel: Smelters"
    s_fec = idees["fec"][63:68]
    s_ued = idees["ued"][63:68]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    # efficiency changes due to transforming all the smelters into methane
    key = "Natural gas and biogas"
    eff_met = s_ued.loc[key] / s_fec.loc[key]

    df.at["methane", sector] += s_ued[subsector] / eff_met

    subsector = "Steel: Electric arc"
    s_fec = idees["fec"][69:70]
    assert s_fec.index[0] == subsector

    df.at["elec", sector] += s_fec[subsector]

    subsector = "Steel: Furnaces, refining and rolling"
    s_fec = idees["fec"][70:77]
    s_ued = idees["ued"][70:77]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "Steel: Furnaces, refining and rolling - Electric"
    eff = s_ued[key] / s_fec[key]

    # assume fully electrified, other processes scaled by used energy
    df.at["elec", sector] += s_ued[subsector] / eff

    subsector = "Steel: Product finishing"
    s_fec = idees["fec"][77:95]
    s_ued = idees["ued"][77:95]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "Steel: Product finishing - Electric"
    eff = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.at["elec", sector] += s_ued[subsector] / eff

    # Process emissions (per physical output)

    s_emi = idees["emi"][52:95]
    assert s_emi.index[0] == sector

    s_out = idees["out"][7:8]
    assert s_out.index[0] == sector

    # tCO2/t material
    df.loc["process emission", sector] += s_emi["Process emissions"] / s_out[sector]

    # final energy consumption MWh/t material
    sel = ["elec", "heat", "methane"]
    df.loc[sel, sector] = df.loc[sel, sector] * toe_to_MWh / s_out[sector]

    ## DRI + Electric arc
    # For primary route: DRI with H2 + EAF

    sector = "DRI + Electric arc"

    df[sector] = df["Electric arc"]

    # add H2 consumption for DRI at 1.7 MWh H2 /ton steel
    df.at["hydrogen", sector] = params["H2_DRI"]

    # add electricity consumption in DRI shaft (0.322 MWh/tSl)
    df.at["elec", sector] += params["elec_DRI"]

    ## Integrated steelworks
    # could be used in combination with CCS)
    # Assume existing fuels are kept, except for furnaces, refining, rolling, finishing
    # Ignore 'derived gases' since these are top gases from furnaces

    sector = "Integrated steelworks"

    df[sector] = 0.0

    s_fec = idees["fec"][3:50]
    assert s_fec.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    subsector = "Steel: Sinter/Pellet-making"

    s_fec = idees["fec"][14:20]
    s_ued = idees["ued"][14:20]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    df.loc["elec", sector] += s_fec["Electricity"]

    sel = ["Natural gas and biogas", "Fuel oil"]
    df.loc["methane", sector] += s_fec[sel].sum()

    df.loc["coal", sector] += s_fec["Solids"]

    subsector = "Steel: Blast /Basic oxygen furnace"

    s_fec = idees["fec"][20:26]
    s_ued = idees["ued"][20:26]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    sel = ["Natural gas and biogas", "Fuel oil"]
    df.loc["methane", sector] += s_fec[sel].sum()

    df.loc["coal", sector] += s_fec["Solids"]

    df.loc["coke", sector] = s_fec["Coke"]

    subsector = "Steel: Furnaces, refining and rolling"

    s_fec = idees["fec"][26:33]
    s_ued = idees["ued"][26:33]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "Steel: Furnaces, refining and rolling - Electric"
    eff = s_ued[key] / s_fec[key]

    # assume fully electrified, other processes scaled by used energy
    df.loc["elec", sector] += s_ued[subsector] / eff

    subsector = "Steel: Product finishing"

    s_fec = idees["fec"][33:50]
    s_ued = idees["ued"][33:50]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "Steel: Product finishing - Electric"
    eff = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.loc["elec", sector] += s_ued[subsector] / eff

    # Process emissions (per physical output)

    s_emi = idees["emi"][3:51]
    assert s_emi.index[0] == sector

    s_out = idees["out"][6:7]
    assert s_out.index[0] == sector

    # tCO2/t material
    df.loc["process emission", sector] = s_emi["Process emissions"] / s_out[sector]

    # final energy consumption MWh/t material
    sel = ["elec", "heat", "methane", "coke", "coal"]
    df.loc[sel, sector] = df.loc[sel, sector] * toe_to_MWh / s_out[sector]

    return df


def chemicals_industry():
    """
    This function calculates the energy consumption and emissions for the
    chemicals industry, focusing on various subsectors such as basic chemicals,
    steam processing, furnaces, and process cooling. The function also accounts
    for specific processes in ammonia, chlorine, methanol production, and other
    chemicals.

    Returns:
        pd.DataFrame: A DataFrame containing the energy consumption (in MWh/t material)
                      and process emissions (in tCO2/t material) for various subsectors
                      within the chemicals industry.
    """
    sector = "Chemicals Industry"
    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    # Basic chemicals

    sector = "Basic chemicals"

    df[sector] = 0.0

    s_fec = idees["fec"][3:9]
    assert s_fec.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    subsector = "Chemicals: Feedstock (energy used as raw material)"
    # There are Solids, Refinery gas, LPG, Diesel oil, Fuel oil,
    # Other liquids, Naphtha, Natural gas for feedstock.
    # Naphta represents 47%, methane 17%. LPG (18%) solids, refinery gas,
    # diesel oil, Fuel oils and other liquids are assimilated to Naphtha

    s_fec = idees["fec"][14:23]
    assert s_fec.index[0] == subsector

    df.loc["naphtha", sector] += s_fec["Naphtha"]

    df.loc["methane", sector] += s_fec["Natural gas"]

    # LPG and other feedstock materials are assimilated to naphtha
    # since they will be produced through Fischer-Tropsch process
    sel = [
        "Solids",
        "Refinery gas",
        "LPG",
        "Diesel oil",
        "Fuel oil",
        "Other liquids",
    ]
    df.loc["naphtha", sector] += s_fec[sel].sum()

    subsector = "Chemicals: Steam processing"
    # All the final energy consumption in the steam processing is
    # converted to methane, since we need >1000 C temperatures here.
    # The current efficiency of methane is assumed in the conversion.

    s_fec = idees["fec"][23:34]
    s_ued = idees["ued"][23:34]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    # efficiency of natural gas
    eff_ch4 = s_ued["Natural gas and biogas"] / s_fec["Natural gas and biogas"]

    # replace all fec by methane
    df.loc["methane", sector] += s_ued[subsector] / eff_ch4

    subsector = "Chemicals: Furnaces"

    s_fec = idees["fec"][34:42]
    s_ued = idees["ued"][34:42]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    # efficiency of electrification
    key = "Chemicals: Furnaces - Electric"
    eff_elec = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.loc["elec", sector] += s_ued[subsector] / eff_elec

    subsector = "Chemicals: Process cooling"

    s_fec = idees["fec"][42:56]
    s_ued = idees["ued"][42:56]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "Chemicals: Process cooling - Electric"
    eff_elec = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.loc["elec", sector] += s_ued[subsector] / eff_elec

    subsector = "Chemicals: Generic electric process"

    s_fec = idees["fec"][56:57]
    assert s_fec.index[0] == subsector

    df.loc["elec", sector] += s_fec[subsector]

    # Process emissions

    # Correct everything by subtracting 2019's ammonia demand and
    # putting in ammonia demand for H2 and electricity separately

    s_emi = idees["emi"][3:58]
    assert s_emi.index[0] == sector

    # convert from MtHVC/a to ktHVC/a
    s_out = params["HVC_production_today"] * 1e3

    # tCO2/t material
    df.loc["process emission", sector] += (
        s_emi["Process emissions"]
        - params["petrochemical_process_emissions"] * 1e3
        - params["NH3_process_emissions"] * 1e3
    ) / s_out

    # emissions originating from feedstock, could be non-fossil origin
    # tCO2/t material
    df.loc["process emission from feedstock", sector] += (
        params["petrochemical_process_emissions"] * 1e3
    ) / s_out

    # convert from ktoe/a to GWh/a
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] *= toe_to_MWh

    # subtract ammonia energy demand (in ktNH3/a)
    ammonia = pd.read_csv(snakemake.input.ammonia_production, index_col=0)
    ammonia_total = ammonia.loc[
        ammonia.index.intersection(eu27), str(max(2018, year))
    ].sum()
    df.loc["methane", sector] -= ammonia_total * params["MWh_CH4_per_tNH3_SMR"]
    df.loc["elec", sector] -= ammonia_total * params["MWh_elec_per_tNH3_SMR"]

    # subtract chlorine demand (in MtCl/a)
    chlorine_total = params["chlorine_production_today"]
    df.loc["hydrogen", sector] -= chlorine_total * params["MWh_H2_per_tCl"] * 1e3
    df.loc["elec", sector] -= chlorine_total * params["MWh_elec_per_tCl"] * 1e3

    # subtract methanol demand (in MtMeOH/a)
    methanol_total = params["methanol_production_today"]
    df.loc["methane", sector] -= methanol_total * params["MWh_CH4_per_tMeOH"] * 1e3
    df.loc["elec", sector] -= methanol_total * params["MWh_elec_per_tMeOH"] * 1e3

    # MWh/t material
    df.loc[sources, sector] = df.loc[sources, sector] / s_out

    df.rename(columns={sector: "HVC"}, inplace=True)

    # HVC mechanical recycling

    sector = "HVC (mechanical recycling)"
    df[sector] = 0.0
    df.loc["elec", sector] = params["MWh_elec_per_tHVC_mechanical_recycling"]

    # HVC chemical recycling

    sector = "HVC (chemical recycling)"
    df[sector] = 0.0
    df.loc["elec", sector] = params["MWh_elec_per_tHVC_chemical_recycling"]

    # Ammonia

    sector = "Ammonia"
    df[sector] = 0.0
    if snakemake.params.ammonia:
        df.loc["ammonia", sector] = params["MWh_NH3_per_tNH3"]
    else:
        df.loc["hydrogen", sector] = params["MWh_H2_per_tNH3_electrolysis"]
        df.loc["elec", sector] = params["MWh_elec_per_tNH3_electrolysis"]

    # Chlorine

    sector = "Chlorine"
    df[sector] = 0.0
    df.loc["hydrogen", sector] = params["MWh_H2_per_tCl"]
    df.loc["elec", sector] = params["MWh_elec_per_tCl"]

    # Methanol

    sector = "Methanol"
    df[sector] = 0.0
    df.loc["methanol", sector] = params["MWh_MeOH_per_tMeOH"]

    # Other chemicals

    sector = "Other chemicals"

    df[sector] = 0.0

    s_fec = idees["fec"][59:65]
    assert s_fec.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    subsector = "Chemicals: High-enthalpy heat processing"

    s_fec = idees["fec"][70:83]
    s_ued = idees["ued"][70:83]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "High-enthalpy heat processing - Electric (microwave)"
    eff_elec = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.loc["elec", sector] += s_ued[subsector] / eff_elec

    subsector = "Chemicals: Furnaces"

    s_fec = idees["fec"][83:92]
    s_ued = idees["ued"][83:92]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "Chemicals: Furnaces - Electric"
    eff_elec = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.loc["elec", sector] += s_ued[subsector] / eff_elec

    subsector = "Chemicals: Process cooling"

    s_fec = idees["fec"][91:105]
    s_ued = idees["ued"][91:105]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "Chemicals: Process cooling - Electric"
    eff = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.loc["elec", sector] += s_ued[subsector] / eff

    subsector = "Chemicals: Generic electric process"

    s_fec = idees["fec"][105:106]
    assert s_fec.index[0] == subsector

    df.loc["elec", sector] += s_fec[subsector]

    # Process emissions

    s_emi = idees["emi"][59:107]
    s_out = idees["out"][9:10]
    assert s_emi.index[0] == sector
    assert sector in str(s_out.index)

    # tCO2/t material
    df.loc["process emission", sector] += s_emi["Process emissions"] / s_out.values

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = df.loc[sources, sector] * toe_to_MWh / s_out.values

    # Pharmaceutical products

    sector = "Pharmaceutical products etc."

    df[sector] = 0.0

    s_fec = idees["fec"][108:114]
    assert s_fec.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    subsector = "Chemicals: High-enthalpy heat processing"

    s_fec = idees["fec"][119:132]
    s_ued = idees["ued"][119:132]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "High-enthalpy heat processing - Electric (microwave)"
    eff_elec = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.loc["elec", sector] += s_ued[subsector] / eff_elec

    subsector = "Chemicals: Furnaces"

    s_fec = idees["fec"][132:140]
    s_ued = idees["ued"][132:140]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "Chemicals: Furnaces - Electric"
    eff = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.loc["elec", sector] += s_ued[subsector] / eff

    subsector = "Chemicals: Process cooling"

    s_fec = idees["fec"][140:154]
    s_ued = idees["ued"][140:154]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    key = "Chemicals: Process cooling - Electric"
    eff_elec = s_ued[key] / s_fec[key]

    # assume fully electrified
    df.loc["elec", sector] += s_ued[subsector] / eff_elec

    subsector = "Chemicals: Generic electric process"

    s_fec = idees["fec"][154:155]
    s_out = idees["out"][10:11]
    assert s_fec.index[0] == subsector
    assert sector in str(s_out.index)

    df.loc["elec", sector] += s_fec[subsector]

    # tCO2/t material
    df.loc["process emission", sector] += 0.0

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = df.loc[sources, sector] * toe_to_MWh / s_out.values

    return df


def nonmetalic_mineral_products():
    """
    This function calculates the energy consumption and emissions for the non-
    metallic mineral products industry, focusing on three main sectors: cement,
    ceramics, and glass production. It takes into account the specific
    processes and their associated energy types and emissions.

    Returns:
        pd.DataFrame: A DataFrame containing the energy consumption (in MWh/t material)
                      and process emissions (in tCO2/t material) for the cement, ceramics,
                      and glass production sectors.
    """

    sector = "Non-metallic mineral products"
    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    # Cement

    # This sector has process-emissions.
    # Includes three subcategories:
    # (a) Grinding, milling of raw material,
    # (b) Pre-heating and pre-calcination,
    # (c) clinker production (kilns),
    # (d) Grinding, packaging.
    # (b)+(c) represent 94% of fec. So (a) is joined to (b) and (d) is joined to (c).
    # Temperatures above 1400C are required for processing limestone and sand into clinker.
    # Everything (except current electricity and heat consumption and existing biomass)
    # is transformed into methane for high T.

    sector = "Cement"

    df[sector] = 0.0

    s_fec = idees["fec"][3:25]
    s_ued = idees["ued"][3:25]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # pre-processing: keep existing elec and biomass, rest to methane
    df.loc["elec", sector] += s_fec["Cement: Grinding, milling of raw material"]
    df.loc["biomass", sector] += s_fec["Biomass and waste"]
    df.loc["methane", sector] += (
        s_fec["Cement: Pre-heating and pre-calcination"] - s_fec["Biomass and waste"]
    )

    subsector = "Cement: Clinker production (kilns)"

    s_fec = idees["fec"][23:32]
    s_ued = idees["ued"][23:32]
    assert s_fec.index[0] == subsector
    assert s_ued.index[0] == subsector

    df.loc["biomass", sector] += s_fec["Biomass and waste"]
    df.loc["methane", sector] += (
        s_fec["Cement: Clinker production (kilns)"] - s_fec["Biomass and waste"]
    )
    df.loc["elec", sector] += s_fec["Cement: Grinding, packaging and precasting"]

    # Process emissions

    # come from calcination of limestone to chemically reactive calcium oxide (lime).
    # Calcium carbonate -> lime + CO2
    # CaCO3  -> CaO + CO2

    s_emi = idees["emi"][3:45]
    assert s_emi.index[0] == sector

    s_out = idees["out"][7:8]
    assert sector in str(s_out.index)

    # tCO2/t material
    df.loc["process emission", sector] += s_emi["Process emissions"] / s_out.values

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = df.loc[sources, sector] * toe_to_MWh / s_out.values

    # Ceramics & other NMM

    # This sector has process emissions.
    # Includes four subcategories:
    # (a) Mixing of raw material,
    # (b) Drying and sintering of raw material,
    # (c) Primary production process,
    # (d) Product finishing.
    # (b) represents 65% of fec and (a) 4%. So (a) is joined to (b).
    # Everything is electrified

    sector = "Ceramics & other NMM"

    df[sector] = 0.0

    s_fec = idees["fec"][46:95]
    s_ued = idees["ued"][46:95]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification
    key = "Ceramics: Microwave drying and sintering"
    # the values are zero in new JRC-data -> assume here value from JRC-2015
    # eff_elec = s_ued[key] / s_fec[key]
    eff_elec = 11.6 / 26

    sel = [
        "Ceramics: Mixing of raw material",
        "Ceramics: Drying and sintering of raw material",
    ]
    df.loc["elec", sector] += s_ued[sel].sum() / eff_elec

    key = "Ceramics: Electric kiln"
    eff_elec = s_ued[key] / s_fec[key]

    df.loc["elec", sector] += s_ued["Ceramics: Primary production process"] / eff_elec

    key = "Ceramics: Electric furnace"
    eff_elec = s_ued[key] / s_fec[key]

    df.loc["elec", sector] += s_ued["Ceramics: Product finishing"] / eff_elec

    s_emi = idees["emi"][46:96]
    assert s_emi.index[0] == sector

    s_out = idees["out"][8:9]
    assert sector in str(s_out.index)

    # tCO2/t material
    df.loc["process emission", sector] += s_emi["Process emissions"] / s_out.values

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = df.loc[sources, sector] * toe_to_MWh / s_out.values

    # Glass production

    # This sector has process emissions.
    # Includes four subcategories:
    # (a) Melting tank
    # (b) Forming
    # (c) Annealing
    # (d) Finishing processes.
    # (a) represents 73%. (b), (d) are joined to (c).
    # Everything is electrified.

    sector = "Glass production"

    df[sector] = 0.0

    s_fec = idees["fec"][97:126]
    s_ued = idees["ued"][97:126]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification
    key = "Glass: Electric melting tank"
    eff_elec = s_ued[key] / s_fec[key]

    df.loc["elec", sector] += s_ued["Glass: Melting tank"] / eff_elec

    key = "Glass: Annealing - electric"
    eff_elec = s_ued[key] / s_fec[key]

    sel = ["Glass: Forming", "Glass: Annealing", "Glass: Finishing processes"]
    df.loc["elec", sector] += s_ued[sel].sum() / eff_elec

    s_emi = idees["emi"][97:127]
    assert s_emi.index[0] == sector

    s_out = idees["out"][9:10]
    assert sector in str(s_out.index)

    # tCO2/t material
    df.loc["process emission", sector] += s_emi["Process emissions"] / s_out.values

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = df.loc[sources, sector] * toe_to_MWh / s_out.values

    return df


def pulp_paper_printing():
    """
    Models the energy consumption for the pulp, paper, and printing sector,
    assuming complete electrification of all processes. This sector does not
    have any process emissions associated with it.

    Returns:
        pd.DataFrame: A DataFrame containing the energy consumption (in MWh/t material)
                      for the pulp, paper, and printing sector.
    """

    sector = "Pulp, paper and printing"
    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    # Pulp production

    # Includes three subcategories:
    # (a) Wood preparation, grinding;
    # (b) Pulping;
    # (c) Cleaning.
    #
    # (b) Pulping is either biomass or electric; left like this (dominated by biomass).
    # (a) Wood preparation, grinding and (c) Cleaning represent only 10% of their current
    # energy consumption is assumed to be electrified without any change in efficiency

    sector = "Pulp production"

    df[sector] = 0.0

    s_fec = idees["fec"][3:29]
    s_ued = idees["ued"][3:29]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Industry-specific
    sel = [
        "Pulp: Wood preparation, grinding",
        "Pulp: Cleaning",
        "Pulp: Pulping electric",
    ]
    df.loc["elec", sector] += s_fec[sel].sum()

    # Efficiency changes due to biomass
    eff_bio = s_ued["Biomass and waste"] / s_fec["Biomass and waste"]
    df.loc["biomass", sector] += s_ued["Pulp: Pulping thermal"] / eff_bio

    s_out = idees["out"][8:9]
    assert sector in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Pulp production (kt)"]
    )

    # Paper production

    # Includes three subcategories:
    # (a) Stock preparation;
    # (b) Paper machine;
    # (c) Product finishing.
    #
    # (b) Paper machine and (c) Product finishing are left electric
    # and thermal is moved to biomass. The efficiency is calculated
    # from the pulping process that is already biomass.
    #
    # (a) Stock preparation represents only 7% and its current energy
    # consumption is assumed to be electrified without any change in efficiency.

    sector = "Paper production"

    df[sector] = 0.0

    s_fec = idees["fec"][30:80]
    s_ued = idees["ued"][30:80]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Industry-specific
    df.loc["elec", sector] += s_fec["Paper: Stock preparation"]

    # add electricity from process that is already electrified
    df.loc["elec", sector] += s_fec["Paper: Paper machine - Electricity"]

    # add electricity from process that is already electrified
    df.loc["elec", sector] += s_fec["Paper: Product finishing - Electricity"]

    s_fec = idees["fec"][55:66]
    s_ued = idees["ued"][55:66]
    assert s_fec.index[0] == "Paper: Paper machine - Steam use"
    assert s_ued.index[0] == "Paper: Paper machine - Steam use"

    # Efficiency changes due to biomass
    eff_bio = s_ued["Biomass and waste"] / s_fec["Biomass and waste"]
    df.loc["biomass", sector] += s_ued["Paper: Paper machine - Steam use"] / eff_bio

    s_fec = idees["fec"][68:79]
    s_ued = idees["ued"][68:79]
    assert s_fec.index[0] == "Paper: Product finishing - Steam use"
    assert s_ued.index[0] == "Paper: Product finishing - Steam use"

    # Efficiency changes due to biomass
    eff_bio = s_ued["Biomass and waste"] / s_fec["Biomass and waste"]
    df.loc["biomass", sector] += s_ued["Paper: Product finishing - Steam use"] / eff_bio

    s_out = idees["out"][9:10]
    assert sector in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = df.loc[sources, sector] * toe_to_MWh / s_out.values

    # Printing and media reproduction

    # (a) Printing and publishing is assumed to be
    # electrified without any change in efficiency.

    sector = "Printing and media reproduction"

    df[sector] = 0.0

    s_fec = idees["fec"][81:93]
    s_ued = idees["ued"][81:93]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()
    df.loc["elec", sector] += s_ued[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]
    df.loc["heat", sector] += s_ued["Low-enthalpy heat"]

    # Industry-specific
    df.loc["elec", sector] += s_fec["Printing and publishing"]
    df.loc["elec", sector] += s_ued["Printing and publishing"]

    s_out = idees["out"][10:11]
    assert sector in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = df.loc[sources, sector] * toe_to_MWh / s_out.values

    return df


def food_beverages_tobacco():
    """
    Calculates the energy consumption for the food, beverages, and tobacco
    sector, assuming complete electrification of all processes. This sector
    does not have any process emissions associated with it.

    Returns:
        pd.DataFrame: A DataFrame containing the energy consumption (in MWh/t material)
                      for the food, beverages, and tobacco sector.
    """

    sector = "Food, beverages and tobacco"
    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    df[sector] = 0.0

    s_fec = idees["fec"][3:79]
    s_ued = idees["ued"][3:79]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification

    key = "Food: Direct Heat - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Food: Oven (direct heat)"] / eff_elec

    key = "Food: Process Heat - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Food: Specific process heat"] / eff_elec

    key = "Food: Electric drying"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Food: Drying"] / eff_elec

    key = "Food: Electric cooling"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += (
        s_ued["Food: Process cooling and refrigeration"] / eff_elec
    )

    # Steam processing goes all to biomass without change in efficiency
    df.loc["biomass", sector] += s_fec["Food: Steam processing"]

    # add electricity from process that is already electrified
    df.loc["elec", sector] += s_fec["Food: Electric machinery"]

    s_out = idees["out"][3:4]
    assert "Physical output" in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Physical output (index)"]
    )

    return df


def non_ferrous_metals():
    sector = "Non Ferrous Metals"
    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    # Alumina

    # High-enthalpy heat is converted to methane.
    # Process heat at T>500C is required here.
    # Refining is electrified.
    # There are no process emissions associated to Alumina manufacturing.

    sector = "Alumina production"

    df[sector] = 0.0

    s_fec = idees["fec"][3:31]
    s_ued = idees["ued"][3:31]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # High-enthalpy heat is transformed into methane

    s_fec = idees["fec"][14:25]
    s_ued = idees["ued"][14:25]
    assert s_fec.index[0] == "Alumina production: High-enthalpy heat"
    assert s_ued.index[0] == "Alumina production: High-enthalpy heat"

    eff_met = s_ued["Natural gas and biogas"] / s_fec["Natural gas and biogas"]
    df.loc["methane", sector] += (
        s_fec["Alumina production: High-enthalpy heat"] / eff_met
    )

    # Efficiency changes due to electrification

    s_fec = idees["fec"][25:31]
    s_ued = idees["ued"][25:31]
    assert s_fec.index[0] == "Alumina production: Refining"
    assert s_ued.index[0] == "Alumina production: Refining"

    eff_elec = s_ued["Electricity"] / s_fec["Electricity"]
    df.loc["elec", sector] += s_ued["Alumina production: Refining"] / eff_elec

    s_out = idees["out"][9:10]
    assert sector in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Alumina production (kt)"]
    )

    # Aluminium primary route

    # Production through the primary route is divided into 50% remains
    # as today and 50% is transformed into secondary route.

    sector = "Aluminium - primary production"

    df[sector] = 0.0

    s_fec = idees["fec"][32:68]
    s_ued = idees["ued"][32:68]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Add aluminium  electrolysis (smelting
    df.loc["elec", sector] += s_fec["Aluminium electrolysis (smelting)"]

    # Efficiency changes due to electrification
    key = "Aluminium processing - Electric"
    eff_elec = s_ued[key] / s_fec[key]

    key = "Aluminium processing  (metallurgy e.g. cast house, reheating)"
    df.loc["elec", sector] += s_ued[key] / eff_elec

    key = "Aluminium finishing - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Aluminium finishing"] / eff_elec

    s_emi = idees["emi"][32:69]
    assert s_emi.index[0] == sector

    s_out = idees["out"][11:12]
    assert sector in str(s_out.index)

    # tCO2/t material
    df.loc["process emission", sector] = (
        s_emi["Process emissions"] / s_out["Aluminium - primary production"]
    )

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Aluminium - primary production"]
    )

    # Aluminium secondary route

    # All is converted into secondary route fully electrified.

    sector = "Aluminium - secondary production"

    df[sector] = 0.0

    s_fec = idees["fec"][70:112]
    s_ued = idees["ued"][70:112]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification
    key = "Secondary aluminium - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    key = "Secondary aluminium (incl. pre-treatment, remelting)"
    df.loc["elec", sector] += s_ued[key] / eff_elec

    key = "Aluminium processing - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    key = "Aluminium processing  (metallurgy e.g. cast house, reheating)"
    df.loc["elec", sector] += s_ued[key] / eff_elec

    key = "Aluminium finishing - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Aluminium finishing"] / eff_elec

    s_out = idees["out"][12:13]
    assert sector in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Aluminium - secondary production"]
    )

    # Other non-ferrous metals

    sector = "Other non-ferrous metals"

    df[sector] = 0.0

    s_fec = idees["fec"][113:156]
    s_ued = idees["ued"][113:156]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification
    key = "Metal production - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Other Metals: production"] / eff_elec

    key = "Metal processing - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    key = "Metal processing  (metallurgy e.g. cast house, reheating)"
    df.loc["elec", sector] += s_ued[key] / eff_elec

    key = "Metal finishing - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Metal finishing"] / eff_elec

    s_emi = idees["emi"][113:157]
    assert s_emi.index[0] == sector

    s_out = idees["out"][13:14]
    assert sector in str(s_out.index)

    # tCO2/t material
    df.loc["process emission", sector] = (
        s_emi["Process emissions"] / s_out["Other non-ferrous metals (kt lead eq.)"]
    )

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector]
        * toe_to_MWh
        / s_out["Other non-ferrous metals (kt lead eq.)"]
    )

    return df


def transport_equipment():
    sector = "Transport equipment"
    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    df[sector] = 0.0

    s_fec = idees["fec"][3:46]
    s_ued = idees["ued"][3:46]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification
    key = "Trans. Eq.: Electric Foundries"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Trans. Eq.: Foundries"] / eff_elec

    key = "Trans. Eq.: Electric connection"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Trans. Eq.: Connection techniques"] / eff_elec

    key = "Trans. Eq.: Heat treatment - Electric"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Trans. Eq.: Heat treatment"] / eff_elec

    df.loc["elec", sector] += s_fec["Trans. Eq.: General machinery"]
    df.loc["elec", sector] += s_fec["Trans. Eq.: Product finishing"]

    # Steam processing is supplied with biomass
    eff_biomass = s_ued["Biomass and waste"] / s_fec["Biomass and waste"]
    df.loc["biomass", sector] += s_ued["Trans. Eq.: Steam processing"] / eff_biomass

    s_out = idees["out"][3:4]
    assert "Physical output" in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Physical output (index)"]
    )

    return df


def machinery_equipment():
    sector = "Machinery equipment"

    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    df[sector] = 0.0

    s_fec = idees["fec"][3:46]
    s_ued = idees["ued"][3:46]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification
    key = "Mach. Eq.: Electric Foundries"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Mach. Eq.: Foundries"] / eff_elec

    key = "Mach. Eq.: Electric connection"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Mach. Eq.: Connection techniques"] / eff_elec

    key = "Mach. Eq.: Heat treatment - Electric"
    eff_elec = s_ued[key] / s_fec[key]

    df.loc["elec", sector] += s_ued["Mach. Eq.: Heat treatment"] / eff_elec

    df.loc["elec", sector] += s_fec["Mach. Eq.: General machinery"]
    df.loc["elec", sector] += s_fec["Mach. Eq.: Product finishing"]

    # Steam processing is supplied with biomass
    eff_biomass = s_ued["Biomass and waste"] / s_fec["Biomass and waste"]
    df.loc["biomass", sector] += s_ued["Mach. Eq.: Steam processing"] / eff_biomass

    s_out = idees["out"][3:4]
    assert "Physical output" in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Physical output (index)"]
    )

    return df


def textiles_and_leather():
    sector = "Textiles and leather"

    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    df[sector] = 0.0

    s_fec = idees["fec"][3:58]
    s_ued = idees["ued"][3:58]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification
    # in new JRC data zero assume old data

    eff_elec = 73.7 / 146.6
    df.loc["elec", sector] += s_ued["Textiles: Drying"] / eff_elec

    df.loc["elec", sector] += s_fec["Textiles: Electric general machinery"]
    df.loc["elec", sector] += s_fec["Textiles: Finishing Electric"]

    # Steam processing is supplied with biomass
    eff_biomass = s_ued[15:26]["Biomass and waste"] / s_fec[15:26]["Biomass and waste"]
    df.loc["biomass", sector] += (
        s_ued["Textiles: Pretreatment with steam"] / eff_biomass
    )
    df.loc["biomass", sector] += (
        s_ued["Textiles: Wet processing with steam"] / eff_biomass
    )

    s_out = idees["out"][3:4]
    assert "Physical output" in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Physical output (index)"]
    )

    return df


def wood_and_wood_products():
    sector = "Wood and wood products"

    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    df[sector] = 0.0

    s_fec = idees["fec"][3:47]
    s_ued = idees["ued"][3:47]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification
    key = "Wood: Electric drying"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Wood: Drying"] / eff_elec

    df.loc["elec", sector] += s_fec["Wood: Electric mechanical processes"]
    df.loc["elec", sector] += s_fec["Wood: Finishing Electric"]

    # Steam processing is supplied with biomass
    eff_biomass = s_ued[15:25]["Biomass and waste"] / s_fec[15:25]["Biomass and waste"]
    df.loc["biomass", sector] += (
        s_ued["Wood: Specific processes with steam"] / eff_biomass
    )

    s_out = idees["out"][3:4]
    assert "Physical output" in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Physical output (index)"]
    )

    return df


def other_industrial_sectors():
    sector = "Other industrial sectors"
    idees = load_idees_data(sector)

    df = pd.DataFrame(index=index)

    df[sector] = 0.0

    s_fec = idees["fec"][3:68]
    s_ued = idees["ued"][3:68]
    assert s_fec.index[0] == sector
    assert s_ued.index[0] == sector

    sel = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    df.loc["elec", sector] += s_fec[sel].sum()

    df.loc["heat", sector] += s_fec["Low-enthalpy heat"]

    # Efficiency changes due to electrification
    key = "Other Industrial sectors: Electric processing"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += (
        s_ued["Other Industrial sectors: Process heating"] / eff_elec
    )

    key = "Other Industries: Electric drying"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += s_ued["Other Industrial sectors: Drying"] / eff_elec

    key = "Other Industries: Electric cooling"
    eff_elec = s_ued[key] / s_fec[key]
    df.loc["elec", sector] += (
        s_ued["Other Industrial sectors: Process Cooling"] / eff_elec
    )

    # Diesel motors are electrified
    key = "Other Industrial sectors: Diesel motors (incl. biofuels)"
    df.loc["elec", sector] += s_fec[key]
    key = "Other Industrial sectors: Electric machinery"
    df.loc["elec", sector] += s_fec[key]

    # Steam processing is supplied with biomass
    eff_biomass = s_ued[15:25]["Biomass and waste"] / s_fec[15:25]["Biomass and waste"]
    df.loc["biomass", sector] += (
        s_ued["Other Industrial sectors: Steam processing"] / eff_biomass
    )

    s_out = idees["out"][3:4]
    assert "Physical output" in str(s_out.index)

    # MWh/t material
    sources = ["elec", "biomass", "methane", "hydrogen", "heat", "naphtha"]
    df.loc[sources, sector] = (
        df.loc[sources, sector] * toe_to_MWh / s_out["Physical output (index)"]
    )

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industry_sector_ratios")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params.industry

    year = params["reference_year"]

    df = pd.concat(
        [
            iron_and_steel(),
            chemicals_industry(),
            nonmetalic_mineral_products(),
            pulp_paper_printing(),
            food_beverages_tobacco(),
            non_ferrous_metals(),
            transport_equipment(),
            machinery_equipment(),
            textiles_and_leather(),
            wood_and_wood_products(),
            other_industrial_sectors(),
        ],
        axis=1,
    )

    df.index.name = "MWh/tMaterial"
    df.to_csv(snakemake.output.industry_sector_ratios)
