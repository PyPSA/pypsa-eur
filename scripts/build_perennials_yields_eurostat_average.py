# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve crop harvest data from the Eurostat API (dataset ``apro_cpshr``) at
NUTS2 and NUTS0 resolution, and compute production-weighted average yields
(t/ha, then converted to MWh/ha for 1G biofuel crops) per NUTS2 region for
1st-generation biofuel feedstocks and perennial grasses.

"Eurostat average" in the name refers to this production-weighted averaging
across crop codes and years - not the geometric area-weighted overlay used
downstream in ``build_perennials_yields.py`` to reproject these NUTS2-level
yields onto clustered network regions; do not confuse the two.

Missing NUTS2 coverage is filled via a three-tier fallback (NUTS2 data ->
NUTS0 country-level data -> spatial neighbor mean or nearest valid NUTS2 for
islands - see ``harmonize_to_nuts2021()``), and results are harmonized to
the NUTS2021 region definitions used by PyPSA-Eur.

Outputs a single CSV with columns for each crop class (cereals, sugar beet,
rapeseed, perennials) indexed by NUTS2 region.

Biofuel conversion efficiencies (t_biofuel / t_feedstock) are read from the
``efficiency`` parameter of the ``ethanol from wheat``, ``ethanol from sugar
beet``, and ``biodiesel from rapeseed`` technologies in the technology-data
cost assumptions, sourced from:

    Banja et al. (2013), "Biofuels in the European Union - A general overview",
    JRC Technical Report, doi:10.2760/69179, Tables 93, 133, 155, 159.
"""

import logging
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from scripts._helpers import load_costs

logger = logging.getLogger(__name__)


def harmonize_to_nuts2021(df, keep_col, nuts2021_n2):
    """
    Df : DataFrame indexed by ['geo', 'TIME_PERIOD', 'mapping']
         contains both NUTS2 and NUTS0 rows
    keep_col : column to harmonize (e.g. 'weighted_YL_(t/ha)')
    nuts2021_n2 : GeoDataFrame with index=NUTS2_ID and geometry
    """

    # NUTS2 target regions
    target = nuts2021_n2.index.sort_values()

    # Split df
    df = df.copy()

    # Identify NUTS2 and NUTS0 by index length
    df["is_nuts2"] = df.index.str.len() == 4

    # Split into NUTS2 and NUTS0
    df_nuts2 = df[df["is_nuts2"]]
    df_nuts0 = df[~df["is_nuts2"]]

    # Prepare working layer for only NUTS2
    df_work = df_nuts2[[keep_col]].reindex(target)

    # Country code extraction (first 2 chars)
    df_work["country"] = df_work.index.str[:2]

    # Fallback 2 — Use NUTS-0 values
    df_work = df_work.join(
        df_nuts0[[keep_col]].rename(columns={keep_col: "fallback"}), on="country"
    )
    df_work[keep_col] = df_work[keep_col].fillna(df_work["fallback"])
    df_work.drop(columns=["fallback"], inplace=True)

    # Join geometry
    nuts_proj = nuts2021_n2.to_crs(epsg=3035)
    gdf = nuts_proj.join(df_work)[[keep_col, "country", "geometry"]]

    # Fallback 3 — Spatial neighbors mean or nearest region if island
    mask_missing = gdf[keep_col].isna() | (gdf[keep_col] <= 0)
    if mask_missing.any():
        logger.warning("Spatial fallback required for %d regions", mask_missing.sum())

        # Pre-calc distance matrix only once
        valid = gdf[gdf[keep_col].notna() & (gdf[keep_col] > 0)]

        for idx in gdf[mask_missing].index:
            region = gdf.loc[idx, "geometry"]

            # Touching neighbors
            neigh_idxs = gdf[gdf.geometry.touches(region)].index.tolist()
            neigh_vals = gdf.loc[neigh_idxs, keep_col].dropna()
            neigh_vals = neigh_vals[neigh_vals > 0]

            if len(neigh_vals) > 0:
                gdf.at[idx, keep_col] = neigh_vals.mean()
            else:
                # Island fallback: nearest valid NUTS2 region
                nearest_idx = valid.distance(region).idxmin()
                gdf.at[idx, keep_col] = valid.at[nearest_idx, keep_col]
                logger.debug("Spatial fallback for %s: nearest = %s", idx, nearest_idx)

    # Final output tidy
    result = gdf[[keep_col]]
    result.index.name = "NUTS2"
    result.sort_index()
    return result


def calculate_yields(
    filepath_nuts2, filepath_nuts0, crops_sel, crops_mapping, biofuel_yields
):
    # filter columns - keep only relevant
    df_crops_raw_nuts2 = pd.read_csv(filepath_nuts2)
    df_crops_raw_nuts2["TIME_PERIOD"] = df_crops_raw_nuts2["TIME_PERIOD"].astype(int)

    df_crops_raw_nuts0 = pd.read_csv(filepath_nuts0)
    df_crops_raw_nuts0["TIME_PERIOD"] = df_crops_raw_nuts0["TIME_PERIOD"].astype(int)

    df_crops_raw = pd.concat(
        [df_crops_raw_nuts0, df_crops_raw_nuts2], ignore_index=True
    )

    # drop empty and irrelevant columns
    columns_to_drop = [
        "Observation value",
        "OBS_FLAG",
        "Observation status (Flag) V2 structure",
        "CONF_STATUS",
        "Confidentiality status (flag)",
        "Time",
        "STRUCTURE_ID",
        "STRUCTURE",
        "STRUCTURE_NAME",
        "Geopolitical entity (reporting)",
        "Time frequency",
    ]

    # useful columns:
    # Index(['freq', 'crops', 'Crops', 'strucpro', 'Structure of production', 'geo','TIME_PERIOD', 'OBS_VALUE'],
    # freq : 'A' meaning annual
    # 'crops' : code of the crops e.g. 'G0000', 'G1000', 'G2000', 'G2100', 'G2900'
    # 'Crops' : name of the corp e.g. Sugar beet (excluding seed
    # 'strucpro' : type of data ['AR', 'MA', 'PR_HU_EU'], where AR is  cultivated area, MA is  and PR_HU_EU is production at standard EU humidity
    # 'OBS_VALUE' : numerical value

    df_crops = df_crops_raw.drop(columns=columns_to_drop, errors="ignore")
    df_crops["OBS_VALUE"] = df_crops["OBS_VALUE"].fillna(0)

    # Step 1: Filter to relevant rows for 2023 and strucpro of interest
    df_sub = df_crops[
        (df_crops["strucpro"].isin(["AR", "PR_HU_EU"]))
        & (df_crops["crops"].isin(crops_sel))
    ][["crops", "geo", "TIME_PERIOD", "strucpro", "OBS_VALUE"]]

    # Pivot so AR and PR_HU_EU are columns for each (crop, geo, year)
    df_pivot = (
        df_sub.pivot_table(
            index=["crops", "geo", "TIME_PERIOD"],
            columns="strucpro",
            values="OBS_VALUE",
        )
        .dropna(subset=["AR", "PR_HU_EU"])
        .reset_index()
    )

    # Compute yield per year (PR_HU_EU / AR)
    df_pivot["YL_(t/ha)"] = np.divide(
        df_pivot["PR_HU_EU"],
        df_pivot["AR"],
        out=np.zeros_like(df_pivot["PR_HU_EU"], dtype=float),
        where=df_pivot["AR"] != 0,
    )

    # Average data across years
    df_avg_yield = df_pivot.groupby(["crops", "geo"], as_index=False)[
        ["AR", "PR_HU_EU", "YL_(t/ha)"]
    ].mean()

    min_year = df_pivot["TIME_PERIOD"].min()
    max_year = df_pivot["TIME_PERIOD"].max()
    df_avg_yield["TIME_PERIOD"] = f"{min_year}-{max_year}"

    # map crops to categories unsustainable biofuels in
    rev_map = {
        code: key
        for key, val in crops_mapping.items()
        for code in (val if isinstance(val, list) else [val])
    }
    df_avg_yield["mapping"] = df_avg_yield["crops"].map(rev_map)

    # calculate weighted production per crop within mapping classes
    df_avg_yield["PR_share"] = df_avg_yield["PR_HU_EU"] / df_avg_yield.groupby(
        ["geo", "TIME_PERIOD", "mapping"]
    )["PR_HU_EU"].transform("sum")
    df_avg_yield["PR_share"] = df_avg_yield["PR_share"].fillna(0)

    # calculated average weighted yield
    df_avg_yield["weighted_YL_(t/ha)"] = (
        df_avg_yield["YL_(t/ha)"] * df_avg_yield["PR_share"]
    )

    # sanity check for very low yields due to small productions
    thresholds = {
        "MINBIOCRP11": 2.0,  # cereals
        "MINBIOCRP21": 50.0,  # sugar beet
        "MINBIORPS1": 1.5,  # rapeseed
        "PERENNIALS": 5.0,  # perennial grasses
    }

    # Apply crop-specific minimum threshold
    df_avg_yield["weighted_YL_(t/ha)"] = df_avg_yield.apply(
        lambda row: (
            row["weighted_YL_(t/ha)"]
            if row["weighted_YL_(t/ha)"] >= thresholds.get(row["mapping"], 0)
            else 0
        ),
        axis=1,
    )

    # weighted yields from current production : applies to unsustainable biofuels
    weighted_yields = df_avg_yield.groupby(["geo", "TIME_PERIOD", "mapping"])[
        "weighted_YL_(t/ha)"
    ].sum()

    unsustainable_biofuels_yields = pd.DataFrame(weighted_yields)

    # unsustainable biomass yield units from t/ha to MWh/ha
    unsustainable_biofuels_yields = unsustainable_biofuels_yields[
        unsustainable_biofuels_yields.index.get_level_values("mapping") != "PERENNIALS"
    ]

    unsustainable_biofuels_yields["energy_yields_(MWh/ha)"] = (
        unsustainable_biofuels_yields["weighted_YL_(t/ha)"]
        * unsustainable_biofuels_yields.index.get_level_values("mapping").map(
            biofuel_yields
        )
    )

    # yields of perennials per hectare in ton/ha
    # standard humidity for perennials = 0.65 (tH2O/t_fresh) -> note production is for fresh until 2025
    std_moist_perennials = 0.65

    perennial_yields = pd.DataFrame(weighted_yields)
    perennial_yields = perennial_yields[
        perennial_yields.index.get_level_values("mapping") == "PERENNIALS"
    ] * (1 - std_moist_perennials)

    # max yields from current production : applies to perennials for green biorefining
    max_yields = df_avg_yield.groupby(["geo", "TIME_PERIOD", "mapping"])[
        "YL_(t/ha)"
    ].max()

    perennial_yields_max = pd.DataFrame(max_yields)
    perennial_yields_max = perennial_yields_max[
        perennial_yields_max.index.get_level_values("mapping") == "PERENNIALS"
    ] * (1 - std_moist_perennials)

    return unsustainable_biofuels_yields, perennial_yields, perennial_yields_max


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_perennials_yields_eurostat_average")

    from scripts._helpers import configure_logging, set_scenario_config

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    CROPS_CSV_NUTS2 = Path(snakemake.input["crops_nuts2"])
    CROPS_CSV_NUTS0 = Path(snakemake.input["crops_nuts0"])
    NUTS2_2021_GEOJSON = Path(snakemake.input["nuts2021"])
    OUT_CSV_YIELDS_ALL = Path(snakemake.output["yields_all"])

    CROPS_CSV_NUTS2.parent.mkdir(parents=True, exist_ok=True)
    OUT_CSV_YIELDS_ALL.parent.mkdir(parents=True, exist_ok=True)

    # G0000: total green plants (kept – often the only G-code with NUTS2 coverage, e.g. DK)
    # G1000: temporary grasses and grazings
    # G2000: aggregate legumes (G2100 + G2900) – included for NUTS0 fallback coverage;
    #         does NOT inflate MAX because it is always ≤ max(G2100, G2900)
    # G2100: lucerne/alfalfa; G2900: clover and other leguminous plants
    perennial_codes = ["G0000", "G1000", "G2000", "G2100", "G2900"]

    crops_mapping = dict(
        MINBIOCRP11=["C0000", "C1000", "C1210", "C1300", "C1310", "C1320"],
        MINBIOCRP21="R2000",
        MINBIORPS1=["I1110", "I1120", "I1130", "I1110-1130", "I0000"],
        PERENNIALS=perennial_codes,
    )
    costs = load_costs(snakemake.input.costs)

    # conv = snakemake.params.biofuel_conversion
    # LHV values are fixed physical constants, not parameters: JRC Technical Report doi:10.2760/69179
    LHV_fuels = {
        "ethanol": 7.447,
        "biodiesel": 10.194,
    }  # MWh/t (26.81 MJ/kg, 36.7 MJ/kg)

    biofuel_yields = {
        "MINBIOCRP11": costs.at["ethanol from wheat", "efficiency"]
        * LHV_fuels["ethanol"],
        "MINBIOCRP21": costs.at["ethanol from sugar beet", "efficiency"]
        * LHV_fuels["ethanol"],
        "MINBIORPS1": costs.at["biodiesel from rapeseed", "efficiency"]
        * LHV_fuels["biodiesel"],
    }

    other_crops_codes = [
        item
        for v in crops_mapping.values()
        for item in (v if isinstance(v, list) else [v])
    ]
    crops_sel = perennial_codes + other_crops_codes

    logger.info("Computing crop yields...")
    unsustainable_biofuels_yields, perennial_yields, perennial_yields_max = (
        calculate_yields(
            filepath_nuts0=CROPS_CSV_NUTS0,
            filepath_nuts2=CROPS_CSV_NUTS2,
            crops_sel=crops_sel,
            crops_mapping=crops_mapping,
            biofuel_yields=biofuel_yields,
        )
    )

    yield_MINBIOCRP11 = unsustainable_biofuels_yields[
        unsustainable_biofuels_yields.index.get_level_values("mapping") == "MINBIOCRP11"
    ]
    yield_MINBIOCRP21 = unsustainable_biofuels_yields[
        unsustainable_biofuels_yields.index.get_level_values("mapping") == "MINBIOCRP21"
    ]
    yield_MINBIORPS1 = unsustainable_biofuels_yields[
        unsustainable_biofuels_yields.index.get_level_values("mapping") == "MINBIORPS1"
    ]

    yield_MINBIOCRP11 = yield_MINBIOCRP11.droplevel(["TIME_PERIOD", "mapping"])
    yield_MINBIOCRP11.index.name = "NUTS2"
    yield_MINBIOCRP21 = yield_MINBIOCRP21.droplevel(["TIME_PERIOD", "mapping"])
    yield_MINBIOCRP21.index.name = "NUTS2"
    yield_MINBIORPS1 = yield_MINBIORPS1.droplevel(["TIME_PERIOD", "mapping"])
    yield_MINBIORPS1.index.name = "NUTS2"
    perennial_yields = perennial_yields.droplevel(["TIME_PERIOD", "mapping"])
    perennial_yields.index.name = "NUTS2"
    perennial_yields_max = perennial_yields_max.droplevel(["TIME_PERIOD", "mapping"])
    perennial_yields_max.index.name = "NUTS2"

    logger.info("Harmonizing to NUTS2021 regions...")
    nuts2021_n2 = (
        gpd.read_file(NUTS2_2021_GEOJSON)
        .loc[:, ["NUTS_ID", "NUTS_NAME", "CNTR_CODE", "geometry"]]
        .set_index("NUTS_ID")
    )

    yield_MINBIOCRP11_full = harmonize_to_nuts2021(
        yield_MINBIOCRP11, "energy_yields_(MWh/ha)", nuts2021_n2
    )
    yield_MINBIOCRP21_full = harmonize_to_nuts2021(
        yield_MINBIOCRP21, "energy_yields_(MWh/ha)", nuts2021_n2
    )
    yield_MINBIORPS1_full = harmonize_to_nuts2021(
        yield_MINBIORPS1, "energy_yields_(MWh/ha)", nuts2021_n2
    )
    yields_perennials_max_full = harmonize_to_nuts2021(
        perennial_yields_max, "YL_(t/ha)", nuts2021_n2
    )
    yields_perennials_full = harmonize_to_nuts2021(
        perennial_yields, "weighted_YL_(t/ha)", nuts2021_n2
    )

    df_yields_all = pd.concat(
        {
            "MINBIOCRP11": yield_MINBIOCRP11_full["energy_yields_(MWh/ha)"],
            "MINBIOCRP21": yield_MINBIOCRP21_full["energy_yields_(MWh/ha)"],
            "MINBIORPS1": yield_MINBIORPS1_full["energy_yields_(MWh/ha)"],
            # Max yield across G-codes: models the best-available perennial crop
            # choice in each region when substituting 1G biofuel crops.
            "PERENNIALS_MAX": yields_perennials_max_full["YL_(t/ha)"],
        },
        axis=1,
    )
    df_yields_all.columns = [
        "Bioethanol barley, wheat, grain maize, oats, other cereals and rye",
        "Sugar from sugar beet",
        "Rape seed",
        "perennials",
    ]
    df_yields_all = df_yields_all.sort_index()

    logger.info("Saving output CSV files...")
    df_yields_all.to_csv(OUT_CSV_YIELDS_ALL, index=True)

    logger.info("Done.")
