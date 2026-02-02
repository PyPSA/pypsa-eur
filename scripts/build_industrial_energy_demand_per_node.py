# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build industrial energy demand per model region.

Description
-------
This rule aggregates the energy demand of the industrial sectors per model region.
For each bus, the following carriers are considered:
- electricity
- coal
- coke
- solid biomass
- methane
- hydrogen
- low-temperature heat
- naphtha
- ammonia
- process emission
- process emission from feedstock

which can later be used as values for the industry load.
"""

import logging

import numpy as np
import pandas as pd
import json
import holidays

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

# Industry category to FfE profile mapping dictionary
# 
# Note: FfE profiles missing "Industry (total)" (#0), "Chemical Industry" (#2) and "Non-ferrous Metals" (#3)
# 
# CHEMICAL SECTORS → "Paper, Pulp and Print"
# Based on Ganz et al. (2021): "Sowohl Papier- als auch Chemieindustrie weisen verhältnismäßig 
# konstante Verbräuche im Wochenverlauf auf, da ein Großteil der Prozesse in diesen Branchen im 
# Dauervollbetrieb fährt" [Both paper and chemical industries show relatively constant weekly 
# consumption due to continuous full-load operation of most processes].
# Applies to: Ammonia (Haber-Bosch), Chlorine (chlor-alkali), Methanol (catalytic synthesis)
# 
# NON-FERROUS METALS → "Iron & steel industry"
# Due to missing "Non-ferrous Metals" profile, all non-ferrous metal production mapped to 
# "Iron & steel industry" based on shared continuous/semi-continuous metal production characteristics.
# Includes: Aluminium (primary/secondary), Alumina, and other non-ferrous metals (Cu, Zn, Pb, etc.)
#
# PHARMACEUTICALS → "Food and Tobacco"
# Mapped due to similar batch production, quality control requirements, climate-controlled 
# environments, and reduced weekend production (not continuous 24/7 like paper/chemicals/metals).
#
# Reference: Ganz, K., Guminski, A., Kolb, M., & von Roon, S. (2021). Wie können europäische 
# Branchen-Lastgänge die Energiewende im Industriesektor unterstützen? ew - Magazin für die 
# Energiewirtschaft, 120(11), 42-45.

INDUSTRY_CATEGORY_TO_PROFILE = {
    "Electric arc": "Iron & steel industry",
    "DRI + Electric arc": "Iron & steel industry",
    "Integrated steelworks": "Iron & steel industry",
    "HVC": "Non-metallic Minerals",
    "HVC (mechanical recycling)": "Non-metallic Minerals",
    "HVC (chemical recycling)": "Non-metallic Minerals",
    "Ammonia": "Paper, Pulp and Print",  
    "Chlorine": "Paper, Pulp and Print",
    "Methanol": "Paper, Pulp and Print",
    "Other chemicals": "Paper, Pulp and Print",
    "Pharmaceutical products etc.": "Food and Tobacco",
    "Cement": "Non-metallic Minerals",
    "Ceramics & other NMM": "Non-metallic Minerals",
    "Glass production": "Non-metallic Minerals",
    "Pulp production": "Paper, Pulp and Print",
    "Paper production": "Paper, Pulp and Print",
    "Printing and media reproduction": "Paper, Pulp and Print",
    "Food, beverages and tobacco": "Food and Tobacco",
    "Alumina production": "Iron & steel industry",
    "Aluminium - primary production": "Iron & steel industry",
    "Aluminium - secondary production": "Iron & steel industry",
    "Other non-ferrous metals": "Non-metallic Minerals",
    "Transport equipment": "Transport Equipment",
    "Machinery equipment": "Machinery",
    "Textiles and leather": "Textile and Leather",
    "Wood and wood products": "Wood and Wood Products",
    "Other industrial sectors": "Non-specified (Industry)",
}


def load_ffe_load_profiles(json_file):
    """Load normalized industry load profiles from FfE JSON file."""

    # Load pre-downloaded JSON data
    with open(json_file, "r") as f:
        data = json.load(f)

    logger.info(f"Loaded FfE data: {data['title']}")

    # Map internal_id to profile names
    id_to_profile = {
        # 0: "Industry (total)", # missing
        1: "Iron & steel industry",
        # 2: "Chemical Industry", # missing
        # 3: "Non-ferrous Metals", # missing
        4: "Non-metallic Minerals",
        5: "Transport Equipment",
        6: "Machinery",
        7: "Mining and Quarrying",
        8: "Food and Tobacco",
        9: "Paper, Pulp and Print",
        10: "Wood and Wood Products",
        11: "Construction",
        12: "Textile and Leather",
        13: "Non-specified (Industry)",
    }

    # Create hourly timestamps for 2017
    timestamps = pd.date_range("2017-01-01", "2017-12-31 23:00:00", freq="h")

    # Parse the data into a DataFrame
    df = pd.json_normalize(data["data"])
    df = df.set_index(df["internal_id"].map(lambda x: x[0]))["values"]
    profiles_df = (
        pd.DataFrame(np.vstack(df), index=df.index, columns=timestamps)
        .rename(index=id_to_profile)
        .T
    )
    logger.info(f"Loaded profiles: {list(profiles_df.columns)}")

    return profiles_df


def create_nodal_electricity_profiles(
    nodal_df, nodal_sector_df, snapshots, ffe_profiles
):
    """Create hourly electricity demand profiles for each node."""

    # Initialize result DataFrame
    nodal_profiles = pd.DataFrame(index=snapshots, columns=nodal_df.index, dtype=float)

    # For each node, create weighted profile
    for node in nodal_df.index:

        # Get electricity demand by sector for this node (TWh/a)
        sector_demands = nodal_sector_df.loc["elec", node]

        # check if sum of sector demands matches total demand
        assert np.isclose(
            sector_demands.sum(), nodal_df.loc[node, "electricity"], rtol=1e-5
        ), f"Sum of sector demands does not match total demand for node {node}.\n"

        # Create weighted profile (8760 hours from reference year)
        node_profile = pd.Series(0.0, index=ffe_profiles.index)

        for sector, demand in sector_demands.items():
            profile_name = INDUSTRY_CATEGORY_TO_PROFILE[sector]
            if profile_name not in ffe_profiles.columns:
                raise ValueError(f"Profile '{profile_name}' for sector '{sector}' not in FfE data")
            node_profile += ffe_profiles[profile_name] * demand

        # Map profile to snapshots with correct day-of-week alignment
        hourly_profile = map_profile_to_snapshots(node_profile, snapshots, node_country=node[:2], tol=0.03) # Tolerance increased to 3% to account for holiday replacement and leap year
        # Save in MW
        nodal_profiles[node] = hourly_profile * 1e6 

    # Check that hourly profiles match annual demand
    assert np.allclose(
        nodal_profiles.sum() / 1e6, nodal_df["electricity"], rtol=1e-5
    ), (
        f"Hourly profiles do not match annual demand.\nMax difference: {(nodal_profiles.sum() / 1e6 - nodal_df['electricity']).abs().max():.6f} TWh"
    )
    logger.info("Hourly profiles verified to match annual demand")

    return nodal_profiles


def map_profile_to_snapshots(reference_profile, snapshots, node_country='DE', tol=0.03):
    """
    Map 2017 reference profile to target snapshots by weekday rotation.
    
    Process:
    1. Compute average day profiles (by dayofweek and hour)
    2. Replace German holidays in 2017 reference with weekday-hour averages
    3. Shift by day drift and remap to target year
    4. Handle leap year if needed
    5. Replace target country holidays with Sunday profile
    6. Optionally drop leap day if not in snapshots
    7. Scale to preserve energy
    
    Parameters
    ----------
    reference_profile : pd.Series
        Normalized load profile from 2017 with DatetimeIndex
    snapshots : pd.DatetimeIndex
        Target snapshots to map to
    node_country : str
        Two-letter country code for holiday replacement (default: 'DE')
    tol : float
        Tolerance for energy deviation (default: 0.03 = 3%)
    
    Returns
    -------
    pd.Series
        Mapped profile with snapshots as index, scaled to preserve total energy
    """
    
    original_energy = reference_profile.sum()
    target_year = snapshots[0].year
    source_year = 2017
    
    s = reference_profile.copy()
    
    # STEP 1: Compute average day profiles (dayofweek, hour)
    average_days = s.groupby([s.index.dayofweek, s.index.hour]).mean()
    
    # STEP 2: Normalize German holidays in source data
    de_holidays = holidays.country_holidays('DE', years=source_year)
    holiday_dates_de = pd.to_datetime(list(de_holidays.keys()))
    is_holiday_source = s.index.normalize().isin(holiday_dates_de)
    
    if is_holiday_source.any():
        s.loc[is_holiday_source] = (
            s.loc[is_holiday_source]
            .index.map(lambda i: average_days.loc[(i.dayofweek, i.hour)])
            .to_numpy()
        )
    
    # STEP 3: Add day drift
    day_drift = (
        pd.Timestamp(source_year, 1, 1).dayofweek
        - pd.Timestamp(target_year, 1, 1).dayofweek
    )
    
    ts = s.shift(day_drift, freq="D")
    ts.index = ts.index.map(lambda t: t.replace(year=target_year))
    ts = ts.sort_index()
    
    # STEP 4: Handle leap year
    is_target_leap = (target_year % 4 == 0 and target_year % 100 != 0) or (target_year % 400 == 0)
    
    if is_target_leap:
        # Shift dates after March back by 1 day
        ts.index = ts.index.where(ts.index.month < 3, ts.index - pd.Timedelta(days=1))
        
        # Add extra day at end of year using appropriate weekday profile
        extra_day_i = pd.date_range(f"{target_year}-12-31", periods=24, freq="h")
        extra_day_dow = extra_day_i[0].dayofweek
        
        extra_day = pd.Series(
            [average_days.loc[(extra_day_dow, h)] for h in range(24)],
            index=extra_day_i
        )
        ts = pd.concat([ts, extra_day]).sort_index()
    
    # STEP 5: Impose target country holidays (replace with Sunday profile)
    try:
        target_holidays = holidays.country_holidays(node_country, years=target_year)
        holiday_dates_target = pd.to_datetime(list(target_holidays.keys()))
        is_holiday_target = ts.index.normalize().isin(holiday_dates_target)
        
        if is_holiday_target.any():
            ts.loc[is_holiday_target] = (
                ts.loc[is_holiday_target]
                .index.map(lambda i: average_days.loc[(6, i.hour)])  # 6 = Sunday
                .to_numpy()
            )
    except NotImplementedError:
        logger.warning(f"Country '{node_country}' not available in holidays library - skipping holiday replacement")
    
    # STEP 6: Drop leap day if not in snapshots
    feb29_in_snapshots = any((snapshots.month == 2) & (snapshots.day == 29))
    
    if not feb29_in_snapshots and is_target_leap:
        ts = ts[~((ts.index.month == 2) & (ts.index.day == 29))]
    
    # STEP 7: Align with provided snapshots (in case of partial year or different resolution)
    # Reindex to match exact snapshots
    ts = ts.reindex(snapshots, method=None)
    
    # Handle any remaining NaNs (shouldn't happen with full year data)
    if ts.isna().any():
        missing_mask = ts.isna()
        ts.loc[missing_mask] = (
            snapshots[missing_mask]
            .map(lambda i: average_days.loc[(i.dayofweek, i.hour)])
            .to_numpy()
        )
    
    # STEP 8: Scale to preserve original energy
    scaling_factor = original_energy / ts.sum()
    
    assert abs(scaling_factor - 1.0) < tol, (
        f"Energy deviation after mapping: {(scaling_factor - 1.0) * 100:.2f}% "
        f"(tolerance: {tol*100:.1f}%)"
    )
    
    return ts * scaling_factor


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_energy_demand_per_node",
            clusters=50,
            planning_horizons=2030,
            run="temporal_industry_load",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # load industry sector ratios
    fn = snakemake.input.industry_sector_ratios
    sector_ratios = pd.read_csv(fn, header=[0, 1], index_col=0)

    # material demand per node and industry (Mton/a)
    fn = snakemake.input.industrial_production_per_node
    nodal_production = pd.read_csv(fn, index_col=0) / 1e3

    # energy demand today to get current electricity
    fn = snakemake.input.industrial_energy_demand_per_node_today
    nodal_today = pd.read_csv(fn, index_col=0)

    nodal_sector_ratios = pd.concat(
        {node: sector_ratios[node[:2]] for node in nodal_production.index}, axis=1
    )

    nodal_production_stacked = nodal_production.stack()
    nodal_production_stacked.index.names = [None, None]

    # final energy consumption per node and industry (TWh/a)
    nodal_df = (
        (nodal_sector_ratios.multiply(nodal_production_stacked))
        .T.groupby(level=0)
        .sum()
    )

    rename_sectors = {
        "elec": "electricity",
        "biomass": "solid biomass",
        "heat": "low-temperature heat",
    }
    nodal_df.rename(columns=rename_sectors, inplace=True)

    nodal_df["current electricity"] = nodal_today["electricity"]

    nodal_df.index.name = "TWh/a (MtCO2/a)"

    # Export annual demand
    fn = snakemake.output.industrial_energy_demand_per_node
    nodal_df.to_csv(fn)

    # Generate hourly electricity profiles
    # final energy consumption per node and industry and industry sector (profiles) (TWh/a)
    nodal_sector_df = nodal_sector_ratios.multiply(nodal_production_stacked)

    logger.info("Creating hourly industry electricity demand profiles...")
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    logger.info("Loading FfE industry load profiles...")
    ffe_profiles = load_ffe_load_profiles(snakemake.input.ffe_profiles)

    nodal_electricity_profiles = create_nodal_electricity_profiles(
        nodal_df, nodal_sector_df, snapshots, ffe_profiles
    )

    # Export hourly profiles
    fn_profiles = snakemake.output.industrial_electricity_demand_per_node_temporal
    nodal_electricity_profiles.to_csv(fn_profiles)
    logger.info(f"Hourly electricity profiles saved to {fn_profiles}")
