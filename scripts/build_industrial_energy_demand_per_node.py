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
INDUSTRY_CATEGORY_TO_PROFILE = {
    "Electric arc": "Iron & steel industry",
    "DRI + Electric arc": "Iron & steel industry",
    "Integrated steelworks": "Iron & steel industry",
    "HVC": "Non-metallic Minerals",
    "HVC (mechanical recycling)": "Non-metallic Minerals",
    "HVC (chemical recycling)": "Non-metallic Minerals",
    "Ammonia": "Paper, Pulp and Print", # approx. profile
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
    nodal_df, nodal_sector_df, snapshots, ffe_profiles, industry_category_to_profile=INDUSTRY_CATEGORY_TO_PROFILE
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
            profile_name = industry_category_to_profile[sector]
            if profile_name not in ffe_profiles.columns:
                raise ValueError(f"Profile '{profile_name}' for sector '{sector}' not in FfE data")
            node_profile += ffe_profiles[profile_name] * demand

        # Map profile to snapshots with correct day-of-week alignment
        hourly_profile = map_profile_to_snapshots(node_profile, snapshots, node_country=node[:2], tol=0.02) # Tolerance increased to 2% to account for holiday replacement and date adjustments
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


def map_profile_to_snapshots(reference_profile, snapshots, node_country='DE', tol=0.02):
    """
    Map 2017 German reference profile to target snapshots by weekday rotation.
    
    Process:
    1. Remove German holidays from 2017 reference (replace with weekday-hour avg)
    2. Rotate and map to target year
    3. Replace target country's holidays with Sunday-hour averages (if available)
    4. Scale once at the end
    """

    original_energy = reference_profile.sum()
    
    target_year = snapshots[0].year
    is_target_leap = (target_year % 4 == 0 and target_year % 100 != 0) or (target_year % 400 == 0)
    feb29_in_snapshots = any((snapshots.month == 2) & (snapshots.day == 29))
    
    # Detect if Feb 29 was dropped from leap year
    feb29_was_dropped = False
    if is_target_leap and not feb29_in_snapshots:
        feb28 = snapshots[(snapshots.month == 2) & (snapshots.day == 28)]
        mar1 = snapshots[(snapshots.month == 3) & (snapshots.day == 1)]
        if len(feb28) > 0 and len(mar1) > 0:
            feb29_was_dropped = (feb28[0].dayofweek + 1) % 7 != mar1[0].dayofweek
    
    # STEP 1: Replace GERMAN holidays in 2017 reference with weekday-hour averages
    reference_profile = reference_profile.copy()
    de_holidays = holidays.country_holidays('DE', years=2017)
    holiday_dates_de = pd.to_datetime(list(de_holidays.keys()))
    holiday_mask_de = reference_profile.index.normalize().isin(holiday_dates_de)
    
    if holiday_mask_de.any():
        non_holiday_mask_de = ~holiday_mask_de
        weekday_hour_avg = reference_profile[non_holiday_mask_de].groupby(
            [reference_profile[non_holiday_mask_de].index.dayofweek,
            reference_profile[non_holiday_mask_de].index.hour]
        ).mean()
        
        holiday_keys = pd.MultiIndex.from_arrays([
            reference_profile[holiday_mask_de].index.dayofweek,
            reference_profile[holiday_mask_de].index.hour
        ])
        reference_profile.loc[holiday_mask_de] = weekday_hour_avg.reindex(holiday_keys).values
    
    # STEP 2: Rotate reference by weekday offset
    ref_jan1_dow = pd.Timestamp('2017-01-01').dayofweek 
    target_jan1_dow = pd.Timestamp(f'{target_year}-01-01').dayofweek
    offset_days = (target_jan1_dow - ref_jan1_dow) % 7
    offset_hours = offset_days * 24
    
    rotated_values = np.roll(reference_profile.values, -offset_hours)
    rotated = pd.Series(rotated_values, index=reference_profile.index)
    
    # Create (month, day, hour) lookup
    ref_lookup = pd.Series(
        rotated.values,
        index=pd.MultiIndex.from_arrays([
            rotated.index.month,
            rotated.index.day,
            rotated.index.hour
        ])
    )
    
    # Build snapshot keys - adjust for dropped Feb 29 if needed
    snap_dates = snapshots.to_series()
    if feb29_was_dropped:
        snap_dates = snap_dates.copy()
        snap_dates[snap_dates.dt.month >= 3] -= pd.Timedelta(days=1)
    
    snap_keys = pd.MultiIndex.from_arrays([
        snap_dates.dt.month,
        snap_dates.dt.day,
        snap_dates.dt.hour
    ])
    
    # Vectorized mapping
    mapped = pd.Series(ref_lookup.reindex(snap_keys).values, index=snapshots)
    
    # Fill NaNs (Feb 29 if present) with weekday-hour averages
    if mapped.isna().any():
        fallback = reference_profile.groupby(
            [reference_profile.index.dayofweek, reference_profile.index.hour]
        ).mean()
        
        missing = mapped.isna()
        mapped.loc[missing] = fallback.reindex(
            pd.MultiIndex.from_arrays([
                snapshots[missing].dayofweek,
                snapshots[missing].hour
            ])
        ).values
    
    # STEP 3: Replace TARGET COUNTRY holidays with Sunday-hour averages
    try:
        target_holidays = holidays.country_holidays(node_country, years=target_year)
        target_holiday_dates = pd.to_datetime(list(target_holidays.keys()))
        target_holiday_mask = snapshots.normalize().isin(target_holiday_dates)
        
        if target_holiday_mask.any():
            non_target_holiday = ~target_holiday_mask
            sunday_mask = snapshots[non_target_holiday].dayofweek == 6
            
            if sunday_mask.any():
                sunday_hour_avg = mapped[non_target_holiday][sunday_mask].groupby(
                    snapshots[non_target_holiday][sunday_mask].hour
                ).mean()
                
                target_holiday_hours = snapshots[target_holiday_mask].hour
                mapped.loc[target_holiday_mask] = sunday_hour_avg.reindex(target_holiday_hours).values
    
    except NotImplementedError:
        logger.warning(f"Country '{node_country}' not available in holidays library - skipping holiday replacement")
    
    # STEP 4: Energy preservation
    scaling_factor = original_energy / mapped.sum()
    assert abs(scaling_factor - 1.0) < tol, (
        f"Energy deviation after mapping: {(scaling_factor - 1.0) * 100:.2f}%"
    )
    
    return mapped * scaling_factor


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
