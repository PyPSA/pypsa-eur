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

# Industry category to FfE profile mapping (updated to match correct profile names)

INDUSTRY_CATEGORY_TO_PROFILE = {
    "Electric arc": "Iron & steel industry",
    "DRI + Electric arc": "Iron & steel industry",
    "Integrated steelworks": "Iron & steel industry",
    "HVC": "Non-specified (Industry)",
    "HVC (mechanical recycling)": "Non-specified (Industry)",
    "HVC (chemical recycling)": "Non-specified (Industry)",
    "Ammonia": "Non-specified (Industry)",
    "Chlorine": "Non-specified (Industry)",
    "Methanol": "Non-specified (Industry)",
    "Other chemicals": "Non-specified (Industry)",
    "Pharmaceutical products etc.": "Non-specified (Industry)",
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
    "Other non-ferrous metals": "Non-specified (Industry)",
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

       # Map internal_id to profile names (internal_id [2,3] missing from json file, reported back to FfE)
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
        nodal_df, nodal_sector_df, snapshots, ffe_profiles,       
):
    """Create hourly electricity demand profiles for each node."""

    # ---- Vectorised combination of sector profiles per node ----
    # Extract electricity sector demands for all nodes at once:
    # nodal_sector_df.loc["elec"] returns Series with MultiIndex (node, sector)
    # Unstack to get DataFrame with index=sectors, columns=nodes
    elec_demands_series = nodal_sector_df.loc["elec"]
    elec_demands = elec_demands_series.unstack(level=0)  # index=sectors, columns=nodes

    # Map each sector to the corresponding FfE profile name
    sector_to_profile = pd.Series(INDUSTRY_CATEGORY_TO_PROFILE)
    sector_to_profile = sector_to_profile.reindex(elec_demands.index)

    if sector_to_profile.isna().any():
        missing = sector_to_profile[sector_to_profile.isna()].index.tolist()
        raise KeyError(
            f"Missing FfE profile mapping for industrial sectors: {missing}"
        )

    # Aggregate sector demands to profile-level demands:
    #   index  -> FfE profile name
    #   columns -> nodes
    profile_demands = elec_demands.groupby(sector_to_profile).sum()

    # Ensure all required profiles exist in the FfE data
    missing_profiles = set(profile_demands.index) - set(ffe_profiles.columns)
    if missing_profiles:
        raise ValueError(
            f"Profiles {missing_profiles} from INDUSTRY_CATEGORY_TO_PROFILE not in FfE data"
        )

    # Compute 2017 reference profiles for all nodes in one matrix multiplication:
    #   ffe_profiles[profiles] : (hours x profiles)
    #   profile_demands        : (profiles x nodes)
    #   -> weighted_profiles   : (hours x nodes)
    weighted_profiles = ffe_profiles[profile_demands.index].dot(profile_demands)

    # ---- Per-node calendar mapping and checks (cannot easily be vectorised) ----
    # Create node-to-country mapping (instead of inferring from node name string)
    # Note: Currently derived from node name prefix, but structured as mapping
    # for easier replacement if proper country data becomes available
    node_to_country = pd.Series(nodal_df.index.str[:2], index=nodal_df.index)
    
    nodal_profiles = pd.DataFrame(index=snapshots, columns=nodal_df.index, dtype=float)

    for node in nodal_df.index:
        # Get electricity demand by sector for this node (TWh/a)
        sector_demands = elec_demands[node]

        # Check if sum of sector demands matches total demand
        assert np.isclose(
            sector_demands.sum(), nodal_df.loc[node, "electricity"], rtol=1e-5
        ), f"Sum of sector demands does not match total demand for node {node}. \n"

        # Node-specific 2017 reference profile (already sector-weighted)
        node_profile = weighted_profiles[node]

        # Map profile to snapshots with correct day-of-week alignment
        hourly_profile = map_profile_to_snapshots(
            node_profile,
            snapshots,
            node_country=node_to_country[node],
            tol=0.02,  # Tolerance increased to 2% for holiday replacement and date adjustment
        )

        # Save in MW
        nodal_profiles[node] = hourly_profile * 1e6
    
    # Check that hourly profiles match annual demand
    assert np.allclose(
        nodal_profiles.sum() / 1e6, nodal_df["electricity"], rtol=1e-5
    ), (
        f"Hourly profiles do not match annual demand. \nMax difference: {(nodal_profiles.sum() / 1e6 - nodal_df['electricity']).abs().max():.6f} TWh"
    )
    logger.info("Hourly profiles verified to match annual demand")

    return nodal_profiles


def create_nodal_electricity_profiles_per_profile(
    nodal_df, nodal_sector_df, snapshots, ffe_profiles, node_to_country
):
    """
    Create hourly electricity demand per FfE profile per node (MW).

    Returns a DataFrame with MultiIndex columns (node, profile) and index = snapshots.
    Used for industry DSR (Option B): one Store per profile per node.
    """
    elec_demands_series = nodal_sector_df.loc["elec"]
    elec_demands = elec_demands_series.unstack(level=0)

    sector_to_profile = pd.Series(INDUSTRY_CATEGORY_TO_PROFILE)
    sector_to_profile = sector_to_profile.reindex(elec_demands.index)
    if sector_to_profile.isna().any():
        missing = sector_to_profile[sector_to_profile.isna()].index.tolist()
        raise KeyError(
            f"Missing FfE profile mapping for industrial sectors: {missing}"
        )

    profile_demands = elec_demands.groupby(sector_to_profile).sum()
    missing_profiles = set(profile_demands.index) - set(ffe_profiles.columns)
    if missing_profiles:
        raise ValueError(
            f"Profiles {missing_profiles} from INDUSTRY_CATEGORY_TO_PROFILE not in FfE data"
        )

    profiles = profile_demands.index.tolist()
    nodes = nodal_df.index.tolist()

    # Build (snapshots x (node, profile)) in MW
    data = {}
    for node in nodes:
        for profile in profiles:
            annual_twh = profile_demands.loc[profile, node]
            if annual_twh <= 0:
                # Zero demand: skip map_profile_to_snapshots (avoids 0/0 -> NaN in scaling)
                data[(node, profile)] = pd.Series(0.0, index=snapshots)
            else:
                ref_profile = ffe_profiles[profile] * annual_twh
                mapped = map_profile_to_snapshots(
                    ref_profile,
                    snapshots,
                    node_country=node_to_country[node],
                    tol=0.02,
                )
                data[(node, profile)] = mapped * 1e6
    per_profile = pd.DataFrame(data, index=snapshots)
    per_profile.columns.names = ["node", "profile"]
    return per_profile


def map_profile_to_snapshots(reference_profile, snapshots, node_country='DE', tol=0.02):
    """
    Map 2017 German reference profile to target snapshots by weekday rotation.

    Process:
    1. Remove German holidays from 2017 reference (replace with weekday-hour avg)
    2. Rotate and map to target year
    3. Replace target country's holidays with Sunday-hour averages (if available)
    4. Drop leap day if needed (handled at the end)
    5. Scale once at the end
    """
    SOURCE_YEAR = 2017
    TARGET_YEAR = snapshots[0].year
    
    # Work with a copy of the reference profile
    s = reference_profile.copy()
    original_energy = s.sum()
    
    # STEP 1: Compute average day profiles
    average_days = s.groupby([s.index.dayofweek, s.index.hour]).mean()
    
    # STEP 2: Normalize German holidays in source year
    de_holidays = holidays.country_holidays('DE', years=SOURCE_YEAR)
    holiday_dates_de = pd.to_datetime(list(de_holidays.keys()))
    is_holiday = s.index.normalize().isin(holiday_dates_de)
    
    if is_holiday.any():
        s.loc[is_holiday] = (
            s.loc[is_holiday]
            .index.map(lambda i: average_days.loc[(i.dayofweek, i.hour)])
            .to_numpy()
        )
    
    # STEP 3: Add day drift (weekday alignment)
    day_drift = (
        pd.Timestamp(SOURCE_YEAR, 1, 1).dayofweek
        - pd.Timestamp(TARGET_YEAR, 1, 1).dayofweek
    )
    
    ts = s.shift(day_drift, freq="D")
    ts.index = ts.index.map(lambda t: t.replace(year=TARGET_YEAR))
    ts.sort_index(inplace=True)
    
    # STEP 4: Handle leap year
    if pd.Timestamp(TARGET_YEAR, 1, 1).is_leap_year:
        # Shift dates from March onward back by 1 day to account for leap year
        mask = ts.index.month >= 3
        ts.index = pd.DatetimeIndex(
            np.where(mask, ts.index - pd.Timedelta(days=1), ts.index)
        )
        
        # Add extra day (Dec 31) for leap year
        extra_day_i = pd.date_range(f"{TARGET_YEAR}-12-31", periods=24, freq="h")
        # Get average profile for the weekday of Dec 31
        weekday_dec31 = extra_day_i[0].dayofweek
        extra_day_values = [average_days.loc[(weekday_dec31, hour)] for hour in range(24)]
        extra_day = pd.Series(extra_day_values, index=extra_day_i)
        ts = pd.concat([ts, extra_day]).sort_index()
    
    # STEP 5: Impose target country holidays
    try:
        target_holidays = holidays.country_holidays(node_country, years=TARGET_YEAR)
        target_holiday_dates = pd.to_datetime(list(target_holidays.keys()))
        is_holiday = ts.index.normalize().isin(target_holiday_dates)
        
        if is_holiday.any():
            ts.loc[is_holiday] = (
                ts.loc[is_holiday].index.map(lambda i: average_days.loc[(6, i.hour)]).to_numpy()
            )
    except NotImplementedError:
        logger.warning(f"Country '{node_country}' not available in holidays library - skipping holiday replacement")
    
    # STEP 6: Drop leap day if needed (snapshots already have it dropped)
    drop_leap_day = not any((snapshots.month == 2) & (snapshots.day == 29))
    if drop_leap_day and ts.index.is_leap_year.any():
        ts = ts[~((ts.index.month == 2) & (ts.index.day == 29))]
    
    # STEP 7: Map to snapshots (reindex to match exact snapshot timestamps)
    mapped = ts.reindex(snapshots)
    
    # Fill any missing values (shouldn't happen, but safety fallback)
    if mapped.isna().any():
        missing = mapped.isna()
        mapped.loc[missing] = (
            mapped.loc[missing].index.map(lambda i: average_days.loc[(i.dayofweek, i.hour)]).to_numpy()
        )
    
    # STEP 8: Scale to preserve energy (with tolerance check)
    mapped_sum = mapped.sum()
    if original_energy == 0 or np.isnan(original_energy) or mapped_sum == 0 or np.isnan(mapped_sum):
        # Avoid 0/0 or nan/... leading to NaN scaling_factor and assertion failure
        return pd.Series(0.0, index=snapshots)
    scaling_factor = original_energy / mapped_sum
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
            run="temporal_industry_load"
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

    # Create node-to-country mapping (instead of inferring from node name string)
    node_to_country = pd.Series(nodal_production.index.str[:2], index=nodal_production.index)
    
    nodal_sector_ratios = pd.concat(
        {node: sector_ratios[node_to_country[node]] for node in nodal_production.index}, axis=1
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

    #Generate hourly electricity profiles

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

    # Export hourly profiles (total per node)
    fn_profiles = snakemake.output.industrial_electricity_demand_per_node_temporal
    nodal_electricity_profiles.to_csv(fn_profiles, float_format="%.2f")
    logger.info(f"Hourly electricity profiles saved to {fn_profiles}")

    # Export hourly profiles per FfE profile per node (for industry DSR Option B)
    per_profile = create_nodal_electricity_profiles_per_profile(
        nodal_df, nodal_sector_df, snapshots, ffe_profiles, node_to_country
    )
    fn_per_profile = snakemake.output.industrial_electricity_demand_per_profile_temporal
    per_profile_flat = per_profile.copy()
    per_profile_flat.columns = [f"{n}|{p}" for n, p in per_profile.columns]
    per_profile_flat.to_csv(fn_per_profile, float_format="%.2f")
    logger.info(f"Hourly electricity profiles per profile saved to {fn_per_profile}")

