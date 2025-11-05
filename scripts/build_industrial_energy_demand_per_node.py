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

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def load_ffe_load_profiles(json_file):
    """Load normalized industry load profiles from FfE JSON file."""
    import json

    # Load pre-downloaded JSON data
    with open(json_file) as f:
        data = json.load(f)

    logger.info(f"Loaded FfE data: {data['title']}")

    # Map internal_id to profile names
    id_to_profile = {
        0: "Industry (total)",
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
    profiles_dict = {}
    for row in data["data"]:
        internal_id = row["internal_id"][0]
        values = row["values"]
        if internal_id in id_to_profile:
            profile_name = id_to_profile[internal_id]
            profiles_dict[profile_name] = values

    profiles_df = pd.DataFrame(profiles_dict, index=timestamps)
    logger.info(f"Loaded profiles: {list(profiles_df.columns)}")

    return profiles_df


def create_nodal_electricity_profiles(
    nodal_df, nodal_sector_df, snapshots, path_to_ffe_json
):
    """Create hourly electricity demand profiles for each node."""

    # Industry category to FfE profile mapping (updated to match correct profile names)
    industry_category_to_profile = {
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

    # Download FfE profiles
    logger.info("Downloading FfE industry load profiles...")
    ffe_profiles = load_ffe_load_profiles(path_to_ffe_json)

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
        total_weight = 0

        for sector, demand in sector_demands.items():
            if demand > 0 and sector in industry_category_to_profile:
                profile_name = industry_category_to_profile[sector]
                if profile_name in ffe_profiles.columns:
                    node_profile += ffe_profiles[profile_name] * demand
                    total_weight += demand
                else:
                    logger.warning(
                        f"Profile '{profile_name}' not found in FfE data for sector '{sector}'. Using flat demand"
                    )
                    flat_profile = pd.Series(
                        1.0 / len(ffe_profiles), index=ffe_profiles.index
                    )
                    node_profile += flat_profile * demand
                    total_weight += demand

        # Map profile to snapshots with correct day-of-week alignment
        hourly_profile = map_profile_to_snapshots(node_profile, snapshots)
        # Save in MW
        nodal_profiles[node] = hourly_profile * 1e6  # TWh/a to MW

    # Check that hourly profiles match annual demand
    assert np.allclose(
        nodal_profiles.sum() / 1e6, nodal_df["electricity"], rtol=1e-5
    ), (
        f"Hourly profiles do not match annual demand.\nMax difference: {(nodal_profiles.sum() / 1e6 - nodal_df['electricity']).abs().max():.6f} TWh"
    )
    logger.info("✓ Hourly profiles verified to match annual demand")

    return nodal_profiles


def map_profile_to_snapshots(reference_profile, snapshots):
    """Map a 2017 reference profile to the snapshot period, matching day-of-week and hour."""

    # Pre-compute lookup tables
    dayofweek_hour_lookup = reference_profile.groupby(
        [reference_profile.index.dayofweek, reference_profile.index.hour]
    ).mean()

    hour_lookup = reference_profile.groupby(reference_profile.index.hour).mean()
    overall_mean = reference_profile.mean()

    # Create keys for vectorized lookup
    keys = pd.MultiIndex.from_arrays([snapshots.dayofweek, snapshots.hour])

    # Primary lookup: match (dayofweek, hour)
    mapped_profile = dayofweek_hour_lookup.reindex(keys)
    mapped_profile.index = snapshots  # Restore original index

    # Fallback 1: For NaN values, try matching hour only
    mask = mapped_profile.isna()
    if mask.any():
        hour_values = hour_lookup.reindex(snapshots[mask].hour)
        mapped_profile[mask] = hour_values.values

        # Fallback 2: For remaining NaN, use overall mean
        mask = mapped_profile.isna()
        if mask.any():
            mapped_profile[mask] = overall_mean

    # Rescale to preserve total energy
    scaling_factor = reference_profile.sum() / mapped_profile.sum()

    # Check scaling is within expected range (day-of-week + leap year effects)
    assert abs(scaling_factor - 1.0) < 0.01, (
        f"Profile mapping scaling factor {scaling_factor:.4f} deviates by >1%"
    )

    mapped_profile = mapped_profile * scaling_factor

    return mapped_profile


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

    # import ratios
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

    nodal_electricity_profiles = create_nodal_electricity_profiles(
        nodal_df, nodal_sector_df, snapshots, snakemake.input.ffe_profiles
    )

    # Export hourly profiles
    fn_profiles = snakemake.output.industrial_electricity_demand_per_node_temporal
    nodal_electricity_profiles.to_csv(fn_profiles)
    logger.info(f"Hourly electricity profiles saved to {fn_profiles}")
