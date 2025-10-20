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


import requests
import numpy as np
from datetime import datetime
import json

import logging

import pandas as pd


from scripts._helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def download_ffe_load_profiles():
    """Download normalized industry load profiles from FfE."""
    url = "https://api.opendata.ffe.de"
    params = {"id_opendata": 59}

    # Check API health
    response = requests.get(url + "/health")
    if response.status_code != 200:
        raise ConnectionError(f"API not available. Status: {response.status_code}")

    # Fetch data
    response = requests.get(url + "/opendata", params=params)
    if response.status_code != 200:
        raise ConnectionError(f"Request failed. Status: {response.status_code}")

    data = response.json()
    logger.info(f"Downloaded FfE data: {data['title']}")

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
        # Data is already parsed as lists, not JSON strings
        internal_id = row["internal_id"][0]  # Extract first element from list [5]
        values = row["values"]  # Already a list

        if internal_id in id_to_profile:
            profile_name = id_to_profile[internal_id]
            profiles_dict[profile_name] = values

    profiles_df = pd.DataFrame(profiles_dict, index=timestamps)

    logger.info(f"Loaded profiles: {list(profiles_df.columns)}")
    return profiles_df


def create_nodal_electricity_profiles(
    nodal_df, nodal_production, sector_ratios, snapshots
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
    ffe_profiles = download_ffe_load_profiles()

    # Initialize result DataFrame
    nodal_profiles = pd.DataFrame(index=snapshots, columns=nodal_df.index, dtype=float)

    # For each node, create weighted profile
    for node in nodal_production.index:
        country = node[:2]

        # Get electricity demand by sector for this node (MWh/tMaterial)
        elec_ratios = sector_ratios[country].loc["elec"]

        # Get production by sector for this node (Mton/a)
        production = nodal_production.loc[node]

        # Calculate electricity demand by sector (TWh/a)
        sector_demands = elec_ratios * production / 1e6  # Convert MWh to TWh

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
                        f"Profile '{profile_name}' not found in FfE data for sector '{sector}'"
                    )

        if total_weight == 0:
            logger.warning(f"No valid profile for node {node}, using flat profile")
            node_profile[:] = 1.0 / len(node_profile)
        else:
            # Normalize profile so it sums to 1
            node_profile = node_profile / node_profile.sum()

        # Map profile to snapshots with correct day-of-week alignment
        hourly_profile = map_profile_to_snapshots(node_profile, snapshots)

        # Scale to annual demand (convert TWh/a to MW hourly values)
        annual_demand_twh = nodal_df.loc[node, "electricity"]
        annual_demand_mwh = annual_demand_twh * 1e6  # TWh to MWh

        # Scale normalized profile to match annual demand
        nodal_profiles[node] = hourly_profile * annual_demand_mwh

    return nodal_profiles


def map_profile_to_snapshots(reference_profile, snapshots):
    """Map a 2017 reference profile to the snapshot period, matching day-of-week and hour."""

    mapped_profile = pd.Series(index=snapshots, dtype=float)

    for timestamp in snapshots:
        # Find hours in reference profile with same day-of-week and hour
        ref_candidates = reference_profile[
            (reference_profile.index.dayofweek == timestamp.dayofweek)
            & (reference_profile.index.hour == timestamp.hour)
        ]

        if len(ref_candidates) > 0:
            # Use mean of matching hours (typically ~52 weeks in a year)
            mapped_profile[timestamp] = ref_candidates.mean()
        else:
            # Fallback: use hour-of-day average
            ref_candidates = reference_profile[
                reference_profile.index.hour == timestamp.hour
            ]
            mapped_profile[timestamp] = (
                ref_candidates.mean()
                if len(ref_candidates) > 0
                else reference_profile.mean()
            )

    # CRITICAL: Renormalize to ensure it sums to 1
    mapped_profile = mapped_profile / mapped_profile.sum()

    return mapped_profile

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_energy_demand_per_node",
            clusters=48,
            planning_horizons=2030,
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
    nodal_df.to_csv(fn, float_format="%.2f")


    # Generate hourly electricity profiles
    logger.info("Creating hourly industry electricity demand profiles...")
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    nodal_electricity_profiles = create_nodal_electricity_profiles(
        nodal_df, nodal_production, sector_ratios, snapshots
    )

    # Check that hourly profiles match annual demand
    hourly_sum_twh = nodal_electricity_profiles.sum() / 1e6
    assert np.allclose(hourly_sum_twh, nodal_df["electricity"], rtol=1e-5), (
        f"Hourly profiles do not match annual demand.\nMax difference: {(hourly_sum_twh - nodal_df['electricity']).abs().max():.6f} TWh"
    )
    logger.info("✓ Hourly profiles verified to match annual demand")

    # Export hourly profiles
    fn_profiles = snakemake.output.industrial_electricity_demand_per_node_temporal
    nodal_electricity_profiles.to_csv(fn_profiles, float_format="%.2f")
    logger.info(f"Hourly electricity profiles saved to {fn_profiles}")
