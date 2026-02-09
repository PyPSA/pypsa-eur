# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build industry DSR checkpoint profile per FfE profile.

Produces a time series (snapshots x FfE profiles) with value 0 at checkpoint hours
and 1.0 elsewhere. Used as e_max_pu / e_min_pu for industry DSR Stores so that
demand is balanced within each period between checkpoints (like residential heat DSM).
Different sectors (FfE profiles) can have different restriction_time (operating hours).
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def build_industry_dsr_profile(snapshots, restriction_time):
    """
    Build DSR checkpoint profile per FfE profile (snapshots x profiles).

    At checkpoint hours for each profile, value is 0 (store must be empty);
    elsewhere value is 1.0. restriction_value is applied in prepare_sector_network.
    Hour of day is taken from snapshots (UTC).

    Parameters
    ----------
    snapshots : pd.DatetimeIndex
        Time index (e.g. n.snapshots).
    restriction_time : dict
        Keys = FfE profile names, values = list of hours (0-23) at which store must be empty.

    Returns
    -------
    pd.DataFrame
        Index = snapshots, columns = profile names, values in {0, 1}.
    """
    profiles = list(restriction_time.keys())
    out = pd.DataFrame(index=snapshots, columns=profiles, dtype=float)
    for profile in profiles:
        hours_check = restriction_time.get(profile, [])
        out[profile] = 1.0
        if hours_check:
            mask = out.index.hour.isin(hours_check)
            out.loc[mask, profile] = 0.0
    return out


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_dsr_profile",
            clusters=50,
            planning_horizons=2030,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Get snapshots and profile names from per-profile temporal demand file
    path = snakemake.input.industrial_electricity_demand_per_profile_temporal
    if isinstance(path, list):
        path = path[0]
    per_profile = pd.read_csv(path, index_col=0, parse_dates=True)
    snapshots = per_profile.index
    # Unique FfE profile names from columns "node|profile"
    profiles = sorted({c.split("|", 1)[1] for c in per_profile.columns})

    restriction_time = snakemake.params.restriction_time or {}
    technology_breakdown = snakemake.params.get("technology_breakdown", {})

    # Check if technology breakdown is provided
    use_technology_breakdown = bool(technology_breakdown)
    
    if use_technology_breakdown:
        # Technology-specific: build profiles for "profile|technology" keys
        logger.info("Building technology-specific DSR profiles")
        
        # Collect all technology keys from restriction_time
        tech_keys = set()
        profile_keys = set()
        
        for key in restriction_time.keys():
            if "|" in key:
                # Technology-specific key: "profile|technology"
                tech_keys.add(key)
                profile_part = key.split("|", 1)[0]
                profile_keys.add(profile_part)
            else:
                # Profile-level key
                profile_keys.add(key)
        
        # Build profiles for technology-specific keys
        tech_restriction_time = {k: restriction_time[k] for k in tech_keys if k in restriction_time}
        profile_restriction_time = {k: restriction_time[k] for k in profile_keys if k in restriction_time and k not in tech_keys}
        
        # Build technology-specific profiles
        all_columns = list(tech_keys) + list(profile_keys)
        dsr_profile = pd.DataFrame(index=snapshots, columns=all_columns, dtype=float)
        
        # Fill technology-specific columns
        for tech_key in tech_keys:
            hours_check = restriction_time.get(tech_key, [])
            dsr_profile[tech_key] = 1.0
            if hours_check:
                mask = dsr_profile.index.hour.isin(hours_check)
                dsr_profile.loc[mask, tech_key] = 0.0
        
        # Fill profile-level columns (for profiles without technology breakdown)
        for profile in profile_keys:
            if profile in tech_keys:
                continue  # Already handled as technology-specific
            hours_check = restriction_time.get(profile, [])
            dsr_profile[profile] = 1.0
            if hours_check:
                mask = dsr_profile.index.hour.isin(hours_check)
                dsr_profile.loc[mask, profile] = 0.0
        
        # Ensure all profiles present in data have a column (fill missing with 1.0)
        for p in profiles:
            if p not in dsr_profile.columns:
                dsr_profile[p] = 1.0
    else:
        # No technology_breakdown: build profile-level profiles (may not be used if no DSR stores are created)
        logger.info("Building profile-level DSR profiles (technology_breakdown not provided)")
        # Restrict to profiles present in the data
        restriction_time = {p: restriction_time[p] for p in profiles if p in restriction_time}
        if not restriction_time:
            logger.warning("No restriction_time keys match profile names; using all 1.0 profile.")
            dsr_profile = pd.DataFrame(1.0, index=snapshots, columns=profiles)
        else:
            dsr_profile = build_industry_dsr_profile(snapshots, restriction_time)
            # Ensure all profiles present in data have a column (fill missing with 1.0)
            for p in profiles:
                if p not in dsr_profile.columns:
                    dsr_profile[p] = 1.0

    dsr_profile.to_csv(snakemake.output.industrial_dsr_profile, float_format="%.4f")
    logger.info(f"Industry DSR profile saved to {snakemake.output.industrial_dsr_profile}")
