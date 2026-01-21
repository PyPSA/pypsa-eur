# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Prepare district heating subnode demand data for prepare_sector_network.

Extends energy/heat totals, district heat shares, and time-series data for subnodes.
"""

import logging

import geopandas as gpd
import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def _to_scalar(val):
    """
    Convert a pandas Series, DataFrame, or array-like to a scalar float.

    Parameters
    ----------
    val : scalar, pandas.Series, pandas.DataFrame, or array-like
        Value to convert

    Returns
    -------
    float
        Scalar float value, or 0.0 if input is None
    """
    if isinstance(val, pd.Series):
        val = val.iloc[0]
    if isinstance(val, pd.DataFrame):
        val = val.iloc[0, 0]
    if hasattr(val, "item"):
        val = val.item()
    return float(val) if val is not None else 0.0


def _get_row(df, idx):
    """
    Get a row from DataFrame as Series, handling duplicate indices.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame to extract row from
    idx : hashable
        Index label of the row to extract

    Returns
    -------
    pandas.Series
        Row data as Series. If index is duplicated, returns first occurrence.
    """
    result = df.loc[idx]
    return result.iloc[0] if isinstance(result, pd.DataFrame) else result


def extend_district_heat_share(
    district_heat_share: pd.DataFrame,
    subnodes: gpd.GeoDataFrame,
    pop_weighted_energy_totals: pd.DataFrame,
    pop_weighted_heat_totals: pd.DataFrame,
    heating_efficiencies: pd.DataFrame,
    district_heating_loss: float,
    space_heat_reduction: float = 0.0,
) -> pd.DataFrame:
    """
    Extend district heat share DataFrame to include subnodes.

    Scales subnode demands to match GIS data. If the sum of subnode demands
    exceeds the cluster's district heating demand, subnodes are scaled
    proportionally and the parent cluster's demand is set to zero.

    Parameters
    ----------
    district_heat_share : pandas.DataFrame
        District heat shares indexed by cluster name
    subnodes : geopandas.GeoDataFrame
        Subnodes with cluster, name, and yearly_heat_demand_MWh columns
    pop_weighted_energy_totals : pandas.DataFrame
        Population-weighted energy totals by cluster
    pop_weighted_heat_totals : pandas.DataFrame
        Population-weighted heat totals by cluster
    heating_efficiencies : pandas.DataFrame
        Heating efficiencies by country code
    district_heating_loss : float
        District heating network loss factor (e.g., 0.15 for 15% loss)
    space_heat_reduction : float, optional
        Fraction of space heating demand reduction, by default 0.0

    Returns
    -------
    pandas.DataFrame
        Extended district heat share with subnodes added
    """
    if len(subnodes) == 0:
        return district_heat_share.copy()

    extended_share = district_heat_share.copy()

    if "district fraction before subnodes" not in extended_share.columns:
        extended_share["district fraction before subnodes"] = extended_share[
            "district fraction of node"
        ].copy()

    extended_share["parent_node"] = extended_share.index.to_series()

    for cluster in subnodes["cluster"].unique():
        cluster_subnodes = subnodes[subnodes["cluster"] == cluster]

        if cluster not in pop_weighted_energy_totals.index:
            logger.warning(f"Cluster {cluster} not in energy totals, skipping")
            continue

        ct = cluster[:2]
        cluster_energy = _get_row(pop_weighted_energy_totals, cluster)
        cluster_heat = _get_row(pop_weighted_heat_totals, cluster)

        # Calculate useful heat: space from heat_totals, water from energy_totals
        total_useful_heat_mwh = 0.0
        for sector in ["residential", "services"]:
            for use in ["water", "space"]:
                col = f"total {sector} {use}"
                eff_col = f"total {sector} {use} efficiency"

                # Get efficiency
                if (
                    eff_col in heating_efficiencies.columns
                    and ct in heating_efficiencies.index
                ):
                    eff = _to_scalar(heating_efficiencies.loc[ct, eff_col])
                else:
                    eff = 1.0

                # Space from heat_totals, water from energy_totals
                if use == "space":
                    if col in cluster_heat.index:
                        val = _to_scalar(cluster_heat[col])
                        total_useful_heat_mwh += (
                            val * eff * (1 - space_heat_reduction) * 1e6
                        )
                else:  # water
                    if col in cluster_energy.index:
                        val = _to_scalar(cluster_energy[col])
                        total_useful_heat_mwh += val * eff * 1e6

        if total_useful_heat_mwh == 0.0:
            logger.warning(
                f"Zero heat demand for cluster {cluster}, skipping its subnodes"
            )
            continue

        # Get cluster's current district heating demand (before subnodes)
        cluster_dist_fraction = extended_share.loc[cluster, "district fraction of node"]
        cluster_dh_demand_mwh = (
            cluster_dist_fraction * total_useful_heat_mwh * (1 + district_heating_loss)
        )

        # Sum of all subnode demands for this cluster
        total_subnode_demand = cluster_subnodes["yearly_heat_demand_MWh"].sum()

        # Check if aggregate subnode demand exceeds cluster's district heating demand
        if total_subnode_demand > cluster_dh_demand_mwh:
            scale_factor = cluster_dh_demand_mwh / total_subnode_demand
            logger.warning(
                f"Cluster {cluster}: Sum of subnode demands ({total_subnode_demand / 1e6:.3f} TWh) "
                f"exceeds cluster's district heating demand ({cluster_dh_demand_mwh / 1e6:.3f} TWh). "
                f"Scaling subnode demands by {scale_factor:.2%}. Parent node demand set to zero."
            )
            # After scaling, all district heating demand goes to subnodes
            remaining_cluster_demand_mwh = 0.0
        else:
            scale_factor = 1.0
            remaining_cluster_demand_mwh = cluster_dh_demand_mwh - total_subnode_demand

        # Process each subnode in the cluster
        for _, subnode in cluster_subnodes.iterrows():
            name = subnode["name"]
            demand_mwh = subnode["yearly_heat_demand_MWh"] * scale_factor

            # Calculate district fraction for subnode
            # Note: For subnodes, district fraction should be 1.0 because:
            # 1. Subnodes are 100% urban central heat (no rural or decentral portion)
            # 2. The subnode's demand is already embedded in pop_weighted_energy_totals
            #    (scaled by extend_pop_weighted_energy_totals)
            # 3. Using dist_fraction < 1.0 would double-scale the demand
            #
            # The "original district heat share" stores the proportion of parent's
            # district heating demand that this subnode represents
            original_share = demand_mwh / (
                total_useful_heat_mwh * (1 + district_heating_loss)
            )

            # Add entry for subnode with parent_node mapping for resource bus lookups
            extended_share.loc[name] = {
                "original district heat share": original_share,
                "district fraction of node": 1.0,  # Subnodes are 100% urban central
                "urban fraction": 1.0,  # Subnodes are 100% urban by definition
                "district fraction before subnodes": 0.0,  # New entry, no pre-subnode value
                "parent_node": cluster,  # Parent cluster for resource bus lookups
            }

        # Update parent cluster's district fraction
        new_cluster_fraction = remaining_cluster_demand_mwh / (
            total_useful_heat_mwh * (1 + district_heating_loss)
        )
        extended_share.loc[cluster, "district fraction of node"] = new_cluster_fraction

    return extended_share


def extend_pop_weighted_energy_totals(
    pop_weighted_energy_totals: pd.DataFrame,
    pop_weighted_heat_totals: pd.DataFrame,
    subnodes: gpd.GeoDataFrame,
    heating_efficiencies: pd.DataFrame,
    district_heating_loss: float,
    space_heat_reduction: float = 0.0,
) -> pd.DataFrame:
    """
    Extend population-weighted energy totals to include subnodes.

    Creates entries for each subnode by scaling the parent cluster's energy
    profile based on the subnode's GIS demand data. Uses heat_totals for space
    heating (due to .update() in prepare_sector_network) and energy_totals for
    water heating when computing useful heat.

    Parameters
    ----------
    pop_weighted_energy_totals : pandas.DataFrame
        Population-weighted energy totals indexed by cluster name
    pop_weighted_heat_totals : pandas.DataFrame
        Population-weighted heat totals indexed by cluster name
    subnodes : geopandas.GeoDataFrame
        Subnodes with cluster, name, and yearly_heat_demand_MWh columns
    heating_efficiencies : pandas.DataFrame
        Heating efficiencies by country code
    district_heating_loss : float
        District heating network loss factor
    space_heat_reduction : float, optional
        Fraction of space heating demand reduction, by default 0.0

    Returns
    -------
    pandas.DataFrame
        Extended energy totals with subnodes added
    """
    if len(subnodes) == 0:
        return pop_weighted_energy_totals.copy()

    extended_totals = pop_weighted_energy_totals.copy()

    for _, subnode in subnodes.iterrows():
        cluster = subnode["cluster"]
        name = subnode["name"]
        demand_mwh = subnode["yearly_heat_demand_MWh"]

        if cluster not in extended_totals.index:
            logger.warning(f"Cluster {cluster} not found, skipping subnode {name}")
            continue

        ct = cluster[:2]
        subnode_totals = _get_row(extended_totals, cluster).copy()

        # Calculate useful heat: space from heat_totals, water from energy_totals
        cluster_useful_heat = 0.0
        cluster_heat_row = _get_row(pop_weighted_heat_totals, cluster)
        cluster_energy_row = _get_row(pop_weighted_energy_totals, cluster)

        for sector in ["residential", "services"]:
            for use in ["water", "space"]:
                col = f"total {sector} {use}"
                eff_col = f"total {sector} {use} efficiency"

                # Get efficiency
                if (
                    eff_col in heating_efficiencies.columns
                    and ct in heating_efficiencies.index
                ):
                    eff = _to_scalar(heating_efficiencies.loc[ct, eff_col])
                else:
                    eff = 1.0

                # Space comes from heat_totals (after .update()), water from energy_totals
                if use == "space":
                    if col in cluster_heat_row.index:
                        val = _to_scalar(cluster_heat_row[col])
                        cluster_useful_heat += (
                            val * eff * (1 - space_heat_reduction) * 1e6
                        )
                else:  # water
                    if col in cluster_energy_row.index:
                        val = _to_scalar(cluster_energy_row[col])
                        cluster_useful_heat += val * eff * 1e6

        if cluster_useful_heat == 0.0:
            continue

        # Scale factor: subnode_useful_heat * (1 + dh_loss) = GIS_demand
        scale_factor = demand_mwh / (cluster_useful_heat * (1 + district_heating_loss))

        # Scale all heat-related columns (water, space, heat)
        heat_cols = [
            col
            for col in subnode_totals.index
            if any(x in col for x in ["water", "space", "heat"])
        ]
        for col in heat_cols:
            cluster_val = _to_scalar(extended_totals.loc[cluster, col])
            subnode_totals[col] = cluster_val * scale_factor

        # Add subnode entry
        extended_totals.loc[name] = subnode_totals

        # District fraction is reduced in extend_district_heat_share, not here

    return extended_totals


def extend_pop_weighted_heat_totals(
    pop_weighted_heat_totals: pd.DataFrame,
    pop_weighted_energy_totals: pd.DataFrame,
    subnodes: gpd.GeoDataFrame,
    heating_efficiencies: pd.DataFrame,
    district_heating_loss: float,
    space_heat_reduction: float = 0.0,
) -> pd.DataFrame:
    """
    Extend population-weighted heat totals to include subnodes.

    Creates entries for each subnode by scaling the parent cluster's heat
    profile. Uses the same scaling methodology as extend_pop_weighted_energy_totals
    since heat_totals.update() overwrites space columns in energy_totals.

    Parameters
    ----------
    pop_weighted_heat_totals : pandas.DataFrame
        Population-weighted heat totals indexed by cluster name
    pop_weighted_energy_totals : pandas.DataFrame
        Population-weighted energy totals indexed by cluster name
    subnodes : geopandas.GeoDataFrame
        Subnodes with cluster, name, and yearly_heat_demand_MWh columns
    heating_efficiencies : pandas.DataFrame
        Heating efficiencies by country code
    district_heating_loss : float
        District heating network loss factor
    space_heat_reduction : float, optional
        Fraction of space heating demand reduction, by default 0.0

    Returns
    -------
    pandas.DataFrame
        Extended heat totals with subnodes added
    """
    if len(subnodes) == 0:
        return pop_weighted_heat_totals.copy()

    extended_totals = pop_weighted_heat_totals.copy()

    for _, subnode in subnodes.iterrows():
        cluster = subnode["cluster"]
        name = subnode["name"]
        demand_mwh = subnode["yearly_heat_demand_MWh"]
        ct = cluster[:2]

        if cluster not in extended_totals.index:
            logger.warning(f"Cluster {cluster} not in heat totals, skipping {name}")
            continue

        # Copy cluster's heat profile
        subnode_totals = _get_row(extended_totals, cluster).copy()

        # Calculate cluster's total useful heat demand using the SAME formula as
        # extend_pop_weighted_energy_totals and prepare_sector_network.build_heat_demand
        # This includes: (space * eff * (1-reduction) + water * eff)
        cluster_useful_heat = 0.0

        # Get energy totals for water heating (not in heat_totals)
        cluster_energy = _get_row(pop_weighted_energy_totals, cluster)

        for sector in ["residential", "services"]:
            for use in ["water", "space"]:
                col = f"total {sector} {use}"
                eff_col = f"total {sector} {use} efficiency"

                # Get efficiency
                if (
                    eff_col in heating_efficiencies.columns
                    and ct in heating_efficiencies.index
                ):
                    eff = _to_scalar(heating_efficiencies.loc[ct, eff_col])
                else:
                    eff = 1.0

                # Space comes from heat_totals, water from energy_totals
                if use == "space":
                    if col in subnode_totals.index:
                        val = _to_scalar(subnode_totals[col])
                        cluster_useful_heat += (
                            val * eff * (1 - space_heat_reduction) * 1e6
                        )
                else:  # water
                    if col in cluster_energy.index:
                        val = _to_scalar(cluster_energy[col])
                        cluster_useful_heat += val * eff * 1e6

        if cluster_useful_heat == 0.0:
            continue

        # Scale factor: subnode_useful_heat * (1 + dh_loss) = GIS_demand
        # So: scale * cluster_useful_heat * (1 + dh_loss) = demand_mwh
        scale_factor = demand_mwh / (cluster_useful_heat * (1 + district_heating_loss))

        # Scale all columns in heat_totals
        for col in subnode_totals.index:
            cluster_val = _to_scalar(extended_totals.loc[cluster, col])
            subnode_totals[col] = cluster_val * scale_factor

        # Add subnode entry
        extended_totals.loc[name] = subnode_totals

        # NOTE: We do NOT reduce the cluster's values here!
        # The cluster's district fraction is reduced in extend_district_heat_share
        # to account for the demand moving to subnodes.

    return extended_totals


def extend_industrial_demand(
    industrial_demand: pd.DataFrame,
    subnodes: gpd.GeoDataFrame,
    district_heat_share: pd.DataFrame,
    pop_weighted_energy_totals: pd.DataFrame,
    heating_efficiencies: pd.DataFrame,
    district_heating_loss: float,
) -> pd.DataFrame:
    """
    Extend industrial demand to include subnodes based on parent cluster's LT heat ratio.

    Allocates low-temperature heat for industry to subnodes proportionally based on
    the ratio of LT industry heat to total district heating demand in the parent cluster.

    Parameters
    ----------
    industrial_demand : pandas.DataFrame
        Industrial demand by sector indexed by cluster name
    subnodes : geopandas.GeoDataFrame
        Subnodes with cluster, name, and yearly_heat_demand_MWh columns
    district_heat_share : pandas.DataFrame
        District heat shares indexed by cluster name
    pop_weighted_energy_totals : pandas.DataFrame
        Population-weighted energy totals indexed by cluster name
    heating_efficiencies : pandas.DataFrame
        Heating efficiencies by country code
    district_heating_loss : float
        District heating network loss factor

    Returns
    -------
    pandas.DataFrame
        Extended industrial demand with subnodes added and parent values reduced
    """
    if len(subnodes) == 0:
        return industrial_demand.copy()

    extended_demand = industrial_demand.copy()

    # Group subnodes by cluster
    for cluster in subnodes["cluster"].unique():
        cluster_subnodes = subnodes[subnodes["cluster"] == cluster]

        if cluster not in industrial_demand.index:
            logger.warning(f"Cluster {cluster} not in industrial demand, skipping")
            continue

        if cluster not in pop_weighted_energy_totals.index:
            logger.warning(f"Cluster {cluster} not in energy totals, skipping")
            continue

        if cluster not in district_heat_share.index:
            logger.warning(f"Cluster {cluster} not in district heat share, skipping")
            continue

        ct = cluster[:2]

        # Get cluster's LT industry heat demand (TWh)
        cluster_lt_heat_twh = _to_scalar(
            industrial_demand.loc[cluster, "low-temperature heat"]
        )

        if cluster_lt_heat_twh == 0.0:
            # No LT industry heat in this cluster - subnodes only get urban central heat
            for _, subnode in cluster_subnodes.iterrows():
                name = subnode["name"]
                # Create subnode entry with zeros for all columns
                extended_demand.loc[name] = 0.0
            continue

        # Calculate cluster's total heat demand (MWh)
        cluster_energy = _get_row(pop_weighted_energy_totals, cluster)
        total_heat_demand_mwh = 0.0
        for sector in ["residential", "services"]:
            for use in ["water", "space"]:
                col = f"total {sector} {use}"
                if col in cluster_energy.index:
                    eff_col = f"total {sector} {use} efficiency"
                    if (
                        eff_col in heating_efficiencies.columns
                        and ct in heating_efficiencies.index
                    ):
                        eff = _to_scalar(heating_efficiencies.loc[ct, eff_col])
                    else:
                        eff = 1.0
                    val = _to_scalar(cluster_energy[col])
                    total_heat_demand_mwh += val * eff * 1e6

        if total_heat_demand_mwh == 0.0:
            continue

        # Get cluster's district heating demand (MWh) including losses
        cluster_dist_fraction = _to_scalar(
            district_heat_share.loc[cluster, "district fraction of node"]
        )
        cluster_dh_demand_mwh = (
            cluster_dist_fraction * total_heat_demand_mwh * (1 + district_heating_loss)
        )

        if cluster_dh_demand_mwh == 0.0:
            continue

        # Compute ratio of LT industry heat to total district heating demand
        cluster_lt_heat_mwh = cluster_lt_heat_twh * 1e6
        lt_industry_ratio = cluster_lt_heat_mwh / cluster_dh_demand_mwh

        # Cap at 1.0 (LT industry heat cannot exceed district heating demand)
        lt_industry_ratio = min(lt_industry_ratio, 1.0)

        logger.info(
            f"Cluster {cluster}: LT industry heat ratio = {lt_industry_ratio:.2%} "
            f"(LT={cluster_lt_heat_twh:.3f} TWh, DH={cluster_dh_demand_mwh / 1e6:.3f} TWh)"
        )

        # Sum of subnode demands for capping
        total_subnode_demand = cluster_subnodes["yearly_heat_demand_MWh"].sum()

        # Scale factor if subnodes exceed cluster district heating demand
        if total_subnode_demand > cluster_dh_demand_mwh:
            scale_factor = cluster_dh_demand_mwh / total_subnode_demand
        else:
            scale_factor = 1.0

        # Process each subnode
        total_subnode_lt_twh = 0.0
        for _, subnode in cluster_subnodes.iterrows():
            name = subnode["name"]
            subnode_dh_mwh = subnode["yearly_heat_demand_MWh"] * scale_factor

            # Split subnode district heating demand: part to LT industry
            subnode_lt_mwh = subnode_dh_mwh * lt_industry_ratio
            subnode_lt_twh = subnode_lt_mwh / 1e6

            # Create subnode entry with LT heat for industry
            # Initialize with zeros for all columns
            extended_demand.loc[name] = 0.0
            extended_demand.loc[name, "low-temperature heat"] = subnode_lt_twh

            total_subnode_lt_twh += subnode_lt_twh

        # Reduce parent cluster's LT industry heat
        remaining_lt_twh = max(0, cluster_lt_heat_twh - total_subnode_lt_twh)
        extended_demand.loc[cluster, "low-temperature heat"] = remaining_lt_twh

        logger.info(
            f"Cluster {cluster}: Moved {total_subnode_lt_twh:.3f} TWh LT heat to subnodes, "
            f"remaining {remaining_lt_twh:.3f} TWh"
        )

    return extended_demand


def extend_hourly_heat_demand(
    hourly_heat_demand: xr.Dataset,
    subnodes: gpd.GeoDataFrame,
) -> xr.Dataset:
    """
    Extend hourly heat demand Dataset to include subnodes.

    Each subnode inherits the hourly heat profile of its parent cluster.

    Parameters
    ----------
    hourly_heat_demand : xarray.Dataset
        Hourly heat demand profiles with 'node' coordinate
    subnodes : geopandas.GeoDataFrame
        Subnodes with name and cluster columns

    Returns
    -------
    xarray.Dataset
        Extended hourly heat demand with subnodes added
    """
    if len(subnodes) == 0:
        return hourly_heat_demand

    # The dataset uses 'node' as the coordinate name
    base_nodes = pd.Index(hourly_heat_demand.coords["node"].values)

    # Build list of new datasets for subnodes
    new_datasets = []
    for _, subnode in subnodes.iterrows():
        name = subnode["name"]
        cluster = subnode["cluster"]
        if cluster in base_nodes:
            # Copy parent's hourly profile for all variables
            parent_data = hourly_heat_demand.sel(node=cluster)
            new_datasets.append(parent_data.assign_coords(node=name))

    if not new_datasets:
        return hourly_heat_demand

    # Concatenate
    new_data = xr.concat(new_datasets, dim="node")
    result = xr.concat([hourly_heat_demand, new_data], dim="node")

    # Ensure clean coordinate
    result_nodes = [str(n) for n in result.coords["node"].values]
    result = result.assign_coords(node=result_nodes)

    return result


def extend_heat_dsm_profile(
    dsm_profile: pd.DataFrame,
    subnodes: gpd.GeoDataFrame,
) -> pd.DataFrame:
    """
    Extend heat demand-side management profile to include subnodes.

    Each subnode inherits the DSM profile of its parent cluster.

    Parameters
    ----------
    dsm_profile : pandas.DataFrame
        Heat DSM profiles with cluster names as columns
    subnodes : geopandas.GeoDataFrame
        Subnodes with name and cluster columns

    Returns
    -------
    pandas.DataFrame
        Extended DSM profile with subnodes added as new columns
    """
    if len(subnodes) == 0:
        return dsm_profile.copy()

    extended = dsm_profile.copy()
    base_nodes = dsm_profile.columns

    for _, subnode in subnodes.iterrows():
        name = subnode["name"]
        cluster = subnode["cluster"]
        if cluster in base_nodes:
            extended[name] = dsm_profile[cluster]

    return extended


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_district_heating_subnodes",
            clusters=48,
            planning_horizons=2050,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load configuration
    district_heating_loss = snakemake.params.get("district_heating_loss", 0.15)
    energy_totals_year = int(snakemake.params.get("energy_totals_year", 2019))

    # Get space heat reduction factor
    investment_year = int(snakemake.wildcards.planning_horizons)
    reduce_space_heat = snakemake.params.get("reduce_space_heat_exogenously", False)
    if reduce_space_heat:
        reduction_factors = snakemake.params.get(
            "reduce_space_heat_exogenously_factor", {}
        )
        # Get the value for the current investment year, or closest available
        if isinstance(reduction_factors, dict):
            if investment_year in reduction_factors:
                space_heat_reduction = reduction_factors[investment_year]
            else:
                # Find the closest year
                years = sorted(reduction_factors.keys())
                closest_year = min(years, key=lambda y: abs(y - investment_year))
                space_heat_reduction = reduction_factors[closest_year]
        else:
            space_heat_reduction = float(reduction_factors)
    else:
        space_heat_reduction = 0.0

    logger.info(f"Using space heat reduction factor: {space_heat_reduction}")
    logger.info(f"Using energy totals year: {energy_totals_year}")

    # Load input data
    subnodes = gpd.read_file(snakemake.input.dh_subnodes)

    if len(subnodes) == 0:
        logger.info("No subnodes found, creating pass-through outputs")
    else:
        logger.info(f"Processing {len(subnodes)} subnodes")

    pop_layout = pd.read_csv(snakemake.input.pop_layout, index_col=0)
    district_heat_share = pd.read_csv(snakemake.input.district_heat_share, index_col=0)
    pop_weighted_energy_totals = pd.read_csv(
        snakemake.input.pop_weighted_energy_totals, index_col=0
    )
    pop_weighted_heat_totals = pd.read_csv(
        snakemake.input.pop_weighted_heat_totals, index_col=0
    )
    # Load heating efficiencies for the energy_totals_year (same as prepare_sector_network)
    heating_efficiencies = pd.read_csv(
        snakemake.input.heating_efficiencies, index_col=[1, 0]
    ).loc[energy_totals_year]
    industrial_demand = pd.read_csv(snakemake.input.industrial_demand, index_col=0)

    # Load time-series data
    hourly_heat_demand = xr.open_dataset(snakemake.input.hourly_heat_demand)
    dsm_profile = pd.read_csv(
        snakemake.input.heat_dsm_profile, header=1, index_col=0, parse_dates=True
    )

    # Extend district heat share (handles demand capping and adds parent_node column)
    district_heat_share_extended = extend_district_heat_share(
        district_heat_share=district_heat_share,
        subnodes=subnodes,
        pop_weighted_energy_totals=pop_weighted_energy_totals,
        pop_weighted_heat_totals=pop_weighted_heat_totals,
        heating_efficiencies=heating_efficiencies,
        district_heating_loss=district_heating_loss,
        space_heat_reduction=space_heat_reduction,
    )

    # Extend population-weighted energy totals
    pop_weighted_energy_totals_extended = extend_pop_weighted_energy_totals(
        pop_weighted_energy_totals=pop_weighted_energy_totals,
        pop_weighted_heat_totals=pop_weighted_heat_totals,
        subnodes=subnodes,
        heating_efficiencies=heating_efficiencies,
        district_heating_loss=district_heating_loss,
        space_heat_reduction=space_heat_reduction,
    )

    # Extend population-weighted heat totals (same scaling as energy_totals)
    pop_weighted_heat_totals_extended = extend_pop_weighted_heat_totals(
        pop_weighted_heat_totals=pop_weighted_heat_totals,
        pop_weighted_energy_totals=pop_weighted_energy_totals,
        subnodes=subnodes,
        heating_efficiencies=heating_efficiencies,
        district_heating_loss=district_heating_loss,
        space_heat_reduction=space_heat_reduction,
    )

    # Extend industrial demand (LT heat for industry in subnodes)
    industrial_demand_extended = extend_industrial_demand(
        industrial_demand=industrial_demand,
        subnodes=subnodes,
        district_heat_share=district_heat_share,  # Use ORIGINAL
        pop_weighted_energy_totals=pop_weighted_energy_totals,
        heating_efficiencies=heating_efficiencies,
        district_heating_loss=district_heating_loss,
    )

    # Extend time-series data
    hourly_heat_demand_extended = extend_hourly_heat_demand(
        hourly_heat_demand, subnodes
    )
    dsm_profile_extended = extend_heat_dsm_profile(dsm_profile, subnodes)

    # Save outputs
    district_heat_share_extended.to_csv(snakemake.output.district_heat_share_subnodes)
    logger.info(
        f"Saved extended district heat share to {snakemake.output.district_heat_share_subnodes}"
    )

    pop_weighted_energy_totals_extended.to_csv(
        snakemake.output.pop_weighted_energy_totals_subnodes
    )
    pop_weighted_heat_totals_extended.to_csv(
        snakemake.output.pop_weighted_heat_totals_subnodes
    )

    industrial_demand_extended.to_csv(snakemake.output.industrial_demand_subnodes)
    logger.info(
        f"Saved extended industrial demand to {snakemake.output.industrial_demand_subnodes}"
    )

    # Save extended time-series
    hourly_heat_demand_extended.to_netcdf(snakemake.output.hourly_heat_demand_subnodes)
    logger.info(
        f"Saved extended hourly heat demand to {snakemake.output.hourly_heat_demand_subnodes}"
    )

    dsm_profile_extended.to_csv(snakemake.output.heat_dsm_profile_subnodes)
    logger.info(
        f"Saved extended DSM profile to {snakemake.output.heat_dsm_profile_subnodes}"
    )

    # Log summary
    if len(subnodes) > 0:
        total_subnode_demand = subnodes["yearly_heat_demand_MWh"].sum()
        logger.info(
            f"Processed {len(subnodes)} subnodes with total demand of {total_subnode_demand / 1e6:.2f} TWh/a"
        )
