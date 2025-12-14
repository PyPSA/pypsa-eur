# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Compose network by combining all electricity and sector components.

This script combines functionality from add_electricity.py, prepare_sector_network.py,
add_existing_baseyear.py, add_brownfield.py, prepare_network.py, and
prepare_perfect_foresight.py into a unified composition step for the streamlined
workflow.

Each major section is tagged with a subtitle indicating which legacy script it
replaces (for example, 'ELECTRICITY COMPONENTS (from add_electricity.py)').

All function calls are made directly in the main section without defining
additional functions, following the additive approach of importing existing
functions and calling them in sequence.
"""

import logging
import os
from collections.abc import Mapping

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

# Import all required functions from existing scripts
from scripts._helpers import (
    PYPSA_V1,
    configure_logging,
    sanitize_custom_columns,
    set_scenario_config,
    update_config_from_wildcards,
    update_p_nom_max,
)
from scripts.add_brownfield import (
    add_brownfield,
    adjust_renewable_profiles,
    disable_grid_expansion_if_limit_hit,
    update_dynamic_ptes_capacity,
    update_heat_pump_efficiency,
)
from scripts.add_electricity import (
    apply_variable_renewable_lifetimes,
    attach_conventional_generators,
    attach_GEM_renewables,
    attach_hydro,
    attach_load,
    attach_storageunits,
    attach_stores,
    attach_wind_and_solar,
    estimate_renewable_capacities,
    finalize_electricity_network,
    load_and_aggregate_powerplants,
    load_costs,
    remove_non_power_buses,
    restrict_electricity_components,
    sanitize_carriers,
    sanitize_locations,
    set_transmission_costs,
)
from scripts.add_existing_baseyear import (
    add_build_year_to_new_assets,
    add_heating_capacities_installed_before_baseyear,
    add_power_capacities_installed_before_baseyear,
)
from scripts.co2_budget import bound_value_for_horizon, co2_budget_for_horizon
from scripts.prepare_network import (
    add_co2limit,
    add_dynamic_emission_prices,
    add_emission_prices,
    add_gaslimit,
    apply_time_segmentation,
    average_every_nhours,
    cap_transmission_capacity,
    set_transmission_limit,
)
from scripts.prepare_perfect_foresight import (
    adjust_stores_for_perfect_foresight,
    apply_investment_period_weightings,
    update_heat_pump_efficiency_for_horizon,
)
from scripts.prepare_sector_network import (
    add_agriculture,
    add_allam_gas,
    add_ammonia,
    add_aviation,
    add_biomass,
    add_carrier_buses,
    add_co2_network,
    add_co2_tracking,
    add_dac,
    add_eu_bus,
    add_generation,
    add_heat,
    add_industry,
    add_land_transport,
    add_methanol,
    add_shipping,
    add_storage_and_grids,
    add_waste_heat,
    cluster_heat_buses,
    co2_emissions_year,
    decentral,
    define_spatial,
    remove_h2_network,
    set_temporal_aggregation,
)

logger = logging.getLogger(__name__)


def extend_snapshot_multiindex(
    existing_snapshots: pd.MultiIndex,
    new_timesteps: pd.Index,
    period: int,
) -> pd.MultiIndex:
    """
    Extend a MultiIndex with period & timestep levels by appending new timesteps for a given period.

    This function takes an existing multiindex with (period, timestep) levels and adds
    a new set of timesteps for a specified period, returning the combined multiindex.

    Parameters
    ----------
    existing_snapshots : pd.MultiIndex
        Existing multiindex with levels ['period', 'timestep']
    new_timesteps : pd.Index
        Single-level index of timesteps to append
    period : int
        Investment period label for the new timesteps

    Returns
    -------
    pd.MultiIndex
        Combined multiindex with both existing and new snapshots

    Raises
    ------
    TypeError
        If existing_snapshots is not a MultiIndex or new_timesteps is a MultiIndex
    ValueError
        If existing_snapshots doesn't have the required ['period', 'timestep'] levels

    Examples
    --------
    >>> existing = pd.MultiIndex.from_product(
    ...     [[2030], pd.date_range('2030-01-01', periods=2, freq='h')],
    ...     names=['period', 'timestep']
    ... )
    >>> new = pd.date_range('2040-01-01', periods=2, freq='h')
    >>> extended = extend_snapshot_multiindex(existing, new, 2040)
    >>> len(extended)
    4
    >>> list(extended.get_level_values('period').unique())
    [2030, 2040]
    """
    # Validate existing_snapshots
    if not isinstance(existing_snapshots, pd.MultiIndex):
        raise TypeError("existing_snapshots must be a MultiIndex")

    if list(existing_snapshots.names) != ["period", "timestep"]:
        raise ValueError(
            f"existing_snapshots must have levels ['period', 'timestep'], "
            f"but has {existing_snapshots.names}"
        )

    # Validate new_timesteps
    if isinstance(new_timesteps, pd.MultiIndex):
        raise TypeError("new_timesteps must be a single-level Index, not a MultiIndex")

    # Create MultiIndex for new period
    new_snapshots = pd.MultiIndex.from_product(
        [[period], new_timesteps], names=["period", "timestep"]
    )

    # Concatenate by converting to tuples and back
    combined = pd.MultiIndex.from_tuples(
        list(existing_snapshots) + list(new_snapshots), names=["period", "timestep"]
    )

    return combined


def concatenate_network_with_previous(
    n_previous: pypsa.Network,
    n_current: pypsa.Network,
    current_horizon: int,
) -> pypsa.Network:
    """
    Concatenate current horizon network with previously composed multi-period network.

    For perfect foresight optimization, this function incrementally builds up a
    multi-period network by adding each planning horizon to the network from
    previous horizons.

    Parameters
    ----------
    n_previous : pypsa.Network
        Previously composed network containing one or more investment periods
    n_current : pypsa.Network
        Network for current planning horizon (single period)
    current_horizon : int
        Current planning horizon year

    Returns
    -------
    pypsa.Network
        Combined multi-period network

    Notes
    -----
    - The previous network may already contain multiple investment periods
    - Static components are merged without duplication
    - Time-varying data is combined with multi-indexed snapshots
    - Investment periods and weightings are recalculated
    """
    logger.info(
        f"Concatenating network for horizon {current_horizon} with previous periods"
    )

    n = n_previous.copy()

    if n.investment_periods.empty:
        raise ValueError("Previous network must have investment periods")
    if n_current.investment_periods.empty:
        raise ValueError("Current network must have investment periods")

    # Extend snapshots using our helper function
    combined_snapshots = list(n.snapshots) + list(n_current.snapshots)
    extended_snapshots = pd.MultiIndex.from_tuples(combined_snapshots)

    # Set snapshots - PyPSA will automatically update investment_periods
    # from the first level of the MultiIndex
    n.set_snapshots(extended_snapshots)
    n.snapshot_weightings.loc[current_horizon] = n_current.snapshot_weightings
    n.set_investment_periods(list(n.investment_periods))

    for c_current in n_current.components.values():
        c = n.c[c_current.name]
        existing_comps = c.static.index
        overlap_comps = c_current.static.index.intersection(existing_comps)
        new_comps = c_current.static.index.difference(existing_comps)
        new_static = c_current.static.loc[new_comps]
        n.add(c_current.name, new_static.index, **new_static)

        for attr, df in c_current.dynamic.items():
            if df.empty:
                continue

            default = c.attrs.default[attr]
            c.dynamic[attr].loc[current_horizon, df.columns] = df.values
            c.dynamic[attr] = c.dynamic[attr].fillna(default)

            # static values from different horizon that are not equal have to be casted to dynamic
            overlap_non_equal_static_only = c_current.static.loc[
                overlap_comps, attr
            ].ne(c.static.loc[overlap_comps, attr]) & ~overlap_comps.isin(
                c.dynamic[attr].columns
            )
            if overlap_non_equal_static_only.any():
                if PYPSA_V1:
                    casted = c._as_dynamic(
                        attr, n.snapshots, overlap_non_equal_static_only
                    )
                else:
                    casted = n.get_switchable_as_dense(
                        c.name, attr, inds=overlap_non_equal_static_only
                    )
                casted.loc[current_horizon] = c_current.static.loc[overlap_comps]

    n.meta = {**n_previous.meta, **n_current.meta}

    # Verify investment periods match snapshots first level before returning
    snapshot_periods = list(n.snapshots.get_level_values("period").unique())
    investment_periods_list = list(n.investment_periods)

    if snapshot_periods != investment_periods_list:
        logger.warning(
            f"Investment periods mismatch detected: snapshots have {snapshot_periods}, "
            f"but investment_periods is {investment_periods_list}. Synchronizing..."
        )
        # Force synchronization by setting investment_periods to match snapshots
        # Use set_investment_periods() to ensure proper validation and storage
        n.set_investment_periods(snapshot_periods)

    logger.info(
        f"Successfully concatenated network: {len(n.investment_periods)} investment periods"
    )
    return n


def validate_network_state(
    n: pypsa.Network, stage_name: str, expected_components: dict[str, int] | None = None
) -> dict[str, int]:
    """
    Validate network state and log component counts after a stage.

    Parameters
    ----------
    n : pypsa.Network
        Network to validate
    stage_name : str
        Name of the stage for logging
    expected_components : dict[str, int] | None
        Optional dict of {component_type: minimum_count} to validate

    Returns
    -------
    dict[str, int]
        Current component counts for reference

    Raises
    ------
    ValueError
        If expected component counts are not met
    """
    counts = {
        "buses": len(n.buses),
        "generators": len(n.generators),
        "loads": len(n.loads),
        "lines": len(n.lines),
        "links": len(n.links),
        "stores": len(n.stores),
        "storage_units": len(n.storage_units),
    }

    logger.debug(
        f"Network state after {stage_name}: "
        f"{counts['buses']} buses, "
        f"{counts['generators']} generators, "
        f"{counts['loads']} loads, "
        f"{counts['lines']} lines, "
        f"{counts['links']} links, "
        f"{counts['stores']} stores, "
        f"{counts['storage_units']} storage_units"
    )

    if expected_components:
        for component, min_count in expected_components.items():
            actual = counts.get(component, 0)
            if actual < min_count:
                raise ValueError(
                    f"After {stage_name}: expected at least {min_count} {component}, "
                    f"but found {actual}"
                )

    return counts


def adjust_renewable_capacity_limits(
    n: pypsa.Network, horizon: str, renewable_carriers: list[str]
) -> None:
    """
    Adjust renewable capacity limits by subtracting existing capacities from previous horizons.

    For each renewable carrier, this function sums the capacity of non-extendable
    generators (representing existing capacities from previous planning horizons)
    and subtracts that value from the p_nom_max of extendable generators in the
    current horizon. This ensures that the technical potential is properly adjusted
    for brownfield scenarios.

    Parameters
    ----------
    n : pypsa.Network
        Network containing renewable generators
    horizon : str
        The current planning horizon year as string
    renewable_carriers : list[str]
        List of renewable carrier names from config

    Notes
    -----
    Modifies n.generators["p_nom_max"] in-place.
    Issues a warning if existing capacities exceed technical potential.
    Clips p_nom_max to non-negative values.
    """
    for carrier in renewable_carriers:
        ext_i = (n.generators.carrier == carrier) & ~n.generators.p_nom_extendable
        grouper = n.generators.loc[ext_i].index.str.replace(
            f" {carrier}.*$", "", regex=True
        )
        existing = n.generators.loc[ext_i, "p_nom"].groupby(grouper).sum()
        existing.index += f" {carrier}-{horizon}"
        n.generators.loc[existing.index, "p_nom_max"] -= existing

    # Check if existing capacities are larger than technical potential
    existing_large = n.generators[
        n.generators["p_nom_min"] > n.generators["p_nom_max"]
    ].index
    if len(existing_large):
        logger.warning(
            f"Existing capacities larger than technical potential for {existing_large}, "
            "adjust technical potential to existing capacities"
        )
        n.generators.loc[existing_large, "p_nom_max"] = n.generators.loc[
            existing_large, "p_nom_min"
        ]

    n.generators["p_nom_max"] = n.generators["p_nom_max"].clip(lower=0)


def adjust_biomass_availability(n: pypsa.Network) -> None:
    """
    Adjust biomass generator availability to meet industrial feedstock demand.

    Scales biomass generator e_sum_max capacity proportionally if the required
    energy for industrial biomass loads exceeds the available biomass energy.

    Parameters
    ----------
    n : pypsa.Network
        Network containing biomass generators and loads

    Notes
    -----
    Modifies n.generators["e_sum_max"] in-place if adjustment is needed.
    Issues a warning log when scaling occurs.
    """
    biomass_generators = n.generators[n.generators.carrier == "solid biomass"]
    if biomass_generators.empty:
        return

    weight_sum = float(n.snapshot_weightings.generators.sum())
    industry_loads = n.loads[n.loads.carrier == "solid biomass for industry"]

    if weight_sum > 0 and not industry_loads.empty:
        required_energy = float((industry_loads.p_set * weight_sum).sum())
        available_energy = float(
            biomass_generators["e_sum_max"].replace([np.inf, -np.inf], pd.NA).sum()
        )

        if pd.isna(available_energy) or available_energy <= 0:
            available_energy = 0.0

        if required_energy > available_energy:
            deficit = required_energy - available_energy
            shares = biomass_generators["e_sum_max"].replace(0.0, pd.NA)

            if shares.notna().any():
                shares = shares.fillna(0.0) / shares.sum()
            else:
                shares = pd.Series(
                    1.0 / len(biomass_generators), index=biomass_generators.index
                )

            for gen, share in shares.items():
                n.generators.at[gen, "e_sum_max"] += deficit * share

            logger.warning(
                "Scaled solid biomass availability by %+0.1f MWh to match industrial feedstock demand.",
                deficit,
            )


def add_electricity_components(
    n: pypsa.Network,
    inputs,
    params,
    costs,
    renewable_carriers: set,
    extendable_carriers: dict,
    conventional_carriers: dict,
    foresight: str,
    sector_mode: bool,
    carriers_to_keep: dict,
) -> None:
    """
    Add all electricity-only components to the network.

    This includes conventional generators, wind/solar, hydro, and storage.
    Corresponds to the functionality from add_electricity.py.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify (in-place)
    inputs : Snakemake input object
        Input file paths from snakemake
    params : object
        Configuration parameters
    costs : pd.DataFrame
        Technology costs
    renewable_carriers : set
        Set of renewable carrier names
    extendable_carriers : dict
        Extendable carrier configuration
    conventional_carriers : dict
        Conventional carrier configuration
    foresight : str
        Foresight mode ('overnight', 'myopic', 'perfect')
    sector_mode : bool
        Whether sector coupling is enabled
    carriers_to_keep : dict
        Carriers to keep when sector_mode is enabled

    Notes
    -----
    Modifies network in-place. Adds generators, storage units, stores, and loads.
    """
    logger.info("Adding electricity components")

    # Load and aggregate powerplants
    ppl = load_and_aggregate_powerplants(
        inputs["powerplants"],
        costs,
        consider_efficiency_classes=params.clustering["consider_efficiency_classes"],
        aggregation_strategies=params.clustering["aggregation_strategies"],
        exclude_carriers=params.electricity["exclude_carriers"],
    )

    # Attach load
    attach_load(
        n,
        inputs["load"],
        inputs["busmap"],
        params.load["scaling_factor"],
    )

    # Set transmission costs
    set_transmission_costs(
        n,
        costs,
        params.lines["length_factor"],
        params.links["length_factor"],
    )

    # Load unit commitment if enabled
    unit_commitment_file = inputs.get("unit_commitment")
    if unit_commitment_file and params.conventional["unit_commitment"]:
        unit_commitment = pd.read_csv(unit_commitment_file, index_col=0)
    else:
        unit_commitment = None

    # Load dynamic fuel prices if enabled
    fuel_price_file = inputs.get("fuel_price")
    if fuel_price_file and params.conventional["dynamic_fuel_price"]:
        fuel_price = pd.read_csv(
            fuel_price_file, index_col=0, header=0, parse_dates=True
        )
        fuel_price = fuel_price.reindex(n.snapshots).ffill()
    else:
        fuel_price = None

    # Gather conventional inputs
    conventional_inputs = {
        k: v for k, v in inputs.items() if k.startswith("conventional_")
    }

    # Attach conventional generators
    attach_conventional_generators(
        n,
        costs,
        ppl,
        conventional_carriers,
        extendable_carriers,
        params.conventional,
        conventional_inputs,
        unit_commitment=unit_commitment,
        fuel_price=fuel_price,
        allowed_carriers=carriers_to_keep.get("Generator") if sector_mode else None,
    )

    # Prepare landfall lengths for renewables
    landfall_lengths = {
        tech: settings["landfall_length"]
        for tech, settings in params.renewable.items()
        if "landfall_length" in settings.keys()
    }

    # Attach wind and solar
    attach_wind_and_solar(
        n,
        costs,
        dict(inputs),
        renewable_carriers,
        extendable_carriers,
        params.lines["length_factor"],
        landfall_lengths,
    )

    if foresight in ["myopic", "perfect"]:
        apply_variable_renewable_lifetimes(n, costs)

    # Attach hydro if included
    profile_hydro = inputs.get("profile_hydro")
    hydro_capacities = inputs.get("hydro_capacities")
    if "hydro" in renewable_carriers and profile_hydro and hydro_capacities:
        hydro_params = params.renewable["hydro"].copy()
        carriers = hydro_params.pop("carriers", [])
        attach_hydro(
            n,
            costs,
            ppl,
            profile_hydro,
            hydro_capacities,
            carriers,
            **hydro_params,
        )

    # Estimate renewable capacities if enabled and overnight mode
    estimate_renewable_caps = params.electricity["estimate_renewable_capacities"]
    if estimate_renewable_caps["enable"]:
        if foresight != "overnight":
            logger.info(
                "Skipping renewable capacity estimation because they are added later "
                "in add_existing_baseyear with foresight mode 'myopic'."
            )
        else:
            tech_map = estimate_renewable_caps["technology_mapping"]
            expansion_limit = estimate_renewable_caps["expansion_limit"]
            year = estimate_renewable_caps["year"]

            if estimate_renewable_caps["from_gem"]:
                attach_GEM_renewables(n, tech_map, dict(inputs))

            estimate_renewable_capacities(
                n, year, tech_map, expansion_limit, params.countries
            )

    # Update p_nom_max
    update_p_nom_max(n)

    # Attach storage units and stores
    max_hours = params.electricity["max_hours"]
    attach_storageunits(
        n,
        costs,
        extendable_carriers,
        max_hours,
        allowed_carriers=carriers_to_keep.get("StorageUnit") if sector_mode else None,
    )
    attach_stores(
        n,
        costs,
        extendable_carriers,
        allowed_carriers=carriers_to_keep.get("Store") if sector_mode else None,
    )

    finalize_electricity_network(n)

    logger.info("Completed electricity components")


def add_sector_components(
    n: pypsa.Network,
    inputs,
    params,
    costs,
    Nyears: float,
    current_horizon: int,
    conventional_carriers: dict,
    spatial,
    foresight: str,
) -> None:
    """
    Add all sector coupling components to the network.

    This includes heat, transport, industry, and other sector-specific components.
    Corresponds to the functionality from prepare_sector_network.py.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify (in-place)
    inputs : Snakemake input object
        Input file paths from snakemake
    params : object
        Configuration parameters
    costs : pd.DataFrame
        Technology costs
    Nyears : float
        Number of years represented in the snapshot weightings
    current_horizon : int
        Current planning horizon year
    conventional_carriers : dict
        Conventional carrier configuration
    spatial : object
        Spatial configuration object from define_spatial
    foresight : str
        Foresight mode (overnight, myopic, perfect)

    Notes
    -----
    Modifies network in-place. Adds sector-specific buses, loads, generators, and links.
    Only called if sector_mode is enabled.
    """
    logger.info("Adding sector components")

    # Load additional data for sectors
    pop_layout = pd.read_csv(inputs["clustered_pop_layout"], index_col=0)
    pop_weighted_energy_totals = (
        pd.read_csv(inputs["pop_weighted_energy_totals"], index_col=0) * Nyears
    )

    # Load heat totals
    pop_weighted_heat_totals = (
        pd.read_csv(inputs["pop_weighted_heat_totals"], index_col=0) * Nyears
    )
    pop_weighted_energy_totals.update(pop_weighted_heat_totals)

    # Gas input nodes
    gas_input_nodes = pd.read_csv(inputs["gas_input_nodes_simplified"], index_col=0)

    # Load heating efficiencies
    year = int(params["energy_totals_year"])
    heating_efficiencies = pd.read_csv(
        inputs["heating_efficiencies"], index_col=[1, 0]
    ).loc[year]

    # Add EU bus
    add_eu_bus(n)

    # Add CO2 tracking
    add_co2_tracking(
        n,
        costs,
        params.sector,
        sequestration_potential_file=inputs.get("sequestration_potential"),
    )
    if foresight in ["myopic", "perfect"]:
        spatial_attrs = vars(spatial)
        for carrier in conventional_carriers:
            if carrier not in spatial_attrs:
                logger.debug(
                    "Skipping carrier bus setup for %s; no spatial configuration available",
                    carrier,
                )
                continue

            add_carrier_buses(
                n=n,
                carrier=carrier,
                costs=costs,
                spatial=spatial,
                options=params.sector,
                cf_industry=params.industry,
            )

    # Add generation
    add_generation(
        n=n,
        costs=costs,
        pop_layout=pop_layout,
        conventionals=params.sector["conventional_generation"],
        spatial=spatial,
        options=params.sector,
        cf_industry=params.industry,
    )

    # Add storage and grids
    add_storage_and_grids(
        n=n,
        costs=costs,
        pop_layout=pop_layout,
        h2_cavern_file=inputs.get("h2_cavern"),
        cavern_types=params.sector["hydrogen_underground_storage_locations"],
        clustered_gas_network_file=inputs.get("clustered_gas_network"),
        gas_input_nodes=gas_input_nodes,
        spatial=spatial,
        options=params.sector,
    )

    # Add transport if enabled
    if params.sector["transport"]["enable"]:
        add_land_transport(
            n=n,
            costs=costs,
            transport_demand_file=inputs.get("transport_demand"),
            transport_data_file=inputs.get("transport_data"),
            avail_profile_file=inputs.get("avail_profile"),
            dsm_profile_file=inputs.get("dsm_profile"),
            temp_air_total_file=inputs.get("temp_air_total"),
            cf_industry=params.industry,
            options=params.sector,
            investment_year=current_horizon,
            nodes=spatial.nodes,
        )

    # Add heating if enabled
    if params.sector["heating"]["enable"]:
        # Gather heat source profiles
        heat_source_profiles = {
            source: inputs[source]
            for source in params.get("limited_heat_sources", [])
            if source in inputs.keys()
        }

        add_heat(
            n=n,
            costs=costs,
            cop_profiles_file=inputs.get("cop_profiles"),
            direct_heat_source_utilisation_profile_file=inputs.get(
                "direct_heat_source_utilisation_profiles"
            ),
            hourly_heat_demand_total_file=inputs.get("hourly_heat_demand_total"),
            ptes_e_max_pu_file=inputs.get("ptes_e_max_pu_profiles"),
            ates_e_nom_max=inputs.get("ates_potentials"),
            ates_capex_as_fraction_of_geothermal_heat_source=params.sector[
                "district_heating"
            ]["ates"]["capex_as_fraction_of_geothermal_heat_source"],
            ates_marginal_cost_charger=params.sector["district_heating"]["ates"][
                "marginal_cost_charger"
            ],
            ates_recovery_factor=params.sector["district_heating"]["ates"][
                "recovery_factor"
            ],
            enable_ates=params.sector["district_heating"]["ates"]["enable"],
            ptes_direct_utilisation_profile=inputs.get(
                "ptes_direct_utilisation_profiles"
            ),
            district_heat_share_file=inputs.get("district_heat_share"),
            solar_thermal_total_file=inputs.get("solar_thermal_total"),
            retro_cost_file=inputs.get("retro_cost"),
            floor_area_file=inputs.get("floor_area"),
            heat_source_profile_files=heat_source_profiles,
            heat_dsm_profile_file=inputs.get("heat_dsm_profile"),
            params=params,
            pop_weighted_energy_totals=pop_weighted_energy_totals,
            heating_efficiencies=heating_efficiencies,
            pop_layout=pop_layout,
            spatial=spatial,
            options=params.sector,
            investment_year=current_horizon,
        )

    # Add biomass if enabled
    if params.sector["biomass"]["enable"]:
        add_biomass(
            n=n,
            costs=costs,
            options=params.sector,
            spatial=spatial,
            cf_industry=params.industry,
            pop_layout=pop_layout,
            biomass_potentials_file=inputs.get("biomass_potentials"),
            biomass_transport_costs_file=inputs.get("biomass_transport_costs"),
            nyears=Nyears,
        )

    # Add ammonia if enabled
    if params.sector["ammonia"]:
        add_ammonia(n, costs, pop_layout, spatial, params.industry)

    # Add methanol if enabled
    if params.sector["methanol"]:
        add_methanol(
            n, costs, options=params.sector, spatial=spatial, pop_layout=pop_layout
        )

    # Add industry if enabled
    if params.sector["industry"]["enable"]:
        add_industry(
            n=n,
            costs=costs,
            industrial_demand_file=inputs.get("industrial_demand"),
            pop_layout=pop_layout,
            pop_weighted_energy_totals=pop_weighted_energy_totals,
            options=params.sector,
            spatial=spatial,
            cf_industry=params.industry,
            investment_year=current_horizon,
        )

    # Add shipping if enabled
    if params.sector["shipping"]["enable"]:
        add_shipping(
            n=n,
            costs=costs,
            shipping_demand_file=inputs.get("shipping_demand"),
            pop_layout=pop_layout,
            pop_weighted_energy_totals=pop_weighted_energy_totals,
            options=params.sector,
            spatial=spatial,
            investment_year=current_horizon,
        )

    # Add aviation if enabled
    if params.sector["aviation"]["enable"]:
        add_aviation(
            n=n,
            costs=costs,
            pop_layout=pop_layout,
            pop_weighted_energy_totals=pop_weighted_energy_totals,
            options=params.sector,
            spatial=spatial,
        )

    # Add waste heat if heating is enabled
    if params.sector["heating"]["enable"]:
        add_waste_heat(n, costs, params.sector, params.industry)

    # Add agriculture if enabled (requires H and I)
    if params.sector["agriculture"]["enable"]:
        add_agriculture(
            n,
            costs,
            pop_layout,
            pop_weighted_energy_totals,
            current_horizon,
            params.sector,
            spatial,
        )

    # Add DAC if enabled
    if params.sector["dac"]:
        add_dac(n, costs)

    # Remove electricity transmission grid if disabled
    if not params.sector["electricity_transmission_grid"]:
        decentral(n)

    # Remove H2 network if disabled
    if not params.sector["H2_network"]:
        remove_h2_network(n)

    # Add CO2 network if enabled
    if params.sector["co2_network"]:
        add_co2_network(
            n,
            costs,
            co2_network_cost_factor=params.sector["co2_network_cost_factor"],
        )

    # Add Allam cycle gas if enabled
    if params.sector["allam_cycle_gas"]:
        add_allam_gas(n, costs, pop_layout=pop_layout, spatial=spatial)

    # Cluster heat buses if enabled
    if params.sector["cluster_heat_buses"]:
        cluster_heat_buses(n)

    # Adjust biomass availability to match industrial demand
    if params.sector["biomass"]["enable"] and params.sector["industry"]["enable"]:
        adjust_biomass_availability(n)

    logger.info("Completed sector components")


def add_existing_capacities(
    n: pypsa.Network,
    inputs,
    params,
    costs,
    renewable_carriers: set,
    baseyear: int,
) -> None:
    """
    Add existing power and heating capacities installed before baseyear.

    Corresponds to the functionality from add_existing_baseyear.py.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify (in-place)
    inputs : Snakemake input object
        Input file paths from snakemake
    params : object
        Configuration parameters
    costs : pd.DataFrame
        Technology costs
    renewable_carriers : set
        Set of renewable carrier names
    baseyear : int
        Base year for existing capacities

    Notes
    -----
    Modifies network in-place. Adds generators and links representing existing capacities.
    Only called for first horizon if existing_capacities is enabled.
    """
    logger.info("Adding existing capacities")

    existing_cfg = params.existing_capacities

    # Add power capacities installed before baseyear
    grouping_years_power = existing_cfg["grouping_years_power"]
    if grouping_years_power:
        add_power_capacities_installed_before_baseyear(
            n=n,
            costs=costs,
            grouping_years=grouping_years_power,
            baseyear=baseyear,
            powerplants_file=inputs.powerplants,
            countries=params.countries,
            capacity_threshold=existing_cfg["threshold_capacity"],
            lifetime_values=params.costs["fill_values"],
            renewable_carriers=list(renewable_carriers),
        )

    # Add heating capacities if sector coupling is enabled
    if params.sector["enabled"] and params.sector["heating"]["enable"]:
        grouping_years_heat = existing_cfg["grouping_years_heat"]
        if grouping_years_heat:
            existing_heating_file = inputs.get("existing_heating_distribution")
            cop_profiles_file = inputs.get("cop_profiles")
            heating_efficiencies_file = inputs.get("heating_efficiencies")

            existing_heat = pd.read_csv(
                existing_heating_file,
                header=[0, 1],
                index_col=0,
            )
            cop_profiles = xr.open_dataarray(cop_profiles_file)

            add_heating_capacities_installed_before_baseyear(
                n=n,
                costs=costs,
                baseyear=baseyear,
                grouping_years=grouping_years_heat,
                existing_capacities=existing_heat,
                heat_pump_cop=cop_profiles,
                heat_pump_source_types=params.heat_pump_sources,
                efficiency_file=heating_efficiencies_file,
                use_time_dependent_cop=params.sector["time_dep_hp_cop"],
                default_lifetime=existing_cfg["default_heating_lifetime"],
                energy_totals_year=int(params.energy_totals_year),
                capacity_threshold=existing_cfg["threshold_capacity"],
                use_electricity_distribution_grid=params.sector[
                    "electricity_distribution_grid"
                ],
            )

    # Add build year to new assets
    add_build_year_to_new_assets(n, baseyear)

    logger.info("Completed adding existing capacities")


def apply_temporal_aggregation(
    n: pypsa.Network,
    inputs,
    params,
) -> None:
    """
    Apply temporal aggregation to the network.

    This includes temporal resolution averaging for electricity and sector components.
    Should be called before adding existing components to ensure proper time resolution.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify (in-place)
    inputs : Snakemake input object
        Input file paths from snakemake
    params : object
        Configuration parameters

    Notes
    -----
    Modifies network in-place by updating snapshot data and weightings.
    """
    logger.info("Applying temporal aggregation")

    clustering_temporal_cfg = params.clustering_temporal
    time_resolution = clustering_temporal_cfg["resolution"]
    time_segmentation = clustering_temporal_cfg["time_segmentation"]

    # Check for conflicting configuration
    has_averaging = isinstance(
        time_resolution, str
    ) and time_resolution.lower().endswith("h")
    has_segmentation = time_segmentation["enable"] and time_segmentation["segments"]

    if has_averaging and has_segmentation:
        raise ValueError(
            "Cannot use both temporal averaging (resolution) and time "
            "segmentation simultaneously. Please configure only one method:\n"
            f"  - resolution: {time_resolution}\n"
            f"  - time_segmentation.segments: {time_segmentation['segments']}\n"
            "Set one to False/empty to use the other."
        )

    # Apply temporal resolution averaging if specified
    if has_averaging:
        logger.info(f"Applying temporal averaging: {time_resolution}")
        if params.sector["enabled"]:
            snapshot_weightings_file = inputs.get("snapshot_weightings")
            n_new = set_temporal_aggregation(
                n,
                time_resolution,
                snapshot_weightings_file,
            )
        else:
            n_new = average_every_nhours(n, time_resolution, params.drop_leap_day)
        # Replace network object content
        n.__dict__.update(n_new.__dict__)

    # Apply time segmentation if specified
    if has_segmentation:
        logger.info("Applying time segmentation")
        apply_time_segmentation(
            n,
            time_segmentation["segments"],
            snakemake.config["solving"]["solver"]["name"],
        )

    logger.info("Completed temporal aggregation")


def prepare_network_for_solving(
    n: pypsa.Network,
    inputs,
    params,
    costs,
    nyears: float,
) -> None:
    """
    Apply final network preparations before solving.

    This includes emission limits, transmission limits, and time segmentation.
    Corresponds to the functionality from prepare_network.py.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify (in-place)
    inputs : Snakemake input object
        Input file paths from snakemake
    params : object
        Configuration parameters
    costs : pd.DataFrame
        Technology costs
    nyears : float
        Number of years represented in the snapshot weightings

    Notes
    -----
    Modifies network in-place. Applies limits, constraints, and time segmentation.
    """
    logger.info("Preparing network for solving")

    electricity_cfg = params.electricity

    # Add gas limit if specified (separate from CO2 budget system)
    if electricity_cfg["gaslimit_enable"]:
        gaslimit = electricity_cfg["gaslimit"]
        add_gaslimit(n, gaslimit, nyears)

    # Set transmission limit if specified
    transmission_limit = electricity_cfg["transmission_limit"]
    if isinstance(transmission_limit, str):
        kind = transmission_limit[0]
        factor = transmission_limit[1:] or "opt"
    elif isinstance(transmission_limit, (list, tuple)) and len(transmission_limit) == 2:
        kind, factor = transmission_limit
    else:
        raise ValueError(
            "transmission_limit must be a string like 'c1.25' or a (kind, factor) pair"
        )

    # Add emission prices
    emission_prices = params.costs["emission_prices"]
    co2_price_file = inputs.get("co2_price")

    if emission_prices["co2_monthly_prices"]:
        add_dynamic_emission_prices(n, co2_price_file)
    elif emission_prices["enable"]:
        add_emission_prices(
            n,
            {"co2": emission_prices["co2"]},
            exclude_co2=False,
        )

    # Add dynamic emission prices if available
    elif co2_price_file and os.path.exists(co2_price_file):
        add_dynamic_emission_prices(n, co2_price_file)

    # Set global transmission limits
    set_transmission_limit(n, kind, factor, costs, nyears)

    # Cap individual transmission capacity for AC lines and DC links
    cap_transmission_capacity(
        n,
        line_max=params.lines["s_nom_max"],
        link_max=params.links["p_nom_max"],
        line_max_extension=params.lines["s_nom_max_extension"],
        link_max_extension=params.links["p_nom_max_extension"],
        line_max_pu=params.lines["s_max_pu"],
        link_max_pu=params.links["p_max_pu"],
    )

    logger.info("Completed network preparation")


def _is_scalar_bound(bound: object) -> bool:
    return isinstance(bound, (int, float))


def _is_mapping_bound(bound: object) -> bool:
    return isinstance(bound, Mapping)


def apply_co2_budget_constraints(
    n: pypsa.Network,
    *,
    inputs,
    params,
    nyears: float,
    foresight: str,
    horizons: list[int],
    current_horizon: int,
    apply_perfect_scalar_now: bool,
) -> None:
    co2_budget = params["co2_budget"]
    upper_cfg = co2_budget["upper"]
    lower_cfg = co2_budget["lower"]

    if upper_cfg is None:
        logger.info(
            f"CO2 budget upper constraint not specified for horizon {current_horizon}. "
            "Skipping CO2 constraint for this horizon."
        )
        return

    upper_is_scalar = _is_scalar_bound(upper_cfg)
    upper_is_mapping = _is_mapping_bound(upper_cfg)

    if not (upper_is_scalar or upper_is_mapping):
        raise TypeError(
            "co2_budget.upper must be null, a number, or a dict mapping year to value. "
            f"Received {type(upper_cfg).__name__}."
        )

    if upper_is_scalar and _is_mapping_bound(lower_cfg):
        raise ValueError(
            "Invalid co2_budget configuration: when co2_budget.upper is a scalar, "
            "co2_budget.lower must be null or a scalar (not a dict)."
        )

    baseline_1990 = None
    if co2_budget["values"] == "fraction":
        upper_raw = bound_value_for_horizon(upper_cfg, current_horizon)
        lower_raw = bound_value_for_horizon(lower_cfg, current_horizon)
        if upper_raw is not None or lower_raw is not None:
            baseline_1990 = co2_emissions_year(
                countries=params.countries,
                input_eurostat=inputs["eurostat"],
                options=params.sector,
                emissions_scope=co2_budget["emissions_scope"],
                input_co2=inputs["co2"],
                year=1990,
            )

    upper, lower = co2_budget_for_horizon(
        co2_budget,
        current_horizon=current_horizon,
        baseline_1990=baseline_1990,
    )

    if upper is None:
        logger.info(
            f"CO2 budget upper constraint not specified for horizon {current_horizon}. "
            "Skipping CO2 constraint for this horizon."
        )
        return

    if upper_is_scalar:
        if foresight == "perfect" and not apply_perfect_scalar_now:
            logger.info(
                f"Deferring scalar CO2 constraint until final perfect-foresight horizon {horizons[-1]}."
            )
            return
        add_co2limit(n, upper, lower, nyears)
        return

    # upper is mapping: apply per-horizon constraint if configured for this horizon
    if foresight == "perfect":
        add_co2limit(
            n,
            upper,
            lower,
            nyears,
            suffix=f"-{current_horizon}",
            investment_period=current_horizon,
        )
    else:
        add_co2limit(n, upper, lower, nyears)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "compose_network",
            run="test-run",
            horizon="2050",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # Extract configuration and parameters
    config = snakemake.config
    params = snakemake.params
    foresight = config["foresight"]
    current_horizon = int(snakemake.wildcards.horizon)
    sector_mode = params.sector["enabled"]
    carriers_to_keep = params.pypsa_eur if sector_mode else {}
    params.temperature_limited_stores = params.sector["district_heating"][
        "temperature_limited_stores"
    ]

    # Handle both single value and list for planning_horizons
    planning_horizons = config["planning_horizons"]
    if isinstance(planning_horizons, (int, str)):
        horizons = [int(planning_horizons)]
    else:
        horizons = [int(h) for h in planning_horizons]

    # Determine if this is the first horizon
    is_first_horizon = current_horizon == horizons[0]
    is_last_horizon = current_horizon == horizons[-1]

    # Load appropriate base network based on foresight mode and horizon
    if is_first_horizon:
        # First horizon always starts from clustered network
        n = pypsa.Network(snakemake.input.clustered)
        logger.info(f"Loading clustered network for first horizon {current_horizon}")
    else:
        # Multi-horizon case: handle based on foresight mode
        if foresight == "myopic":
            # For myopic: load previous solved network and apply brownfield constraints
            n = pypsa.Network(snakemake.input.clustered)
            adjust_renewable_profiles(
                n,
                snakemake.input,
                params,
                current_horizon,
            )
        elif foresight == "perfect":
            # For perfect foresight: compose current horizon first, then concatenate
            # at the end with previously composed multi-period network
            logger.info(
                f"Loading clustered network for perfect foresight horizon {current_horizon}"
            )
            n = pypsa.Network(snakemake.input.clustered)
        else:  # overnight
            # Should not reach here for overnight with multiple horizons
            n = pypsa.Network(snakemake.input.clustered)

    if sector_mode:
        restrict_electricity_components(n, carriers_to_keep)
        remove_non_power_buses(n)

    # Validate base network state
    validate_network_state(
        n, "base network load", expected_components={"buses": 1, "lines": 0}
    )

    # Calculate year weighting
    nyears = n.snapshot_weightings.objective.sum() / 8760.0

    # Load costs with proper horizon
    costs = load_costs(snakemake.input.tech_costs)
    # Define carrier sets
    renewable_carriers = set(params.renewable_carriers)
    extendable_carriers = params.electricity["extendable_carriers"]
    conventional_carriers = params.electricity["conventional_carriers"]

    # Use snakemake inputs directly
    inputs = snakemake.input

    # ========== ELECTRICITY COMPONENTS (from add_electricity.py) ==========
    add_electricity_components(
        n,
        inputs,
        params,
        costs,
        renewable_carriers,
        extendable_carriers,
        conventional_carriers,
        foresight,
        sector_mode,
        carriers_to_keep,
    )

    # Validate after electricity components
    counts_after_elec = validate_network_state(
        n,
        "electricity components",
        expected_components={"buses": 1, "generators": 1, "loads": 1},
    )

    # Clean up orphaned components in sector mode (generators referencing non-existent buses)
    if sector_mode:
        # Remove generators with undefined buses
        valid_buses = set(n.buses.index)
        for c in n.iterate_components(n.one_port_components | n.branch_components):
            if "bus" in c.df.columns:
                invalid = ~c.df["bus"].isin(valid_buses)
                if invalid.any():
                    logger.warning(
                        f"Removing {invalid.sum()} {c.name} with undefined buses: {c.df.index[invalid].tolist()}"
                    )
                    n.remove(c.name, c.df.index[invalid])
            if "bus0" in c.df.columns:
                invalid0 = ~c.df["bus0"].isin(valid_buses)
                invalid1 = ~c.df["bus1"].isin(valid_buses)
                invalid = invalid0 | invalid1
                if invalid.any():
                    logger.warning(
                        f"Removing {invalid.sum()} {c.name} with undefined buses: {c.df.index[invalid].tolist()}"
                    )
                    n.remove(c.name, c.df.index[invalid])

    # ========== SECTOR COMPONENTS (from prepare_sector_network.py) ==========
    if sector_mode:
        # Define spatial scope for sector components
        pop_layout = pd.read_csv(inputs["clustered_pop_layout"], index_col=0)
        spatial = define_spatial(pop_layout.index, params.sector)

        add_sector_components(
            n,
            inputs,
            params,
            costs,
            nyears,
            current_horizon,
            conventional_carriers,
            spatial,
            foresight,
        )

        # Validate after sector components
        counts_after_sector = validate_network_state(
            n, "sector components", expected_components={"buses": 2}
        )

    # ========== TEMPORAL AGGREGATION (from prepare_sector_network.py / prepare_network.py) ==========
    # Apply temporal aggregation before adding existing components
    apply_temporal_aggregation(n, inputs, params)

    # Validate after temporal aggregation
    validate_network_state(n, "temporal aggregation")

    # ========== EXISTING CAPACITIES (from add_existing_baseyear.py) ==========
    if params.existing_capacities["enabled"] and is_first_horizon:
        baseyear = params.existing_capacities["baseyear"]
        add_existing_capacities(
            n,
            inputs,
            params,
            costs,
            renewable_carriers,
            baseyear,
        )

        # Validate after existing capacities
        validate_network_state(n, "existing capacities")

    # ========== BROWNFIELD FOR MYOPIC (from add_brownfield.py) ==========
    # Apply brownfield constraints for myopic mode (non-first horizon)
    if foresight == "myopic" and not is_first_horizon:
        # Add build year to new assets (those from this planning horizon)
        add_build_year_to_new_assets(n, current_horizon)

        # Load previous solved network
        n_previous = pypsa.Network(snakemake.input.network_previous)

        # Update heat pump efficiency from previous iteration
        update_heat_pump_efficiency(n, n_previous, current_horizon)

        # Update dynamic pit storage capacity if enabled
        if params.tes and params.dynamic_ptes_capacity:
            update_dynamic_ptes_capacity(n, n_previous, current_horizon)

        logger.info(
            f"Applying brownfield constraints from horizon {horizons[horizons.index(current_horizon) - 1]}"
        )

        # Apply brownfield constraints
        add_brownfield(
            n,
            n_previous,
            current_horizon,
            h2_retrofit=params["h2_retrofit"],
            h2_retrofit_capacity_per_ch4=params["h2_retrofit_capacity_per_ch4"],
            capacity_threshold=params["capacity_threshold"],
        )
        disable_grid_expansion_if_limit_hit(n)

        # Adjust renewable capacity limits based on existing capacities
        adjust_renewable_capacity_limits(
            n, str(current_horizon), params.renewable_carriers
        )

    # ========== NETWORK PREPARATION (from prepare_network.py) ==========

    prepare_network_for_solving(
        n,
        inputs,
        params,
        costs,
        nyears,
    )

    # ========== CO2 BUDGET CONSTRAINTS ==========
    # Dict-based upper/lower are applied per horizon.
    # Scalar upper/lower are applied as a single timeless constraint:
    # - overnight/myopic: applied for each composed single-period network
    # - perfect foresight: applied only in the last horizon after concatenation
    apply_co2_budget_constraints(
        n,
        inputs=inputs,
        params=params,
        nyears=nyears,
        foresight=foresight,
        horizons=horizons,
        current_horizon=current_horizon,
        apply_perfect_scalar_now=not (foresight == "perfect"),
    )

    # Validate before finalization
    validate_network_state(n, "network preparation")

    # ========== PERFECT FORESIGHT CONCATENATION (from prepare_perfect_foresight.py) ==========
    # For perfect foresight, concatenate with previous composed network
    if foresight == "perfect":
        # Set build year for all new assets in the current horizon
        add_build_year_to_new_assets(n, current_horizon)

        # Apply per-horizon adjustments BEFORE concatenation:
        # Update heat pump efficiency for current horizon (if sector coupling enabled)
        if sector_mode:
            update_heat_pump_efficiency_for_horizon(n, current_horizon)

        # Adjust store cycling behavior for perfect foresight
        adjust_stores_for_perfect_foresight(n)

        logger.info(f"Converting horizon {current_horizon} to multi-period structure")
        n.set_investment_periods([current_horizon])

        if not is_first_horizon:
            logger.info(
                "Concatenating with previous composed network for perfect foresight"
            )
            n_previous = pypsa.Network(snakemake.input.network_previous)
            n = concatenate_network_with_previous(
                n_previous=n_previous,
                n_current=n,
                current_horizon=current_horizon,
            )

            # Apply investment period weightings after concatenation
            social_discountrate = params.costs["social_discountrate"]
            apply_investment_period_weightings(n, social_discountrate)
            logger.info(
                f"Investment period weightings applied for {len(n.investment_periods)} periods"
            )

        # Sanitize energy storages
        n.stores.e_initial_per_period = ~n.stores.e_cyclic
        n.stores.e_cyclic_per_period = n.stores.e_cyclic

    # Perfect foresight scalar CO2 constraint: add only once, after concatenation.
    if (
        foresight == "perfect"
        and is_last_horizon
        and _is_scalar_bound(params["co2_budget"]["upper"])
    ):
        nyears_total = n.snapshot_weightings.objective.sum() / 8760.0
        apply_co2_budget_constraints(
            n,
            inputs=inputs,
            params=params,
            nyears=nyears_total,
            foresight=foresight,
            horizons=horizons,
            current_horizon=current_horizon,
            apply_perfect_scalar_now=True,
        )
    # ========== FINAL CLEANUP ==========

    # Sanitize carriers
    sanitize_custom_columns(n)
    sanitize_carriers(n, config)

    # Sanitize locations if present
    if "location" in n.buses:
        sanitize_locations(n)

    if "marginal_cost" in n.generators_t:
        n.generators_t["marginal_cost"] = n.generators_t["marginal_cost"].fillna(
            n.generators.marginal_cost
        )

    # Add metadata
    n.meta = dict(config, **dict(wildcards=dict(snakemake.wildcards)))
    # Store horizon in metadata for concatenation tracking
    n.meta["horizon"] = current_horizon

    # Export composed network
    logger.info(f"Exporting composed network for horizon {current_horizon}")
    logger.info(
        f"Network investment_periods before export: {list(n.investment_periods)}"
    )
    logger.info(
        f"Network snapshots periods: {list(n.snapshots.get_level_values('period').unique()) if isinstance(n.snapshots, pd.MultiIndex) else 'Not MultiIndex'}"
    )
    try:
        n.consistency_check()
    except ValueError as e:
        if "cannot include dtype 'M' in a buffer" in str(e):
            logger.warning(
                "Skipping consistency check due to pandas/PyPSA compatibility issue with MultiIndex snapshots. "
                "This is a known issue with pandas 2.3+ and PyPSA 1.0 for multi-period networks."
            )
        else:
            raise
    n.export_to_netcdf(snakemake.output[0])
