# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Perfect foresight utility functions for multi-period optimization.

This module provides functions for perfect foresight optimization in the
streamlined workflow. Functions are used by compose_network.py for incremental
network composition across planning horizons.

The main script portion (concat_networks, if __name__ == "__main__") is DEPRECATED
and not used in the streamlined workflow, as networks are now composed incrementally
rather than batch-concatenated.

Functions currently used in streamlined workflow:
- get_investment_weighting(): Calculate social discount weightings
- adjust_stores_for_perfect_foresight(): Configure store cycling behavior

Functions not yet integrated (TODO):
- set_all_phase_outs(): National energy policy phase-out constraints
- update_heat_pump_efficiency(): Heat pump COP updates across all periods
- apply_time_segmentation_perfect(): Time segmentation for multi-period networks
"""

import logging

import numpy as np
import pandas as pd
import pypsa
from six import iterkeys

from scripts._helpers import (
    PYPSA_V1,
    configure_logging,
    sanitize_custom_columns,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.add_electricity import sanitize_carriers
from scripts.add_existing_baseyear import add_build_year_to_new_assets

# Allow for PyPSA versions <0.35
if PYPSA_V1:
    from pypsa.common import expand_series
else:
    from pypsa.descriptors import expand_series

logger = logging.getLogger(__name__)

logger.warning(
    "Running perfect foresight is not properly tested and may not work as expected. "
    "Use at your own risk!"
)


# helper functions ---------------------------------------------------
def get_missing(df: pd.DataFrame, n: pypsa.Network, c: str) -> pd.DataFrame:
    """
    Get missing assets in network n compared to df for component c.

    Parameters
    ----------
    df : pd.DataFrame
        Static values of pypsa components
    n : pypsa.Network
        Network to which new assets should be added
    c : str
        pypsa component.list_name (e.g. "generators")

    Returns
    -------
    pd.DataFrame
        Static values of missing assets
    """
    df_final = getattr(n, c)
    missing_i = df.index.difference(df_final.index)
    return df.loc[missing_i]


def get_investment_weighting(time_weighting: pd.Series, r: float = 0.01) -> pd.Series:
    """
    Calculate investment period objective weightings based on social discounting.

    This function computes the present value weightings for each investment period
    using a social discount rate. The weightings ensure that costs in different
    periods are properly compared in net present value terms.

    Uses closed-form formula for efficiency instead of iterative summation.

    Parameters
    ----------
    time_weighting : pd.Series
        Time weightings (duration in years) for each investment period.
        Index should match investment periods.
    r : float, default 0.01
        Social discount rate per unit (e.g., 0.02 for 2% discount rate)

    Returns
    -------
    pd.Series
        Objective weightings for each investment period that account for
        social discounting over the planning horizon

    Notes
    -----
    The weighting formula accounts for:
    - The time value of money through discounting
    - The duration of each investment period
    - The cumulative time from present to each period

    For a period starting at year `start` and ending at year `end`,
    the weighting is the sum of discount factors for all years in that period,
    calculated using a closed-form geometric series formula.

    Examples
    --------
    >>> periods = pd.Series([10, 10, 10], index=[2030, 2040, 2050])
    >>> weights = get_investment_weighting(periods, r=0.02)
    >>> # Earlier periods get higher weights due to discounting
    """
    end = time_weighting.cumsum()
    start = time_weighting.cumsum().shift().fillna(0)

    # Calculate social discount factor for each period using closed-form formula
    # This sums the discount factors from start to end of each period
    if r == 0:
        # With zero discount rate, weight equals the duration
        delta = end - start
    else:
        # Geometric series sum: sum of (1/(1+r)^t) for t from start to end-1
        # Using closed form: (1/(1+r)^start) * (1 - (1/(1+r))^duration) / (1 - 1/(1+r))
        # Simplified: ((1/(1+r)^start) - (1/(1+r)^end)) / (1 - 1/(1+r))
        # Further: ((1/(1+r)^start) - (1/(1+r)^end)) * (1+r) / r
        discount_start = 1 / ((1 + r) ** start)
        discount_end = 1 / ((1 + r) ** end)
        delta = (discount_start - discount_end) * (1 + r) / r

    return delta


def add_year_to_constraints(n: pypsa.Network, baseyear: int) -> None:
    """
    Add investment period to global constraints and rename index.

    Parameters
    ----------
    n : pypsa.Network
    baseyear : int
        year in which optimized assets are built
    """

    for c in n.iterate_components(["GlobalConstraint"]):
        c.df["investment_period"] = baseyear
        c.df.rename(index=lambda x: x + "-" + str(baseyear), inplace=True)


def hvdc_transport_model(n: pypsa.Network) -> None:
    """
    Convert AC lines to DC links for multi-decade optimisation with line
    expansion.

    Losses of DC links are assumed to be 3% per 1000km
    """

    logger.info("Convert AC lines to DC links to perform multi-decade optimisation.")

    n.add(
        "Link",
        n.lines.index,
        bus0=n.lines.bus0,
        bus1=n.lines.bus1,
        p_nom_extendable=True,
        p_nom=n.lines.s_nom,
        p_nom_min=n.lines.s_nom,
        p_min_pu=-1,
        efficiency=1 - 0.03 * n.lines.length / 1000,
        marginal_cost=0,
        carrier="DC",
        length=n.lines.length,
        capital_cost=n.lines.capital_cost,
    )

    # Remove AC lines
    logger.info("Removing AC lines")
    lines_rm = n.lines.index
    n.remove("Line", lines_rm)

    # Set efficiency of all DC links to include losses depending on length
    n.links.loc[n.links.carrier == "DC", "efficiency"] = (
        1 - 0.03 * n.links.loc[n.links.carrier == "DC", "length"] / 1000
    )


def adjust_electricity_grid(n: pypsa.Network, year: int, years: list[int]) -> None:
    """
    Adjust electricity grid for multi-decade optimization.

    Parameters
    ----------
    n : pypsa.Network
        Network to adjust
    year : int
        Year in which optimized assets are built
    years : list[int]
        List of investment periods

    Returns
    -------
    None
        Modifies network in place
    """
    n.lines["carrier"] = "AC"
    links_i = n.links[n.links.carrier == "DC"].index
    if n.lines.s_nom_extendable.any() or n.links.loc[links_i, "p_nom_extendable"].any():
        hvdc_transport_model(n)
        links_i = n.links[n.links.carrier == "DC"].index
        n.links.loc[links_i, "lifetime"] = 100
        if year != years[0]:
            n.links.loc[links_i, "p_nom_min"] = 0
            n.links.loc[links_i, "p_nom"] = 0


# --------------------------------------------------------------------
def concat_networks(
    years: list[int], network_paths: list[str], social_discountrate: float
) -> pypsa.Network:
    """
    Concat given pypsa networks and add build years.

    Parameters
    ----------
    years : list[int]
        List of years representing investment periods
    network_paths : list[str]
        List of paths to network files for each investment period

    Returns
    -------
    pypsa.Network
        Network for the whole planning horizon
    """
    n = pypsa.Network()

    # Loop over each input network file and its corresponding investment year
    for i, network_path in enumerate(network_paths):
        year = years[i]
        network = pypsa.Network(network_path)
        adjust_electricity_grid(network, year, years)
        add_build_year_to_new_assets(network, year)

        # static ----------------------------------
        for component in network.iterate_components(
            [
                "Bus",
                "Carrier",
                "Generator",
                "Link",
                "Store",
                "Load",
                "Line",
                "StorageUnit",
            ]
        ):
            df_year = component.df.copy()
            missing = get_missing(df_year, n, component.list_name)

            n.add(component.name, missing.index, **missing)

        # time variant --------------------------------------------------
        network_sns = pd.MultiIndex.from_product([[year], network.snapshots])
        snapshots = n.snapshots.drop("now", errors="ignore").union(network_sns)
        n.set_snapshots(snapshots)

        # Iterate all component types in the loaded network
        for component in network.iterate_components():
            pnl = getattr(n, component.list_name + "_t")
            for k in iterkeys(component.pnl):
                pnl_year = component.pnl[k].copy().reindex(snapshots, level=1)
                if pnl_year.empty and (not (component.name == "Load" and k == "p_set")):
                    continue
                if k not in pnl:
                    # TODO: for some reason efficiency2 isn't available, used this workaround:
                    #  initialize an empty time-series DataFrame for any missing key in pnl (e.g., 'efficiency2')
                    pnl[k] = pd.DataFrame(index=snapshots)
                if component.name == "Load":
                    static_load = network.loads.loc[network.loads.p_set != 0]
                    static_load_t = expand_series(static_load.p_set, network_sns).T
                    pnl_year = pd.concat(
                        [pnl_year.reindex(network_sns), static_load_t], axis=1
                    )
                    columns = (pnl[k].columns.union(pnl_year.columns)).unique()
                    pnl[k] = pnl[k].reindex(columns=columns)
                    pnl[k].loc[pnl_year.index, pnl_year.columns] = pnl_year

                else:
                    # For components that aren't new, we just extend
                    # time-varying data from the previous investment
                    # period.
                    if i > 0:
                        pnl[k].loc[(year,)] = pnl[k].loc[(years[i - 1],)].values

                    # Now, add time-varying data for new components.
                    cols = pnl_year.columns.difference(pnl[k].columns)
                    pnl[k] = pd.concat([pnl[k], pnl_year[cols]], axis=1)

        n.snapshot_weightings.loc[year, :] = network.snapshot_weightings.values

        # (3) global constraints
        for component in network.iterate_components(["GlobalConstraint"]):
            add_year_to_constraints(network, year)
            n.add(component.name, component.df.index, **component.df)

    # set investment periods
    n.investment_periods = n.snapshots.levels[0]
    # weighting of the investment period -> assuming last period same weighting as the period before
    time_w = n.investment_periods.to_series().diff().shift(-1).ffill()
    n.investment_period_weightings["years"] = time_w
    # set objective weightings
    objective_w = get_investment_weighting(
        n.investment_period_weightings["years"], social_discountrate
    )
    n.investment_period_weightings["objective"] = objective_w
    # all former static loads are now time-dependent -> set static = 0
    n.loads["p_set"] = 0
    n.loads_t.p_set.fillna(0, inplace=True)

    return n


def apply_investment_period_weightings(
    n: pypsa.Network, social_discountrate: float
) -> None:
    """
    Apply investment period weightings for perfect foresight optimization.

    Sets both time-based weightings (years) and objective weightings (with social
    discounting) for each investment period in a multi-period network.

    Parameters
    ----------
    n : pypsa.Network
        Multi-period network with investment periods defined
    social_discountrate : float
        Social discount rate per unit (e.g., 0.02 for 2%)

    Notes
    -----
    Modifies n.investment_period_weightings in-place by setting:
    - "years": Duration of each investment period in years
    - "objective": Social-discounted objective weighting for optimization

    The function assumes the last period has the same duration as the
    period before it (using fillna forward fill).

    Only applies to multi-invest networks with defined investment periods.

    Raises
    ------
    ValueError
        If network has no investment periods defined

    Examples
    --------
    >>> n = pypsa.Network()
    >>> n.set_snapshots(pd.MultiIndex.from_product([[2030, 2040, 2050],
    ...     pd.date_range('2030-01-01', periods=3, freq='h')]))
    >>> apply_investment_period_weightings(n, 0.02)
    >>> n.investment_period_weightings
              years  objective
    2030       10.0   9.526279
    2040       10.0   7.829213
    2050       10.0   6.436916
    """
    if n.investment_periods.empty:
        raise ValueError(
            "Cannot apply investment period weightings: network has no investment periods"
        )

    # Calculate time weightings (duration of each period in years)
    # For single period, use 5 years as conservative default
    # For multiple periods, each period lasts until next, last uses same as previous
    if len(n.investment_periods) == 1:
        time_w = pd.Series([5], index=n.investment_periods)
    else:
        time_w = n.investment_periods.to_series().diff().shift(-1).ffill()
    n.investment_period_weightings["years"] = time_w

    # Calculate objective weightings with social discounting
    objective_w = get_investment_weighting(time_w, social_discountrate)
    n.investment_period_weightings["objective"] = objective_w

    logger.info(
        f"Applied investment period weightings with social discount rate {social_discountrate:.1%}"
    )
    logger.debug(f"Time weightings (years): {time_w.to_dict()}")
    logger.debug(f"Objective weightings: {objective_w.to_dict()}")


def adjust_stores_for_perfect_foresight(n: pypsa.Network) -> None:
    """
    Adjust store cycling behavior for perfect foresight multi-period optimization.

    For perfect foresight, stores need different cycling constraints than single-period
    optimization:
    - Most stores should cycle per investment period (not over entire horizon)
    - Some stores (CO2, biomass, biogas, EV) should not cycle at all
    - Biomass, biogas, and CO2 storage stores should reset initial energy at each period

    Parameters
    ----------
    n : pypsa.Network
        Multi-period network to adjust

    Notes
    -----
    Modifies n.stores attributes in-place:
    - Sets e_cyclic_per_period=True for most cyclic stores
    - Sets e_cyclic=False for stores that cycle per period
    - Sets e_cyclic_per_period=False for non-cyclic stores (CO2, biomass, etc.)
    - Sets e_initial_per_period=True for biomass, biogas, and co2 stored stores

    This ensures proper storage behavior across investment periods in
    perfect foresight optimization.
    """
    if n.stores.empty:
        logger.debug("No stores in network, skipping store adjustments")
        return

    # Stores that should cycle within each period (not over entire horizon)
    # These are typical energy storage that should start/end each period balanced
    cyclic_i = n.stores[n.stores.e_cyclic].index
    if not cyclic_i.empty:
        n.stores.loc[cyclic_i, "e_cyclic_per_period"] = True
        n.stores.loc[cyclic_i, "e_cyclic"] = False
        logger.info(
            f"Set {len(cyclic_i)} stores to cycle per period instead of over entire horizon"
        )

    # Stores that should NOT cycle (accumulate over entire horizon)
    # - CO2: accumulates or depletes over planning horizon
    # - solid biomass/biogas: limited annual resource availability
    # - EV battery: daily charging patterns, not cyclic over long periods
    non_cyclic_carriers = [
        "co2",
        "co2 stored",
        "solid biomass",
        "biogas",
        "EV battery",
    ]

    for carrier in non_cyclic_carriers:
        store_i = n.stores[n.stores.carrier == carrier].index
        if not store_i.empty:
            n.stores.loc[store_i, "e_cyclic"] = False
            n.stores.loc[store_i, "e_cyclic_per_period"] = False
            logger.debug(f"Set {len(store_i)} '{carrier}' stores to non-cyclic")

    # Stores that should reset to initial energy at start of each period
    # Biomass and biogas: annual resource availability is renewed each year
    # CO2 stored: sequestration storage resets each period
    e_initial_carriers = ["solid biomass", "biogas", "co2 stored"]
    for carrier in e_initial_carriers:
        if carrier == "solid biomass":
            store_i = n.stores[n.stores.carrier.str.contains("solid biomass")].index
        else:
            store_i = n.stores[n.stores.carrier == carrier].index

        if not store_i.empty:
            n.stores.loc[store_i, "e_initial_per_period"] = True
            logger.info(
                f"Set {len(store_i)} '{carrier}' stores to reset initial energy per period"
            )

    logger.info("Completed store cycling adjustments for perfect foresight")


def set_phase_out(
    n: pypsa.Network, carrier: list[str], ct: str, phase_out_year: int
) -> None:
    """
    Set planned phase outs for given assets.

    Parameters
    ----------
    n : pypsa.Network
        Network to adjust
    carrier : list[str]
        List of carrier names to phase out
    ct : str
        Two-letter country code
    phase_out_year : int
        Year when phase out should be complete

    Returns
    -------
    None
        Modifies network in place
    """
    df = n.links[(n.links.carrier.isin(carrier)) & (n.links.bus1.str[:2] == ct)]
    # assets which are going to be phased out before end of their lifetime
    assets_i = df[df[["build_year", "lifetime"]].sum(axis=1) > phase_out_year].index
    build_year = n.links.loc[assets_i, "build_year"]
    # adjust lifetime
    n.links.loc[assets_i, "lifetime"] = (phase_out_year - build_year).astype(float)


def set_all_phase_outs(n: pypsa.Network, years: list[int]) -> None:
    # TODO move this to a csv or to the config
    planned = [
        (["nuclear"], "DE", 2022),
        (["nuclear"], "BE", 2025),
        (["nuclear"], "ES", 2027),
        (["coal", "lignite"], "DE", 2030),
        (["coal", "lignite"], "ES", 2027),
        (["coal", "lignite"], "FR", 2022),
        (["coal", "lignite"], "GB", 2024),
        (["coal", "lignite"], "IT", 2025),
        (["coal", "lignite"], "DK", 2030),
        (["coal", "lignite"], "FI", 2030),
        (["coal", "lignite"], "HU", 2030),
        (["coal", "lignite"], "SK", 2030),
        (["coal", "lignite"], "GR", 2030),
        (["coal", "lignite"], "IE", 2030),
        (["coal", "lignite"], "NL", 2030),
        (["coal", "lignite"], "RS", 2030),
    ]
    for carrier, ct, phase_out_year in planned:
        set_phase_out(n, carrier, ct, phase_out_year)
    # remove assets which are already phased out
    remove_i = n.links[n.links[["build_year", "lifetime"]].sum(axis=1) < years[0]].index
    n.remove("Link", remove_i)


def set_carbon_constraints(
    n: pypsa.Network, co2_budget: float | None, sector_opts: str
) -> None:
    """
    Add global constraints for carbon emissions.

    Parameters
    ----------
    n : pypsa.Network
        Network to add constraints to
    co2_budget : float or None
        CO2 budget in Gt CO2; if float, converted to t CO2
    sector_opts : str
        String of sector options separated by "-"

    Returns
    -------
    pypsa.Network
        Network with carbon constraints added
    """
    if co2_budget and isinstance(co2_budget, float):
        budget = co2_budget * 1e9  # convert to t CO2

        logger.info(f"add carbon budget of {budget}")
        n.add(
            "GlobalConstraint",
            "Budget",
            type="Co2Budget",
            carrier_attribute="co2_emissions",
            sense="<=",
            constant=budget,
            investment_period=n.investment_periods[-1],
        )

        # drop other CO2 limits
        drop_i = n.global_constraints[n.global_constraints.type == "co2_limit"].index
        n.remove("GlobalConstraint", drop_i)

        n.add(
            "GlobalConstraint",
            "carbon_neutral",
            type="co2_limit",
            carrier_attribute="co2_emissions",
            sense="<=",
            constant=0,
            investment_period=n.investment_periods[-1],
        )

    # set minimum CO2 emission constraint to avoid too fast reduction
    if "co2min" in sector_opts.split("-"):
        emissions_1990 = 4.53693
        emissions_2019 = 3.344096
        target_2030 = 0.45 * emissions_1990
        annual_reduction = (emissions_2019 - target_2030) / 11
        first_year = n.snapshots.levels[0][0]
        time_weightings = n.investment_period_weightings.loc[first_year, "years"]
        co2min = emissions_2019 - ((first_year - 2019) * annual_reduction)
        logger.info(f"add minimum emissions for {first_year} of {co2min} t CO2/a")
        n.add(
            "GlobalConstraint",
            f"Co2Min-{first_year}",
            type="Co2min",
            carrier_attribute="co2_emissions",
            sense=">=",
            investment_period=first_year,
            constant=co2min * 1e9 * time_weightings,
        )


def adjust_lvlimit(n: pypsa.Network) -> None:
    """
    Convert global constraints for single investment period to one uniform if
    all attributes stay the same.
    """
    c = "GlobalConstraint"
    cols = ["carrier_attribute", "sense", "constant", "type"]
    glc_type = "transmission_volume_expansion_limit"
    if (n.df(c)[n.df(c).type == glc_type][cols].nunique() == 1).all():
        glc = n.df(c)[n.df(c).type == glc_type][cols].iloc[[0]]
        glc.index = pd.Index(["lv_limit"])
        remove_i = n.df(c)[n.df(c).type == glc_type].index
        n.remove(c, remove_i)
        n.add(c, glc.index, **glc)


def adjust_CO2_glc(n: pypsa.Network) -> None:
    c = "GlobalConstraint"
    glc_name = "CO2Limit"
    glc_type = "primary_energy"
    mask = (n.df(c).index.str.contains(glc_name)) & (n.df(c).type == glc_type)
    n.df(c).loc[mask, "type"] = "co2_limit"


def add_H2_boilers(n: pypsa.Network) -> None:
    """
    Gas boilers can be retrofitted to run with H2.

    Add H2 boilers for heating for all existing gas boilers.
    """
    c = "Link"
    logger.info("Add H2 boilers.")
    # existing gas boilers
    mask = n.links.carrier.str.contains("gas boiler") & ~n.links.p_nom_extendable
    gas_i = n.links[mask].index
    df = n.links.loc[gas_i]
    # adjust bus 0
    df["bus0"] = df.bus1.map(n.buses.location) + " H2"
    # rename carrier and index
    df["carrier"] = df.carrier.apply(
        lambda x: x.replace("gas boiler", "retrofitted H2 boiler")
    )
    df.rename(
        index=lambda x: x.replace("gas boiler", "retrofitted H2 boiler"), inplace=True
    )
    # todo, costs for retrofitting
    df["capital_costs"] = 100
    # set existing capacity to zero
    df["p_nom"] = 0
    df["p_nom_extendable"] = True
    # add H2 boilers to network
    n.add(c, df.index, **df)


def apply_time_segmentation_perfect(
    n: pypsa.Network, segments: int, solver_name: str = "cbc"
) -> None:
    """
    Aggregate time series to segments with different lengths.

    Parameters
    ----------
    n : pypsa.Network
        Network to segment
    segments : int
        Number of segments for typical period subdivision
    solver_name : str, default "cbc"
        Name of solver to use for segmentation

    Returns
    -------
    pypsa.Network
        Network with segmented time series
    """
    try:
        import tsam.timeseriesaggregation as tsam
    except ImportError:
        raise ModuleNotFoundError(
            "Optional dependency 'tsam' not found.Install via 'pip install tsam'"
        )

    # get all time-dependent data
    columns = pd.MultiIndex.from_tuples([], names=["component", "key", "asset"])
    raw = pd.DataFrame(index=n.snapshots, columns=columns)
    for c in n.iterate_components():
        for attr, pnl in c.pnl.items():
            # exclude e_min_pu which is used for SOC of EVs in the morning
            if not pnl.empty and attr != "e_min_pu":
                df = pnl.copy()
                df.columns = pd.MultiIndex.from_product([[c.name], [attr], df.columns])
                raw = pd.concat([raw, df], axis=1)
    raw = raw.dropna(axis=1)
    sn_weightings = {}

    for year in raw.index.levels[0]:
        logger.info(f"Find representative snapshots for {year}.")
        raw_t = raw.loc[year]
        # normalise all time-dependent data
        annual_max = raw_t.max().replace(0, 1)
        raw_t = raw_t.div(annual_max, level=0)
        # get representative segments
        agg = tsam.TimeSeriesAggregation(
            raw_t,
            hoursPerPeriod=len(raw_t),
            noTypicalPeriods=1,
            noSegments=int(segments),
            segmentation=True,
            solver=solver_name,
        )
        segmented = agg.createTypicalPeriods()

        weightings = segmented.index.get_level_values("Segment Duration")
        offsets = np.insert(np.cumsum(weightings[:-1]), 0, 0)
        timesteps = [raw_t.index[0] + pd.Timedelta(f"{offset}h") for offset in offsets]
        snapshots = pd.DatetimeIndex(timesteps)
        sn_weightings[year] = pd.Series(
            weightings, index=snapshots, name="weightings", dtype="float64"
        )

    sn_weightings = pd.concat(sn_weightings)
    n.set_snapshots(sn_weightings.index)
    n.snapshot_weightings = n.snapshot_weightings.mul(sn_weightings, axis=0)


def update_heat_pump_efficiency_for_horizon(
    n: pypsa.Network, current_horizon: int
) -> None:
    """
    Update efficiency of all heat pumps to use current horizon's technology (incremental).

    For incremental perfect foresight workflow: after preparing a single horizon,
    update all heat pumps in that horizon to use current technology.

    Parameters
    ----------
    n : pypsa.Network
        Single-period network converted to multi-period [current_horizon]
    current_horizon : int
        The horizon year (e.g., 2040)

    Notes
    -----
    Modifies n.links_t["efficiency"] and n.links_t["efficiency2"] in-place.
    Called BEFORE concatenation for each horizon.
    """
    if not isinstance(n.snapshots, pd.MultiIndex):
        logger.debug(
            "Network is not multi-period, skipping heat pump efficiency update"
        )
        return

    heat_pump_idx = n.links.index[n.links.index.str.contains("heat pump")]
    if heat_pump_idx.empty:
        logger.debug("No heat pumps found in network")
        return

    if current_horizon not in n.investment_periods:
        logger.warning(f"Horizon {current_horizon} not in investment periods, skipping")
        return

    # Heat pumps with current tech (e.g., "bus AB heat pump-2040")
    current_tech_hps = heat_pump_idx.str[:-4] + str(current_horizon)

    if "efficiency" not in n.links_t or n.links_t["efficiency"].empty:
        logger.debug("No dynamic efficiency data")
        return

    try:
        correct_efficiency = n.links_t["efficiency"].loc[
            (current_horizon, slice(None)), current_tech_hps
        ]
        n.links_t["efficiency"].loc[(current_horizon, slice(None)), heat_pump_idx] = (
            correct_efficiency.values
        )
        logger.info(
            f"Updated {len(heat_pump_idx)} heat pumps to {current_horizon} technology"
        )

        # Also update efficiency2 if exists
        if "efficiency2" in n.links_t and not n.links_t["efficiency2"].empty:
            hps_with_eff2 = heat_pump_idx.intersection(n.links_t["efficiency2"].columns)
            tech_with_eff2 = current_tech_hps.intersection(
                n.links_t["efficiency2"].columns
            )
            if not hps_with_eff2.empty and not tech_with_eff2.empty:
                correct_eff2 = n.links_t["efficiency2"].loc[
                    (current_horizon, slice(None)), tech_with_eff2
                ]
                n.links_t["efficiency2"].loc[
                    (current_horizon, slice(None)), hps_with_eff2
                ] = correct_eff2.values
                logger.debug(f"Updated efficiency2 for {len(hps_with_eff2)} heat pumps")
    except KeyError as e:
        logger.warning(f"Could not update heat pump efficiency: {e}")


def update_heat_pump_efficiency(n: pypsa.Network, years: list[int]) -> None:
    """
    Update the efficiency of heat pumps for all periods (batch version).

    DEPRECATED in streamlined workflow: Use update_heat_pump_efficiency_for_horizon()
    for incremental updates instead.

    Parameters
    ----------
    n : pypsa.Network
        The concatenated network.
    years : list[int]
        List of planning horizon years.

    Notes
    -----
    This batch function is only used in the old concat_networks() workflow.
    """
    logger.warning(
        "Using batch update_heat_pump_efficiency(). "
        "Consider update_heat_pump_efficiency_for_horizon() for incremental updates."
    )

    heat_pump_idx = n.links.index[n.links.index.str.contains("heat pump")]
    for year in years:
        correct_efficiency = n.links_t["efficiency"].loc[
            (year, slice(None)), heat_pump_idx.str[:-4] + str(year)
        ]
        n.links_t["efficiency"].loc[(year, slice(None)), heat_pump_idx] = (
            correct_efficiency.values
        )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_perfect_foresight",
            opts="",
            clusters="37",
            sector_opts="1p7-4380H-T-H-B-I-A-dist1",
        )
    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    update_config_from_wildcards(snakemake.config, snakemake.wildcards)
    # parameters -----------------------------------------------------------
    years = snakemake.config["scenario"]["planning_horizons"]
    social_discountrate = snakemake.params.costs["social_discountrate"]

    logger.info(
        f"Concat networks of investment period {years} with social discount rate of {social_discountrate * 100}%"
    )

    # concat prepared networks of planning horizon to single network ------------
    network_paths = [snakemake.input.brownfield_network] + [
        snakemake.input[f"network_{year}"] for year in years[1:]
    ]
    n = concat_networks(years, network_paths, social_discountrate)

    # temporal aggregate
    solver_name = snakemake.config["solving"]["solver"]["name"]
    segments = snakemake.params.time_resolution
    if isinstance(segments, (int, float)):
        apply_time_segmentation_perfect(n, segments, solver_name=solver_name)

    # adjust global constraints lv limit if the same for all years
    adjust_lvlimit(n)
    # adjust global constraints CO2 limit
    adjust_CO2_glc(n)
    # adjust stores to multi period investment
    adjust_stores_for_perfect_foresight(n)

    # set phase outs
    set_all_phase_outs(n, years)

    # add H2 boiler
    add_H2_boilers(n)

    # set carbon constraints
    set_carbon_constraints(
        n,
        co2_budget=snakemake.config["co2_budget"],
        sector_opts=snakemake.wildcards.sector_opts,
    )

    # update meta
    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))

    # update heat pump efficiency
    update_heat_pump_efficiency(n=n, years=years)

    # export network
    sanitize_custom_columns(n)
    sanitize_carriers(n, snakemake.config)
    n.export_to_netcdf(snakemake.output[0])
