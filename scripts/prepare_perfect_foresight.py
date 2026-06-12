# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Perfect foresight utility functions for multi-period optimization.
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import (
    PYPSA_V1,
)
from scripts.add_existing_baseyear import add_build_year_to_new_assets

logger = logging.getLogger(__name__)

logger.warning(
    "Running perfect foresight is not properly tested and may not work as expected. "
    "Use at your own risk!"
)


# helper functions ---------------------------------------------------


def extend_snapshot_multiindex(
    existing_snapshots: pd.MultiIndex,
    new_timesteps: pd.Index,
    period: int,
) -> pd.MultiIndex:
    if not isinstance(existing_snapshots, pd.MultiIndex):
        raise TypeError("existing_snapshots must be a MultiIndex")

    if list(existing_snapshots.names) != ["period", "timestep"]:
        raise ValueError(
            f"existing_snapshots must have levels ['period', 'timestep'], "
            f"but has {existing_snapshots.names}"
        )

    if isinstance(new_timesteps, pd.MultiIndex):
        raise TypeError("new_timesteps must be a single-level Index, not a MultiIndex")

    new_snapshots = pd.MultiIndex.from_product(
        [[period], new_timesteps], names=["period", "timestep"]
    )

    combined = pd.MultiIndex.from_tuples(
        list(existing_snapshots) + list(new_snapshots), names=["period", "timestep"]
    )

    return combined


def concatenate_network_with_previous(
    n_previous: pypsa.Network,
    n_current: pypsa.Network,
    current_horizon: int,
) -> pypsa.Network:
    logger.info(
        f"Concatenating network for horizon {current_horizon} with previous periods"
    )

    n = n_previous.copy()

    if n.investment_periods.empty:
        raise ValueError("Previous network must have investment periods")
    if n_current.investment_periods.empty:
        raise ValueError("Current network must have investment periods")

    combined_snapshots = list(n.snapshots) + list(n_current.snapshots)
    extended_snapshots = pd.MultiIndex.from_tuples(combined_snapshots)

    n.set_snapshots(extended_snapshots)
    n.snapshot_weightings.loc[current_horizon] = n_current.snapshot_weightings
    n.set_investment_periods(list(n.investment_periods))

    previous_horizon = n_previous.investment_periods[-1]

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
            c.dynamic[attr].loc[current_horizon] = (
                c.dynamic[attr].loc[previous_horizon].values
            )
            c.dynamic[attr].loc[current_horizon, df.columns] = df.values

            persisting = c.dynamic[attr].columns.difference(df.columns)
            twins = persisting.str.replace(
                r"-\d{4}$", f"-{current_horizon}", regex=True
            )
            has_twin = twins.isin(df.columns)
            if has_twin.any():
                c.dynamic[attr].loc[current_horizon, persisting[has_twin]] = df[
                    twins[has_twin]
                ].values

            c.dynamic[attr] = c.dynamic[attr].fillna(default)

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

    snapshot_periods = list(n.snapshots.get_level_values("period").unique())
    investment_periods_list = list(n.investment_periods)

    assert snapshot_periods == investment_periods_list, (
        "Investment periods do not match snapshot periods after concatenation"
    )

    logger.info(
        f"Successfully concatenated network: {len(n.investment_periods)} investment periods"
    )
    return n


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
        "co2 sequestered",
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
    e_initial_carriers = ["solid biomass", "biogas", "co2 stored", "co2 sequestered"]
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


def apply_phase_outs(
    n: pypsa.Network, phase_outs: list[dict], horizons: list[int]
) -> None:
    """
    Cap conventional asset lifetimes to enforce national policy phase-outs.

    For each phase-out rule, assets of the listed carriers located in the listed
    countries that would otherwise operate beyond the phase-out year have their
    lifetime shortened so they retire by that year. Assets retired before the
    first planning horizon are removed. Applies to both generators
    (electricity-only) and links (sector-coupled).
    """

    def cap_lifetimes(static: pd.DataFrame, bus: str) -> None:
        for rule in phase_outs:
            country = static[bus].str[:2]
            in_scope = static.carrier.isin(rule["carriers"]) & country.isin(
                rule["countries"]
            )
            retires_late = static.build_year + static.lifetime > rule["year"]
            assets = static.index[in_scope & retires_late]
            static.loc[assets, "lifetime"] = (
                rule["year"] - static.loc[assets, "build_year"]
            ).astype(float)

    cap_lifetimes(n.generators, "bus")
    cap_lifetimes(n.links, "bus1")

    for component, static in [("Generator", n.generators), ("Link", n.links)]:
        retired = static.index[static.build_year + static.lifetime < horizons[0]]
        n.remove(component, retired)


def main(
    n: pypsa.Network,
    n_previous: pypsa.Network | None,
    params,
    current_horizon: int,
) -> pypsa.Network:
    is_first_horizon = current_horizon == params.horizons[0]

    add_build_year_to_new_assets(n, current_horizon)

    adjust_stores_for_perfect_foresight(n)

    logger.info(f"Converting horizon {current_horizon} to multi-period structure")
    n.set_investment_periods([current_horizon])

    if not is_first_horizon:
        logger.info(
            "Concatenating with previous composed network for perfect foresight"
        )
        n = concatenate_network_with_previous(
            n_previous=n_previous,
            n_current=n,
            current_horizon=current_horizon,
        )

        social_discountrate = params.costs["social_discountrate"]
        apply_investment_period_weightings(n, social_discountrate)
        logger.info(
            f"Investment period weightings applied for {len(n.investment_periods)} periods"
        )

    return n
