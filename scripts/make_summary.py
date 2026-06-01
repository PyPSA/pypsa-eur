# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create summary CSV files for all scenario runs including costs, capacities,
capacity factors, curtailment, energy balances, prices and other metrics.
"""

import logging
from functools import wraps
from typing import TypeAlias

import numpy as np
import pandas as pd
import pypsa

try:
    from numpy import trapezoid
except ImportError:
    from numpy import trapz as trapezoid
from pypsa import NetworkCollection

from scripts._helpers import configure_logging, set_scenario_config

idx = pd.IndexSlice
logger = logging.getLogger(__name__)

NetworkLike: TypeAlias = pypsa.Network | NetworkCollection


def _loop_over_collection(func):
    """Decorator to handle NetworkCollection by looping over individual networks."""

    @wraps(func)
    def wrapper(obj: NetworkLike) -> pd.Series | pd.DataFrame:
        if not isinstance(obj, NetworkCollection):
            return func(obj)
        results = []
        for horizon, n in zip(obj.index, obj.networks):
            result = func(n)
            if isinstance(result, pd.DataFrame):
                result = result.iloc[:, 0]
            results.append(result.to_frame(name=horizon))
        return pd.concat(results, axis=1)

    return wrapper


OUTPUTS = [
    "costs",
    "capacities",
    "energy",
    "energy_balance",
    "capacity_factors",
    "metrics",
    "curtailment",
    "prices",
    "weighted_prices",
    "market_values",
    "nodal_costs",
    "nodal_capacities",
    "nodal_energy_balance",
    "nodal_capacity_factors",
    "cumulative_costs",
]


def assign_carriers(n: pypsa.Network) -> None:
    if "carrier" not in n.lines:
        n.lines["carrier"] = "AC"


def assign_locations(n: pypsa.Network) -> None:
    for c in n.components[n.one_port_components]:
        if c.static.empty:
            continue
        c.static["location"] = c.static.bus.map(n.buses.location)

    for c in n.components[n.branch_components]:
        if c.static.empty:
            continue
        c_bus_cols = c.static.filter(regex="^bus")
        locs = c_bus_cols.apply(lambda c: c.map(n.buses.location)).sort_index(axis=1)
        # Use first location that is not "EU"; take "EU" if nothing else available
        c.static["location"] = locs.apply(
            lambda row: next(
                (loc for loc in row.dropna() if loc != "EU"),
                "EU",
            ),
            axis=1,
        )


@_loop_over_collection
def calculate_nodal_capacity_factors(n: pypsa.Network) -> pd.Series:
    """
    Calculate the regional dispatched capacity factors / utilisation rates for each technology carrier based on location bus attribute.
    """
    comps = n.one_port_components ^ {"Store"} | n.passive_branch_components
    return n.statistics.capacity_factor(comps=comps, groupby=["location", "carrier"])


@_loop_over_collection
def calculate_capacity_factors(n: pypsa.Network) -> pd.Series:
    """
    Calculate the average dispatched capacity factors / utilisation rates for each technology carrier.
    """
    comps = n.one_port_components ^ {"Store"} | n.passive_branch_components
    return n.statistics.capacity_factor(comps=comps).sort_index()


@_loop_over_collection
def calculate_nodal_costs(n: pypsa.Network) -> pd.Series:
    """
    Calculate optimized regional costs for each technology split by marginal and capital costs and based on location bus attribute.
    """
    grouper = ["location", "carrier"]
    costs = pd.concat(
        {
            "capital": n.statistics.capex(groupby=grouper),
            "marginal": n.statistics.opex(groupby=grouper),
        }
    )
    costs.index.names = ["cost", "component", "location", "carrier"]
    return costs


@_loop_over_collection
def calculate_costs(n: pypsa.Network) -> pd.Series:
    """
    Calculate optimized total costs for each technology split by marginal and capital costs.
    """
    costs = pd.concat(
        {
            "capital": n.statistics.capex(),
            "marginal": n.statistics.opex(),
        }
    )
    costs.index.names = ["cost", "component", "carrier"]
    return costs


def calculate_nodal_capacities(obj: NetworkLike) -> pd.Series | pd.DataFrame:
    """
    Calculate optimized regional capacities for each technology relative to bus/bus0 based on location bus attribute.
    """
    result = obj.statistics.optimal_capacity(groupby=["location", "carrier"])
    if isinstance(obj, NetworkCollection):
        result = result.unstack(level="horizon")
        result.columns.name = None
    return result


def calculate_capacities(obj: NetworkLike) -> pd.Series | pd.DataFrame:
    """
    Calculate optimized total capacities for each technology relative to bus/bus0.
    """
    result = obj.statistics.optimal_capacity()
    if isinstance(obj, NetworkCollection):
        result = result.unstack(level="horizon")
        result.columns.name = None
    return result


@_loop_over_collection
def calculate_curtailment(n: pypsa.Network) -> pd.Series:
    """
    Calculate the curtailment of electricity generation technologies in percent.
    """
    carriers = ["AC", "low voltage"]
    duration = n.snapshot_weightings.generators.sum()
    curtailed_abs = n.statistics.curtailment(
        bus_carrier=carriers, aggregate_across_components=True
    )
    available = (
        n.statistics.optimal_capacity("Generator", bus_carrier=carriers) * duration
    )
    curtailed_rel = curtailed_abs / available * 100
    return curtailed_rel.sort_index()


@_loop_over_collection
def calculate_energy(n: pypsa.Network) -> pd.Series:
    """
    Calculate the net energy supply (positive) and consumption (negative) by technology carrier across all ports.
    """
    energy = n.statistics.energy_balance(groupby="carrier")
    return energy.sort_values(ascending=False) if isinstance(energy, pd.Series) else energy


@_loop_over_collection
def calculate_energy_balance(n: pypsa.Network) -> pd.Series:
    """
    Calculate the energy supply (positive) and consumption (negative) by technology carrier for each bus carrier.
    """
    balance = n.statistics.energy_balance()
    return (
        balance.sort_values(ascending=False)
        if isinstance(balance, pd.Series)
        else balance
    )


@_loop_over_collection
def calculate_nodal_energy_balance(n: pypsa.Network) -> pd.Series:
    """
    Calculate the regional energy balances (positive values for supply, negative values for consumption) for each technology carrier and bus carrier based on the location bus attribute.
    """
    return n.statistics.energy_balance(groupby=["carrier", "location", "bus_carrier"])


@_loop_over_collection
def calculate_metrics(n: pypsa.Network) -> pd.Series:
    """
    Calculate system-level metrics, e.g. shadow prices, grid expansion, total costs.
    """
    capex = n.statistics.capex()
    opex = n.statistics.opex()
    buses_i = n.buses.query("carrier == 'AC'").index

    def metrics_for(period: int | None) -> pd.Series:
        metrics = {}
        dc_links = n.links.query("carrier == 'DC'")
        metrics["line_volume_DC"] = dc_links.eval("length * p_nom_opt").sum()
        metrics["line_volume_AC"] = n.lines.eval("length * s_nom_opt").sum()
        metrics["line_volume"] = metrics["line_volume_AC"] + metrics["line_volume_DC"]

        total_capex = capex[period].sum() if period is not None else capex.sum()
        total_opex = opex[period].sum() if period is not None else opex.sum()
        metrics["total costs"] = total_capex + total_opex

        prices = n.buses_t.marginal_price[buses_i]
        if period is not None:
            prices = prices.loc[period]
        metrics["electricity_price_zero_hours"] = (
            prices.where(prices < 0.1).count().sum() / prices.size
        )
        flat_prices = prices.stack()
        metrics["electricity_price_mean"] = flat_prices.mean()
        metrics["electricity_price_std"] = flat_prices.std()

        if "lv_limit" in n.global_constraints.index:
            metrics["line_volume_limit"] = n.global_constraints.at[
                "lv_limit", "constant"
            ]
            metrics["line_volume_shadow"] = n.global_constraints.at["lv_limit", "mu"]
        if "CO2Limit" in n.global_constraints.index:
            metrics["co2_shadow"] = n.global_constraints.at["CO2Limit", "mu"]
        if "co2_sequestration_limit" in n.global_constraints.index:
            metrics["co2_storage_shadow"] = n.global_constraints.at[
                "co2_sequestration_limit", "mu"
            ]
        return pd.Series(metrics)

    if n.has_investment_periods:
        columns = {p: metrics_for(p) for p in n.investment_periods}
        return pd.DataFrame(columns).sort_index()
    return metrics_for(None).sort_index()


@_loop_over_collection
def calculate_prices(n: pypsa.Network) -> pd.Series:
    """
    Calculate time-averaged prices per carrier.
    """
    marginal_price = n.buses_t.marginal_price
    if n.has_investment_periods:
        by_bus = marginal_price.groupby(level=0).mean().T
        return by_bus.groupby(n.buses.carrier).mean().sort_index()
    return marginal_price.mean().groupby(n.buses.carrier).mean().sort_index()


@_loop_over_collection
def calculate_weighted_prices(n: pypsa.Network) -> pd.Series:
    """
    Calculate load-weighted prices per bus carrier.
    """

    def weighted_for(period: int | None) -> pd.Series:
        weighted_prices = {}
        for carrier in n.buses.carrier.unique():
            load = n.statistics.withdrawal(
                groupby="bus",
                aggregate_time=False,
                bus_carrier=carrier,
                aggregate_across_components=True,
            ).T
            if load.empty:
                continue
            price = n.buses_t.marginal_price.loc[:, n.buses.carrier == carrier]
            weights = n.snapshot_weightings.generators
            if period is not None:
                load = load.loc[period]
                price = price.loc[period]
                weights = weights.loc[period]
            if load.sum().sum() <= 0:
                continue
            price = price.reindex(columns=load.columns, fill_value=1)
            a = weights @ (load * price).sum(axis=1)
            b = weights @ load.sum(axis=1)
            weighted_prices[carrier] = a / b
        return pd.Series(weighted_prices)

    if n.has_investment_periods:
        columns = {p: weighted_for(p) for p in n.investment_periods}
        return pd.DataFrame(columns).sort_index()
    return weighted_for(None).sort_index()


@_loop_over_collection
def calculate_market_values(n: pypsa.Network) -> pd.Series:
    """
    Calculate market values for electricity.
    """
    market_value = n.statistics.market_value(
        bus_carrier="AC", aggregate_across_components=True
    )
    if isinstance(market_value, pd.Series):
        return market_value.sort_values().dropna()
    return market_value.dropna(how="all")


def calculate_cumulative_costs(
    costs_df: pd.DataFrame, planning_horizons: pd.Index
) -> pd.DataFrame:
    """
    Calculate cumulative costs with social discounting for myopic/perfect foresight.

    Only applicable when multiple horizons exist. Applies social discount rates
    from 0% to 10% and integrates costs over the transition path.
    """
    if len(planning_horizons) <= 1:
        return pd.DataFrame()

    total_costs = costs_df.sum()
    discount_rates = pd.Index(
        data=np.arange(0, 0.11, 0.01), name="social_discount_rate"
    )

    cumulative_cost = pd.DataFrame(
        index=discount_rates, columns=planning_horizons, dtype=float
    )

    for rate in discount_rates:
        for horizon in planning_horizons:
            years_diff = horizon - planning_horizons[0]
            cumulative_cost.loc[rate, horizon] = total_costs[horizon] / (
                (1 + rate) ** years_diff
            )

    integrated_costs = pd.Series(
        index=discount_rates, name="cumulative_cost", dtype=float
    )
    for rate in discount_rates:
        integrated_costs[rate] = trapezoid(
            cumulative_cost.loc[rate, :].values, x=planning_horizons.values
        )

    cumulative_cost["cumulative_cost"] = integrated_costs
    return cumulative_cost


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "make_summary",
            configfiles="config/test/config.overnight.yaml",
            horizon=2030,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    pypsa.options.params.statistics.nice_names = False
    pypsa.options.params.statistics.drop_zero = False

    foresight = snakemake.config["foresight"]
    planning_horizons = pd.Index(snakemake.config["planning_horizons"])
    network_files = snakemake.input.networks

    logger.debug(f"Processing {foresight} mode with {len(network_files)} network(s)")

    is_overnight = foresight == "overnight" or len(network_files) == 1

    if len(network_files) > 1 and foresight != "perfect":
        logger.debug(f"Loading {len(network_files)} networks for myopic mode")
        networks = []
        for i, network_file in enumerate(network_files):
            logger.debug(
                f"Loading network {i + 1}/{len(network_files)}: {network_file}"
            )
            n = pypsa.Network(network_file)
            assign_carriers(n)
            assign_locations(n)
            networks.append(n)
        obj: NetworkLike = NetworkCollection(
            networks, index=pd.Index(planning_horizons, name="horizon")
        )
        logger.debug(
            f"Created NetworkCollection with horizons: {list(planning_horizons)}"
        )
    else:
        logger.debug(f"Loading single network (foresight: {foresight})")
        obj = pypsa.Network(network_files[0])
        assign_carriers(obj)
        assign_locations(obj)

    costs_df = None
    for output in OUTPUTS:
        if output == "cumulative_costs":
            continue
        logger.debug(f"Calculating {output}")
        result = globals()["calculate_" + output](obj)
        if is_overnight and isinstance(result, pd.Series):
            result = pd.DataFrame({planning_horizons[0]: result})
        if output == "costs":
            costs_df = result
        result.to_csv(snakemake.output[output])

    if costs_df is not None and len(planning_horizons) > 1:
        logger.debug("Calculating cumulative_costs")
        cumulative_costs = calculate_cumulative_costs(costs_df, planning_horizons)
        cumulative_costs.to_csv(snakemake.output.cumulative_costs)
    else:
        pd.DataFrame().to_csv(snakemake.output.cumulative_costs)

    logger.debug("Summary calculation completed successfully")
