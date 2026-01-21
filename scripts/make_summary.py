# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create summary CSV files for all scenario runs including costs, capacities,
capacity factors, curtailment, energy balances, prices and other metrics.
"""

import logging

import numpy as np
import pandas as pd
import pypsa
from numpy import atleast_1d

try:
    from numpy import trapezoid
except ImportError:
    # before numpy 2.0
    from numpy import trapz as trapezoid
from pypsa import NetworkCollection

from scripts._helpers import configure_logging, set_scenario_config

idx = pd.IndexSlice
logger = logging.getLogger(__name__)

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
    for c in n.iterate_components(n.one_port_components):
        c.df["location"] = c.df.bus.map(n.buses.location)

    for c in n.iterate_components(n.branch_components):
        c_bus_cols = c.df.filter(regex="^bus")
        locs = c_bus_cols.apply(lambda c: c.map(n.buses.location)).sort_index(axis=1)
        # Use first location that is not "EU"; take "EU" if nothing else available
        c.df["location"] = locs.apply(
            lambda row: next(
                (loc for loc in row.dropna() if loc != "EU"),
                "EU",
            ),
            axis=1,
        )


def calculate_nodal_capacity_factors(n: pypsa.Network) -> pd.Series:
    """
    Calculate the regional dispatched capacity factors / utilisation rates for each technology carrier based on location bus attribute.
    """
    comps = n.one_port_components ^ {"Store"} | n.passive_branch_components
    return n.statistics.capacity_factor(comps=comps, groupby=["location", "carrier"])


def calculate_capacity_factors(n: pypsa.Network) -> pd.Series:
    """
    Calculate the average dispatched capacity factors / utilisation rates for each technology carrier.

    Returns
    -------
    pd.Series
        MultiIndex Series with levels ["component", "carrier"]
    """

    comps = n.one_port_components ^ {"Store"} | n.passive_branch_components
    return n.statistics.capacity_factor(comps=comps).sort_index()


def calculate_nodal_costs(n: pypsa.Network) -> pd.Series:
    """
    Calculate optimized regional costs for each technology split by marginal and capital costs and based on location bus attribute.

    Returns
    -------
    pd.Series
        MultiIndex Series with levels ["cost", "component", "location", "carrier"]
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


def calculate_costs(n: pypsa.Network) -> pd.Series:
    """
    Calculate optimized total costs for each technology split by marginal and capital costs.

    Returns
    -------
    pd.Series
        MultiIndex Series with levels ["cost", "component", "carrier"]
    """
    costs = pd.concat(
        {
            "capital": n.statistics.capex(),
            "marginal": n.statistics.opex(),
        }
    )
    costs.index.names = ["cost", "component", "carrier"]

    return costs


def calculate_nodal_capacities(n: pypsa.Network) -> pd.Series:
    """
    Calculate optimized regional capacities for each technology relative to bus/bus0 based on location bus attribute.

    Returns
    -------
    pd.Series
        MultiIndex Series with levels ["component", "location", "carrier"]
    """
    return n.statistics.optimal_capacity(groupby=["location", "carrier"])


def calculate_capacities(n: pypsa.Network) -> pd.Series:
    """
    Calculate optimized total capacities for each technology relative to bus/bus0.

    Returns
    -------
    pd.Series
        MultiIndex Series with levels ["component", "carrier"]
    """
    return n.statistics.optimal_capacity()


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


def calculate_energy(n: pypsa.Network) -> pd.Series:
    """
    Calculate the net energy supply (positive) and consumption (negative) by technology carrier across all ports.

    Returns
    -------
    pd.Series or pd.DataFrame
        Single-period: Series with MultiIndex ["component", "carrier"]
        Multi-period: DataFrame with periods as columns
    """
    result = n.statistics.energy_balance(groupby="carrier")
    if isinstance(result, pd.DataFrame):
        # Multi-period network returns DataFrame with periods as columns - keep it!
        return result
    return result.sort_values(ascending=False)


def calculate_energy_balance(n: pypsa.Network) -> pd.Series:
    """
    Calculate the energy supply (positive) and consumption (negative) by technology carrier for each bus carrier.

    Returns
    -------
    pd.Series or pd.DataFrame
        Single-period: Series with MultiIndex ["component", "carrier", "bus_carrier"]
        Multi-period: DataFrame with periods as columns

    Examples
    --------
    >>> eb = calculate_energy_balance(n)
    >>> eb.xs("methanol", level='bus_carrier')
    """
    result = n.statistics.energy_balance()
    if isinstance(result, pd.DataFrame):
        # Multi-period network returns DataFrame with periods as columns - keep it!
        return result
    return result.sort_values(ascending=False)


def calculate_nodal_energy_balance(n: pypsa.Network) -> pd.Series:
    """
    Calculate the regional energy balances (positive values for supply, negative values for consumption) for each technology carrier and bus carrier based on the location bus attribute.

    Returns
    -------
    pd.Series
        MultiIndex Series with levels ["component", "carrier", "location", "bus_carrier"]

    Examples
    --------
    >>> eb = calculate_nodal_energy_balance(n)
    >>> eb.xs(("AC", "BE0 0"), level=["bus_carrier", "location"])
    """
    return n.statistics.energy_balance(groupby=["carrier", "location", "bus_carrier"])


def calculate_metrics(n: pypsa.Network) -> pd.Series:
    """
    Calculate system-level metrics, e.g. shadow prices, grid expansion, total costs.
    Also calculate average, standard deviation and share of zero hours for electricity prices.
    """

    metrics = {}

    dc_links = n.links.query("carrier == 'DC'")
    metrics["line_volume_DC"] = dc_links.eval("length * p_nom_opt").sum()
    metrics["line_volume_AC"] = n.lines.eval("length * s_nom_opt").sum()
    metrics["line_volume"] = metrics["line_volume_AC"] + metrics["line_volume_DC"]

    metrics["total costs"] = n.statistics.capex().sum() + n.statistics.opex().sum()

    buses_i = n.buses.query("carrier == 'AC'").index
    prices = n.buses_t.marginal_price[buses_i]

    # threshold higher than marginal_cost of VRE
    zero_hours = prices.where(prices < 0.1).count().sum()
    metrics["electricity_price_zero_hours"] = zero_hours / prices.size
    metrics["electricity_price_mean"] = prices.unstack().mean()
    metrics["electricity_price_std"] = prices.unstack().std()

    if "lv_limit" in n.global_constraints.index:
        metrics["line_volume_limit"] = n.global_constraints.at["lv_limit", "constant"]
        metrics["line_volume_shadow"] = n.global_constraints.at["lv_limit", "mu"]

    if "CO2Limit" in n.global_constraints.index:
        metrics["co2_shadow"] = n.global_constraints.at["CO2Limit", "mu"]

    if "co2_sequestration_limit" in n.global_constraints.index:
        metrics["co2_storage_shadow"] = n.global_constraints.at[
            "co2_sequestration_limit", "mu"
        ]

    return pd.Series(metrics).sort_index()


def calculate_prices(n: pypsa.Network) -> pd.Series:
    """
    Calculate time-averaged prices per carrier.
    """
    return n.buses_t.marginal_price.mean().groupby(n.buses.carrier).mean().sort_index()


def calculate_weighted_prices(n: pypsa.Network) -> pd.Series:
    """
    Calculate load-weighted prices per bus carrier.
    """
    carriers = n.buses.carrier.unique()

    weighted_prices = {}

    for carrier in carriers:
        load = n.statistics.withdrawal(
            groupby="bus",
            aggregate_time=False,
            bus_carrier=carrier,
            aggregate_across_components=True,
        ).T

        if not load.empty and load.sum().sum() > 0:
            price = n.buses_t.marginal_price.loc[:, n.buses.carrier == carrier]
            price = price.reindex(columns=load.columns, fill_value=1)

            weights = n.snapshot_weightings.generators
            a = weights @ (load * price).sum(axis=1)
            b = weights @ load.sum(axis=1)
            weighted_prices[carrier] = a / b

    return pd.Series(weighted_prices).sort_index()


def calculate_market_values(n: pypsa.Network) -> pd.Series:
    """
    Calculate market values for electricity.

    Returns
    -------
    pd.Series or pd.DataFrame
        Single-period: Series
        Multi-period: DataFrame with periods as columns
    """
    result = n.statistics.market_value(
        bus_carrier="AC", aggregate_across_components=True
    )
    if isinstance(result, pd.DataFrame):
        # Multi-period network returns DataFrame with periods as columns - keep it!
        return result.dropna()
    return result.sort_values().dropna()


def calculate_cumulative_costs(
    costs_df: pd.DataFrame, planning_horizons: pd.Index
) -> pd.DataFrame:
    """
    Calculate cumulative costs with social discounting for myopic/perfect foresight.

    Only applicable when multiple horizons exist. Applies social discount rates
    from 0% to 10% and integrates costs over the transition path.

    Parameters
    ----------
    costs_df : pd.DataFrame
        Costs DataFrame with index (cost, component, carrier) and horizon columns
    planning_horizons : pd.Index
        Planning horizons (e.g., [2030, 2040, 2050])

    Returns
    -------
    pd.DataFrame
        Cumulative costs for different social discount rates with horizons as columns
        and an additional 'cumulative_cost' column with integrated values
    """
    if len(planning_horizons) <= 1:
        # Not applicable for single horizon (overnight)
        return pd.DataFrame()

    # Sum costs across all components/carriers for each horizon
    total_costs = costs_df.sum()

    # Create discount rate index (0% to 10% in 1% increments)
    discount_rates = pd.Index(
        data=np.arange(0, 0.11, 0.01), name="social_discount_rate"
    )

    cumulative_cost = pd.DataFrame(
        index=discount_rates, columns=planning_horizons, dtype=float
    )

    # Apply social discounting: express costs in money value of first planning horizon
    for rate in discount_rates:
        for horizon in planning_horizons:
            years_diff = horizon - planning_horizons[0]
            cumulative_cost.loc[rate, horizon] = total_costs[horizon] / (
                (1 + rate) ** years_diff
            )

    # Integrate costs over transition path using trapezoidal rule
    integrated_costs = pd.Series(
        index=discount_rates, name="cumulative_cost", dtype=float
    )
    for rate in discount_rates:
        integrated_costs[rate] = trapezoid(
            cumulative_cost.loc[rate, :].values, x=planning_horizons.values
        )

    # Add cumulative cost as additional column
    cumulative_cost["cumulative_cost"] = integrated_costs

    return cumulative_cost


def calculate_nodal_capacity_factors_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate the regional dispatched capacity factors for each technology carrier based on location bus attribute using NetworkCollection.
    """
    comps = (
        nc.networks.iloc[0].one_port_components ^ {"Store"}
        | nc.networks.iloc[0].passive_branch_components
    )
    result = nc.statistics.capacity_factor(comps=comps, groupby=["location", "carrier"])
    result = result.unstack(level="horizon")
    result.columns.name = None
    return result


def calculate_capacity_factors_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate the average dispatched capacity factors for each technology carrier using NetworkCollection.
    """
    comps = (
        nc.networks.iloc[0].one_port_components ^ {"Store"}
        | nc.networks.iloc[0].passive_branch_components
    )
    result = nc.statistics.capacity_factor(comps=comps).sort_index()
    result = result.unstack(level="horizon")
    result.columns.name = None
    return result


def calculate_nodal_costs_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate optimized regional costs for each technology split by marginal and capital costs using NetworkCollection.
    """
    grouper = ["location", "carrier"]
    costs = pd.concat(
        {
            "capital": nc.statistics.capex(groupby=grouper),
            "marginal": nc.statistics.opex(groupby=grouper),
        }
    )
    costs.index.names = ["cost", "component", "horizon", "location", "carrier"]
    costs = costs.unstack(level="horizon")
    costs.columns.name = None
    return costs


def calculate_costs_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate optimized total costs for each technology split by marginal and capital costs using NetworkCollection.
    """
    costs = pd.concat(
        {
            "capital": nc.statistics.capex(),
            "marginal": nc.statistics.opex(),
        }
    )
    costs.index.names = ["cost", "component", "horizon", "carrier"]
    costs = costs.unstack(level="horizon")
    costs.columns.name = None
    return costs


def calculate_nodal_capacities_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate optimized regional capacities for each technology using NetworkCollection.
    """
    result = nc.statistics.optimal_capacity(groupby=["location", "carrier"])
    result = result.unstack(level="horizon")
    result.columns.name = None
    return result


def calculate_capacities_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate optimized total capacities for each technology using NetworkCollection.
    """
    result = nc.statistics.optimal_capacity()
    result = result.unstack(level="horizon")
    result.columns.name = None
    return result


def calculate_curtailment_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate the curtailment of electricity generation technologies in percent using NetworkCollection.
    """
    carriers = ["AC", "low voltage"]

    # Calculate for each network and aggregate
    curtailment_series = []
    for horizon, n in zip(nc.index, nc.networks):
        duration = n.snapshot_weightings.generators.sum()

        curtailed_abs = n.statistics.curtailment(
            bus_carrier=carriers, aggregate_across_components=True
        )
        available = (
            n.statistics.optimal_capacity("Generator", bus_carrier=carriers) * duration
        )

        curtailed_rel = curtailed_abs / available * 100
        curtailed_rel.name = horizon
        curtailment_series.append(curtailed_rel)

    result = pd.concat(curtailment_series, axis=1).sort_index()
    return result


def calculate_energy_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate the net energy supply and consumption by technology carrier using NetworkCollection.
    """
    result = nc.statistics.energy_balance(groupby="carrier").sort_values(
        ascending=False
    )
    result = result.unstack(level="horizon")
    result.columns.name = None
    return result


def calculate_energy_balance_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate the energy supply and consumption by technology carrier for each bus carrier using NetworkCollection.
    """
    result = nc.statistics.energy_balance().sort_values(ascending=False)
    result = result.unstack(level="horizon")
    result.columns.name = None
    return result


def calculate_nodal_energy_balance_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate the regional energy balances for each technology carrier and bus carrier using NetworkCollection.
    """
    result = nc.statistics.energy_balance(
        groupby=["carrier", "location", "bus_carrier"]
    )
    result = result.unstack(level="horizon")
    result.columns.name = None
    return result


def calculate_metrics_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate system-level metrics for each horizon using NetworkCollection.
    """
    metrics_list = []

    for horizon, n in zip(nc.index, nc.networks):
        metrics = {}

        dc_links = n.links.query("carrier == 'DC'")
        metrics["line_volume_DC"] = dc_links.eval("length * p_nom_opt").sum()
        metrics["line_volume_AC"] = n.lines.eval("length * s_nom_opt").sum()
        metrics["line_volume"] = metrics["line_volume_AC"] + metrics["line_volume_DC"]

        metrics["total costs"] = n.statistics.capex().sum() + n.statistics.opex().sum()

        buses_i = n.buses.query("carrier == 'AC'").index
        prices = n.buses_t.marginal_price[buses_i]

        zero_hours = prices.where(prices < 0.1).count().sum()
        metrics["electricity_price_zero_hours"] = zero_hours / prices.size
        metrics["electricity_price_mean"] = prices.unstack().mean()
        metrics["electricity_price_std"] = prices.unstack().std()

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

        metrics_series = pd.Series(metrics).sort_index()
        metrics_series.name = horizon
        metrics_list.append(metrics_series)

    result = pd.concat(metrics_list, axis=1)
    return result


def calculate_prices_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate time-averaged prices per carrier using NetworkCollection.
    """
    prices_list = []

    for horizon, n in zip(nc.index, nc.networks):
        prices = (
            n.buses_t.marginal_price.mean().groupby(n.buses.carrier).mean().sort_index()
        )
        prices.name = horizon
        prices_list.append(prices)

    result = pd.concat(prices_list, axis=1)
    return result


def calculate_weighted_prices_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate load-weighted prices per bus carrier using NetworkCollection.
    """
    weighted_prices_list = []

    for horizon, n in zip(nc.index, nc.networks):
        carriers = n.buses.carrier.unique()
        weighted_prices = {}

        for carrier in carriers:
            load = n.statistics.withdrawal(
                groupby="bus",
                aggregate_time=False,
                bus_carrier=carrier,
                aggregate_across_components=True,
            ).T

            if not load.empty and load.sum().sum() > 0:
                price = n.buses_t.marginal_price.loc[:, n.buses.carrier == carrier]
                price = price.reindex(columns=load.columns, fill_value=1)

                weights = n.snapshot_weightings.generators
                a = weights @ (load * price).sum(axis=1)
                b = weights @ load.sum(axis=1)
                weighted_prices[carrier] = a / b

        wp_series = pd.Series(weighted_prices).sort_index()
        wp_series.name = horizon
        weighted_prices_list.append(wp_series)

    result = pd.concat(weighted_prices_list, axis=1)
    return result


def calculate_market_values_collection(nc: NetworkCollection) -> pd.Series:
    """
    Calculate market values for electricity using NetworkCollection.
    """
    market_values_list = []

    for horizon, n in zip(nc.index, nc.networks):
        mv = n.statistics.market_value(
            bus_carrier="AC", aggregate_across_components=True
        )
        # Should be Series for single-period networks in myopic mode
        mv = mv.sort_values().dropna()
        mv.name = horizon
        market_values_list.append(mv)

    result = pd.concat(market_values_list, axis=1)
    return result


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
    planning_horizons = pd.Index(atleast_1d(snakemake.config["planning_horizons"]))
    network_files = snakemake.input.networks

    logger.debug(f"Processing {foresight} mode with {len(network_files)} network(s)")

    if foresight == "perfect":
        # Perfect foresight: Single multi-period network
        logger.debug("Loading multi-period network for perfect foresight")
        n = pypsa.Network(network_files[0])
        assign_carriers(n)
        assign_locations(n)

        costs_df = None
        for output in OUTPUTS:
            if output == "cumulative_costs":
                continue  # Handle separately after costs are calculated
            logger.debug(f"Calculating {output}")
            result = globals()["calculate_" + output](n)
            if output == "costs":
                costs_df = result  # Save for cumulative costs calculation
            result.to_csv(snakemake.output[output])

        # Calculate cumulative costs from costs DataFrame
        if costs_df is not None and len(planning_horizons) > 1:
            logger.debug("Calculating cumulative_costs")
            cumulative_costs = calculate_cumulative_costs(costs_df, planning_horizons)
            cumulative_costs.to_csv(snakemake.output.cumulative_costs)

    elif len(network_files) == 1:
        # Overnight mode: Single network, single horizon
        logger.debug(
            f"Loading single network for overnight mode (horizon: {planning_horizons[0]})"
        )
        n = pypsa.Network(network_files[0])
        assign_carriers(n)
        assign_locations(n)

        costs_df = None
        for output in OUTPUTS:
            if output == "cumulative_costs":
                continue  # Not applicable for single horizon
            logger.debug(f"Calculating {output}")
            result = globals()["calculate_" + output](n)
            if output == "costs":
                costs_df = result  # Save for cumulative costs calculation
            # Wrap in DataFrame with horizon column
            if isinstance(result, pd.Series):
                result = pd.DataFrame({planning_horizons[0]: result})
            result.to_csv(snakemake.output[output])

        # Cumulative costs not applicable for overnight (single horizon)
        pd.DataFrame().to_csv(snakemake.output.cumulative_costs)

    else:
        # Myopic mode: Multiple networks via NetworkCollection
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

        nc = NetworkCollection(
            networks, index=pd.Index(planning_horizons, name="horizon")
        )
        logger.debug(
            f"Created NetworkCollection with horizons: {list(planning_horizons)}"
        )

        costs_df = None
        for output in OUTPUTS:
            if output == "cumulative_costs":
                continue  # Handle separately after costs are calculated
            logger.debug(f"Calculating {output}")
            result = globals()["calculate_" + output + "_collection"](nc)
            if output == "costs":
                costs_df = result  # Save for cumulative costs calculation
            result.to_csv(snakemake.output[output])

        # Calculate cumulative costs from costs DataFrame
        if costs_df is not None and len(planning_horizons) > 1:
            logger.debug("Calculating cumulative_costs")
            cumulative_costs = calculate_cumulative_costs(costs_df, planning_horizons)
            cumulative_costs.to_csv(snakemake.output.cumulative_costs)

    logger.debug("Summary calculation completed successfully")
