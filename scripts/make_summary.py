# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create summary CSV files for all scenario runs including costs, capacities,
capacity factors, curtailment, energy balances, prices and other metrics.
"""

import logging

import pandas as pd
import pypsa

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
    pd.Series
        MultiIndex Series with levels ["component", "carrier"]
    """
    return n.statistics.energy_balance(groupby="carrier").sort_values(ascending=False)


def calculate_energy_balance(n: pypsa.Network) -> pd.Series:
    """
    Calculate the energy supply (positive) and consumption (negative) by technology carrier for each bus carrier.

    Returns
    -------
    pd.Series
        MultiIndex Series with levels ["component", "carrier", "bus_carrier"]

    Examples
    --------
    >>> eb = calculate_energy_balance(n)
    >>> eb.xs("methanol", level='bus_carrier')
    """
    return n.statistics.energy_balance().sort_values(ascending=False)


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
    """
    return (
        n.statistics.market_value(bus_carrier="AC", aggregate_across_components=True)
        .sort_values()
        .dropna()
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "make_summary",
            clusters="5",
            opts="",
            sector_opts="",
            planning_horizons="2030",
            configfiles="config/test/config.overnight.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    assign_carriers(n)
    assign_locations(n)

    pypsa.options.set_option("params.statistics.nice_names", False)
    pypsa.options.set_option("params.statistics.drop_zero", False)

    for output in OUTPUTS:
        globals()["calculate_" + output](n).to_csv(snakemake.output[output])
