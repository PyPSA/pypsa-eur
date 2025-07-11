# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Prepare and extend default cost data with custom cost modifications. Custom costs can target all planning horizons and / or technologies using the 'all' identifier.

Inputs
------

- ``resources/costs_{planning_horizons}.csv``: Default cost data for specified planning horizon
- ``data/custom_costs.csv``: Custom cost modifications

Outputs
-------

- ``resources/costs_{planning_horizons}_prepped.csv``: Prepared cost data with custom modifications applied
"""

import logging
import warnings

import pandas as pd
from _helpers import get_snapshots

from scripts.add_electricity import calculate_annuity

logger = logging.getLogger(__name__)


def expand_all_technologies(
    costs: pd.DataFrame, custom_costs: pd.DataFrame
) -> pd.DataFrame:
    """
    Expand 'all' technology entries in custom_costs to all available technologies from costs.

    Parameters
    ----------
    costs : pd.DataFrame
        DataFrame containing default costs with all available technologies
    custom_costs : pd.DataFrame
        DataFrame containing custom costs with potential 'all' technology entries

    Returns
    -------
    pd.DataFrame
        Expanded custom_costs DataFrame with 'all' entries replaced by specific technologies
    """
    all_entries = (
        custom_costs[custom_costs["technology"] == "all"]
        .drop("technology", axis=1)
        .assign(_merge_key=1)
    )
    specific_entries = custom_costs[custom_costs["technology"] != "all"]

    if all_entries.empty:
        return specific_entries

    # Perform cartesian product merge to expand all technologies
    available_technologies = (
        costs[["technology"]].drop_duplicates().assign(_merge_key=1)
    )
    expanded_df = pd.merge(
        all_entries,
        available_technologies,
        on="_merge_key",
    ).drop("_merge_key", axis=1)

    # Remove duplicates based on technology and parameter, keeping the last (specific entries)
    custom_costs_expanded = pd.concat(
        [expanded_df, specific_entries], ignore_index=True
    ).drop_duplicates(subset=["technology", "parameter"], keep="last")

    return custom_costs_expanded


def prepare_costs(
    costs: pd.DataFrame, config: dict, max_hours: dict = None, nyears: float = 1.0
) -> pd.DataFrame:
    """
    Standardize and prepare extended costs data.

    Parameters
    ----------
    costs : pd.DataFrame
        DataFrame containing extended costs
    config : dict
        Dictionary containing cost-related configuration parameters
    max_hours : dict, optional
        Dictionary specifying maximum hours for storage technologies
    nyears : float, optional
        Number of years for investment, by default 1.0

    Returns
    -------
    costs : pd.DataFrame
        DataFrame containing the prepared cost data

    """
    # Copy marginal_cost and capital_cost for backward compatibility
    for key in ("marginal_cost", "capital_cost"):
        if key in config:
            config["overwrites"][key] = config[key]

    # set all asset costs and other parameters
    costs = costs.set_index(["technology", "parameter"]).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.loc[costs.unit.str.contains("/GW"), "value"] /= 1e3

    costs.unit = costs.unit.str.replace("/kW", "/MW")
    costs.unit = costs.unit.str.replace("/GW", "/MW")

    # min_count=1 is important to generate NaNs which are then filled by fillna
    costs = costs.value.unstack(level=1).groupby("technology").sum(min_count=1)
    costs = costs.fillna(config["fill_values"])

    # Process overwrites for various attributes
    for attr in ("investment", "lifetime", "FOM", "VOM", "efficiency", "fuel"):
        overwrites = config["overwrites"].get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites
            warnings.warn(
                "Config-based cost overwrites is deprecated. Use external 'data/custom_costs.csv' file instead.",
                DeprecationWarning,
            )
            logger.info(f"Overwriting {attr} with:\n{overwrites}")

    annuity_factor = calculate_annuity(costs["lifetime"], costs["discount rate"])
    annuity_factor_fom = annuity_factor + costs["FOM"] / 100.0
    costs["capital_cost"] = annuity_factor_fom * costs["investment"] * nyears

    costs.at["OCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["CCGT", "fuel"] = costs.at["gas", "fuel"]

    costs["marginal_cost"] = costs["VOM"] + costs["fuel"] / costs["efficiency"]

    costs.at["OCGT", "CO2 intensity"] = costs.at["gas", "CO2 intensity"]
    costs.at["CCGT", "CO2 intensity"] = costs.at["gas", "CO2 intensity"]

    costs.at["solar", "capital_cost"] = costs.at["solar-utility", "capital_cost"]
    costs = costs.rename({"solar-utility single-axis tracking": "solar-hsat"})

    # Calculate storage costs if max_hours is provided
    if max_hours is not None:

        def costs_for_storage(store, link1, link2=None, max_hours=1.0):
            capital_cost = link1["capital_cost"] + max_hours * store["capital_cost"]
            if link2 is not None:
                capital_cost += link2["capital_cost"]
            return pd.Series(
                {
                    "capital_cost": capital_cost,
                    "marginal_cost": 0.0,
                    "CO2 intensity": 0.0,
                }
            )

        costs.loc["battery"] = costs_for_storage(
            costs.loc["battery storage"],
            costs.loc["battery inverter"],
            max_hours=max_hours["battery"],
        )
        costs.loc["H2"] = costs_for_storage(
            costs.loc["hydrogen storage underground"],
            costs.loc["fuel cell"],
            costs.loc["electrolysis"],
            max_hours=max_hours["H2"],
        )

    for attr in ("marginal_cost", "capital_cost"):
        overwrites = config["overwrites"].get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            idx = overwrites.index.intersection(costs.index)
            costs.loc[idx, attr] = overwrites.loc[idx]
            warnings.warn(
                "Config-based cost overwrites is deprecated. Use external 'data/custom_costs.csv' file instead.",
                DeprecationWarning,
            )
            logger.info(f"Overwriting {attr} with:\n{overwrites}")

    return costs


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_cost_data", planning_horizons=2030)

    config = snakemake.params.costs

    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day, tz="UTC"
    )
    nyears = len(snapshots) / 8760
    planning_horizon = str(snakemake.wildcards.planning_horizons)

    # Retrieve costs assumptions
    costs = pd.read_csv(snakemake.input.costs)
    custom_costs = (
        pd.read_csv(snakemake.input.custom_costs, dtype={"planning_horizon": "object"})
        .query("planning_horizon in [@planning_horizon, 'all']")
        .drop("planning_horizon", axis=1)
    )

    if config.get("custom_costs", False) and not custom_costs.empty:
        # Expand "all" technologies across all available technologies from default costs
        custom_costs_expanded = expand_all_technologies(costs, custom_costs)

        # Combine default costs assumptions with custom costs superseding default values
        costs_extended = pd.concat(
            [costs, custom_costs_expanded], ignore_index=True
        ).drop_duplicates(subset=["technology", "parameter"], keep="last")
    else:
        costs_extended = costs

    # Prepare costs
    costs_prepped = prepare_costs(
        costs_extended, config, snakemake.params.max_hours, nyears
    )

    costs_prepped.to_csv(snakemake.output[0])
