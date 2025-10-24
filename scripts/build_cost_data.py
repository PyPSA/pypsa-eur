# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Prepare and extend default cost data with custom cost modifications. Custom costs can target all planning horizons and / or technologies using the 'all' identifier.

Preparing the cost data includes:
- aligning all units to conventional units (i.e. MW / MWh),
- filling in missing data,
- computing 'capital_cost' parameter (annualised investment costs and FOM),
- computing 'marginal_cost' parameter (fuel costs and VOM),
- computing storage costs for batteries and hydrogen,
- (deprecated) overwriting attributes using config-based modifications.

Inputs
------

- ``resources/costs_{planning_horizons}.csv``: Default cost data for specified planning horizon
- (by default) ``data/custom_costs.csv``: Custom cost modifications (can be configured with `costs:custom_costs:file`

Outputs
-------

- ``resources/costs_{planning_horizons}_processed.csv``: Prepared cost data with custom modifications applied
"""

import logging
import warnings

import pandas as pd
from _helpers import get_snapshots

from scripts.add_electricity import calculate_annuity

logger = logging.getLogger(__name__)

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
    for attr in (
        "investment",
        "lifetime",
        "FOM",
        "VOM",
        "efficiency",
        "fuel",
        "standing losses",
    ):
        overwrites = config["overwrites"].get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites
            warnings.warn(
                "Config-based cost overwrites is deprecated. Use external file instead (by default 'data/custom_costs.csv').",
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

    costs = costs.rename(columns={"standing losses": "standing_losses"})

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
                    "standing_losses": 0.0,
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
                "Config-based cost overwrites is deprecated. Use external file instead (by default 'data/custom_costs.csv').",
                DeprecationWarning,
            )
            logger.info(f"Overwriting {attr} with:\n{overwrites}")

    return costs


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_cost_data", planning_horizons=2030)

    config = snakemake.config["costs"]

    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day, tz="UTC"
    )
    nyears = len(snapshots) / 8760
    planning_horizon = str(snakemake.wildcards.planning_horizons)

    # Retrieve costs assumptions
    # If the index is the unique pair (technology, parameter),
    # it can be easily overwritten.
    costs = pd.read_csv(snakemake.input.costs, index_col=[0,1])
    if snakemake.input.costs is not None:
        custom_costs = pd.read_csv(
                snakemake.params.costs["custom_cost_fn"],
                dtype={"planning_horizon": "object"},
                index_col=[1,2] # 0 is planning_horizon
        ).query("planning_horizon in [@planning_horizon, 'all']")
        custom_costs = custom_costs.drop("planning_horizon", axis=1)

        # Overwrite unique pairs of (technology, paramter)
        missing_idx = custom_costs.index.difference(costs.index)
        if len(missing_idx) > 0:
            costs = pd.concat([costs, custom_costs.loc[missing_idx]])

        costs.loc[custom_costs.index] = custom_costs
        costs = costs.reset_index()

    # Prepare costs
    costs_processed = prepare_costs(
        costs, config, snakemake.params.max_hours, nyears
    )

    costs_processed.to_csv(snakemake.output[0])
