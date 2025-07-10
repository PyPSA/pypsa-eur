# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Extends cost data with custom cost modifications.

Inputs
------

- ``resources/costs_{planning_horizons}.csv``: Default cost data for specified planning horizon
- ``data/custom_costs.csv``: Custom cost modifications

Outputs
-------

- ``resources/costs_{planning_horizons}_extended.csv``: Extended cost data with custom modifications applied
"""

import pandas as pd


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


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_cost_data", planning_horizons=2030)

    # Retrieve costs assumptions
    planning_horizon = str(snakemake.wildcards.planning_horizons)
    costs = pd.read_csv(snakemake.input.costs)
    custom_costs = (
        pd.read_csv(snakemake.input.custom_costs, dtype={"planning_horizon": "object"})
        .query("planning_horizon in [@planning_horizon, 'all']")
        .drop("planning_horizon", axis=1)
    )

    if snakemake.params.custom_costs and not custom_costs.empty:
        # Expand "all" technologies across all available technologies from default costs
        custom_costs_expanded = expand_all_technologies(costs, custom_costs)

        # Combine costs assumptions with custom_costs superseding default costs values
        costs_extended = pd.concat(
            [costs, custom_costs_expanded], ignore_index=True
        ).drop_duplicates(subset=["technology", "parameter"], keep="last")
    else:
        costs_extended = costs

    costs_extended.to_csv(snakemake.output[0], index=False)
