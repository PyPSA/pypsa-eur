# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates plots from summary CSV files.
"""

import logging

import matplotlib.pyplot as plt
import pandas as pd
from _helpers import configure_logging, set_scenario_config
from plot_summary import preferred_order, rename_techs

logger = logging.getLogger(__name__)
plt.style.use("ggplot")


def plot_costs(cost_df, drop=None):

    df = cost_df.groupby(cost_df.index.get_level_values(2)).sum()

    # convert to billions
    df = df / 1e9

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.max(axis=1) < snakemake.config["plotting"]["costs_threshold"]]

    logger.info(
        f"Dropping technology with costs below {snakemake.config['plotting']['costs_threshold']} EUR billion per year"
    )
    logger.debug(df.loc[to_drop])

    df = df.drop(to_drop)

    logger.info(f"Total system cost of {round(df.sum().iloc[0])} EUR billion per year")

    new_index = preferred_order.intersection(df.index).append(
        df.index.difference(preferred_order)
    )

    # new_columns = df.sum().sort_values().index
    if drop != None:
        df = df.droplevel([1, 2, 3], axis=1)

    planning_horizons = (
        df.columns.get_level_values("planning_horizon").unique().sort_values()
    )
    scenarios = df.columns.get_level_values(0).unique()

    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(planning_horizons),
        figsize=(12, 8),
        sharey=True,  # This ensures that all subplots share the same y-axis
    )

    if len(planning_horizons) == 1:
        axes = [axes]

    for ax, year in zip(axes, planning_horizons):
        subset = df.xs(year, level="planning_horizon", axis=1)

        subset.T[new_index].plot(
            kind="bar",
            ax=ax,
            stacked=True,
            legend=False,
            color=[snakemake.config["plotting"]["tech_colors"][i] for i in new_index],
        )

        # Set title and x-label
        ax.set_title(year)

        ax.set_xlabel("")

        # Set x-ticks as scenario names (level=0)
        ax.set_xticks(range(len(scenarios)))

        if ax == axes[0]:
            ax.set_ylabel("System Cost [EUR billion per year]")

        ax.grid(axis="x")

        ax.set_ylim([0, snakemake.config["plotting"]["costs_max"]])

    handles, labels = ax.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()

    axes[-1].legend(
        handles,
        labels,
        ncol=1,
        loc="upper left",
        bbox_to_anchor=[1, 1],
        frameon=False,
    )

    plt.tight_layout()

    fig.savefig(snakemake.output.costs, bbox_inches="tight")

    df.sum().unstack().T.plot()
    plt.ylabel("System Cost [EUR billion per year]")
    plt.xlabel("planning horizon")
    plt.legend(bbox_to_anchor=(1, 1))
    plt.savefig(
        snakemake.output.costs.split(".svg")[0] + "-total.svg", bbox_inches="tight"
    )


def plot_balances(balances_df, drop=None):

    co2_carriers = ["co2", "co2 stored", "process emissions"]
    balances = {i.replace(" ", "_"): [i] for i in balances_df.index.levels[0]}
    balances["energy"] = [
        i for i in balances_df.index.levels[0] if i not in co2_carriers
    ]

    for k, v in balances.items():
        df = balances_df.loc[v]
        df = df.groupby(df.index.get_level_values(2)).sum()

        # convert MWh to TWh
        df = df / 1e6

        # remove trailing link ports
        df.index = [
            (
                i[:-1]
                if (
                    (i not in ["co2", "NH3", "H2"])
                    and (i[-1:] in ["0", "1", "2", "3", "4"])
                )
                else i
            )
            for i in df.index
        ]

        df = df.groupby(df.index.map(rename_techs)).sum()

        to_drop = df.index[
            df.abs().max(axis=1) < snakemake.config["plotting"]["energy_threshold"] / 10
        ]

        units = "MtCO2/a" if v[0] in co2_carriers else "TWh/a"
        logger.debug(
            f"Dropping technology energy balance smaller than {snakemake.config['plotting']['energy_threshold']/10} {units}"
        )
        logger.debug(df.loc[to_drop])

        df = df.drop(to_drop)

        logger.debug(
            f"Total energy balance for {v} of {round(df.sum().iloc[0],2)} {units}"
        )

        if df.empty:
            continue

        new_index = preferred_order.intersection(df.index).append(
            df.index.difference(preferred_order)
        )

        if drop != None:
            df = df.droplevel([1, 2, 3], axis=1)

        planning_horizons = (
            df.columns.get_level_values("planning_horizon").unique().sort_values()
        )
        scenarios = df.columns.get_level_values(0).unique()

        fig, axes = plt.subplots(
            nrows=1,
            ncols=len(planning_horizons),
            figsize=(12, 8),
            sharey=True,  # This ensures that all subplots share the same y-axis
        )

        if len(planning_horizons) == 1:
            axes = [axes]

        for ax, year in zip(axes, planning_horizons):
            subset = df.xs(year, level="planning_horizon", axis=1)

            subset.T[new_index].plot(
                kind="bar",
                ax=ax,
                stacked=True,
                legend=False,
                color=[
                    snakemake.config["plotting"]["tech_colors"][i] for i in new_index
                ],
            )

            # Set title and x-label
            ax.set_title(year)

            # Set x-ticks as scenario names (level=0)
            ax.set_xticks(range(len(scenarios)))

            if ax == axes[0]:
                if v[0] in co2_carriers:
                    ax.set_ylabel("CO2 [MtCO2/a]")
                else:
                    ax.set_ylabel("Energy [TWh/a]")

            ax.grid(axis="x")

        handles, labels = ax.get_legend_handles_labels()
        handles.reverse()
        labels.reverse()

        axes[-1].legend(
            handles,
            labels,
            ncol=1,
            loc="upper left",
            bbox_to_anchor=[1, 1],
            frameon=False,
        )

        plt.tight_layout()

        fig.savefig(snakemake.output.balances[:-10] + k + ".svg", bbox_inches="tight")


# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("plot_summary_all")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n_header = 4

    prefix = snakemake.config["run"]["prefix"]
    path = snakemake.output[0].split("graphs")[0]
    scenarios = snakemake.config["run"]["name"]

    costs = {}
    balances = {}
    for scenario in scenarios:
        try:
            costs[scenario] = pd.read_csv(
                f"{path}/{scenario}/csvs/costs.csv",
                index_col=list(range(3)),
                header=list(range(n_header)),
            )
            balances[scenario] = pd.read_csv(
                f"{path}/{scenario}/csvs/supply_energy.csv",
                index_col=list(range(3)),
                header=list(range(n_header)),
            )
        except FileNotFoundError:
            logger.info(f"{scenario} not solved yet.")

    costs = pd.concat(costs, axis=1)
    balances = pd.concat(balances, axis=1)

    plot_costs(costs, drop=True)

    plot_balances(balances, drop=True)
