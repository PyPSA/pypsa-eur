# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Creates plots from summary CSV files.
"""

import logging
import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pandas as pd

from scripts._helpers import (
    configure_logging,
    create_placeholder_plot,
    rename_techs,
    set_scenario_config,
)
from scripts.prepare_sector_network import determine_emission_sectors

logger = logging.getLogger(__name__)
plt.style.use("bmh")


# consolidate and rename

preferred_order = pd.Index(
    [
        "transmission lines",
        "hydroelectricity",
        "hydro reservoir",
        "run of river",
        "pumped hydro storage",
        "solid biomass",
        "biogas",
        "onshore wind",
        "offshore wind",
        "offshore wind (AC)",
        "offshore wind (DC)",
        "solar PV",
        "solar thermal",
        "solar rooftop",
        "solar",
        "building retrofitting",
        "ground heat pump",
        "air heat pump",
        "heat pump",
        "resistive heater",
        "power-to-heat",
        "gas-to-power/heat",
        "CHP",
        "OCGT",
        "gas boiler",
        "gas",
        "natural gas",
        "methanation",
        "ammonia",
        "hydrogen storage",
        "power-to-gas",
        "power-to-liquid",
        "battery storage",
        "hot water storage",
        "CO2 sequestration",
    ]
)


def check_tech_colors(tech_colors, keys):
    """
    Check if all keys exist in tech_colors mapping, otherwise raise KeyError.
    """
    missing = [k for k in keys if k not in tech_colors]
    if missing:
        raise KeyError(
            f"The following technology carrier(s) do not have a defined color in the plotting configuration: {missing}"
        )


def plot_costs():
    cost_df = pd.read_csv(
        snakemake.input.costs, index_col=list(range(3)), header=list(range(n_header))
    )

    df = cost_df.groupby("carrier").sum()

    # convert to billions
    df = df / 1e9

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.max(axis=1) < snakemake.params.plotting["costs_threshold"]]

    logger.debug(
        f"Dropping technology with costs below {snakemake.params['plotting']['costs_threshold']} EUR billion per year"
    )
    logger.debug(df.loc[to_drop])

    df = df.drop(to_drop)

    total_cost = df.sum().iloc[0]
    logger.debug(f"Total system cost of {total_cost:.2f} EUR billion per year")

    # Check if there's any data left to plot
    if df.empty:
        logger.warning(
            "No cost data to plot after filtering. Creating placeholder plot."
        )
        create_placeholder_plot(
            snakemake.output.costs,
            "No cost data available\n(all costs below threshold)",
            ylabel="System Cost [EUR billion per year]",
        )
        return  # Early return is OK here since other functions will still be called

    new_index = preferred_order.intersection(df.index).append(
        df.index.difference(preferred_order)
    )

    # new_columns = df.sum().sort_values().index

    check_tech_colors(snakemake.params.plotting["tech_colors"], new_index)

    fig, ax = plt.subplots(figsize=(12, 8))

    df.loc[new_index].T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[snakemake.params.plotting["tech_colors"][i] for i in new_index],
    )

    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    # Handle auto-scaling if configured
    costs_max = snakemake.params.plotting["costs_max"]
    if costs_max == "auto":
        costs_max = None
        logger.debug("Auto-scaling y-axis (costs_max='auto')")

    ax.set_ylim([0, costs_max])

    ax.set_ylabel("System Cost [EUR billion per year]")

    ax.set_xlabel("")

    ax.grid(axis="x")

    ax.legend(
        handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1, 1], frameon=False
    )

    fig.savefig(snakemake.output.costs, bbox_inches="tight")
    plt.close(fig)


def plot_energy():
    energy_df = pd.read_csv(
        snakemake.input.energy, index_col=list(range(2)), header=list(range(n_header))
    )

    df = energy_df.groupby("carrier").sum()

    # convert MWh to TWh
    df = df / 1e6

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[
        df.abs().max(axis=1) < snakemake.params.plotting["energy_threshold"]
    ]

    logger.debug(
        f"Dropping all technology with energy consumption or production below {snakemake.params['plotting']['energy_threshold']} TWh/a"
    )
    logger.debug(df.loc[to_drop])

    df = df.drop(to_drop)

    total_energy = df.sum().iloc[0]
    logger.debug(f"Total energy of {total_energy:.2f} TWh/a")

    if df.empty:
        logger.warning(
            "No energy data to plot after filtering. Creating placeholder plot."
        )
        create_placeholder_plot(
            snakemake.output.energy,
            "No energy data available\n(all values below threshold)",
            ylabel="Energy [TWh/a]",
        )
        return

    new_index = preferred_order.intersection(df.index).append(
        df.index.difference(preferred_order)
    )

    # new_columns = df.columns.sort_values()

    check_tech_colors(snakemake.params.plotting["tech_colors"], new_index)

    fig, ax = plt.subplots(figsize=(12, 8))

    logger.debug(df.loc[new_index])

    df.loc[new_index].T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[snakemake.params.plotting["tech_colors"][i] for i in new_index],
    )

    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    # Handle auto-scaling if configured
    energy_min = snakemake.params.plotting["energy_min"]
    energy_max = snakemake.params.plotting["energy_max"]

    if energy_max == "auto":
        energy_max = None
        logger.debug("Auto-scaling y-axis max (energy_max='auto')")

    if energy_min == "auto":
        energy_min = None
        logger.debug("Auto-scaling y-axis min (energy_min='auto')")

    ax.set_ylim([energy_min, energy_max])

    ax.set_ylabel("Energy [TWh/a]")

    ax.set_xlabel("")

    ax.grid(axis="x")

    ax.legend(
        handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1, 1], frameon=False
    )

    fig.savefig(snakemake.output.energy, bbox_inches="tight")
    plt.close(fig)


def plot_balances():
    co2_carriers = ["co2", "co2 stored", "process emissions"]

    balances_df = pd.read_csv(
        snakemake.input.balances, index_col=list(range(3)), header=list(range(n_header))
    )

    balances = {k: df for k, df in balances_df.groupby("bus_carrier")}
    balances["energy"] = balances_df.groupby(["component", "carrier"]).sum()

    for bus_carrier, df in balances.items():
        df = df.groupby("carrier").sum()

        # convert MWh to TWh
        df = df / 1e6

        df = df.groupby(df.index.map(rename_techs)).sum()

        to_drop = df.index[
            df.abs().max(axis=1) < snakemake.params.plotting["energy_threshold"] / 10
        ]

        units = "MtCO2/a" if bus_carrier in co2_carriers else "TWh/a"
        logger.debug(
            f"Dropping technology energy balance smaller than {snakemake.params['plotting']['energy_threshold'] / 10} {units}"
        )
        logger.debug(df.loc[to_drop])

        df = df.drop(to_drop)

        logger.debug(
            f"Total energy balance for {bus_carrier} of {round(df.sum().iloc[0], 2)} {units}"
        )

        if df.empty:
            continue

        new_index = preferred_order.intersection(df.index).append(
            df.index.difference(preferred_order)
        )

        new_columns = df.columns.sort_values()

        check_tech_colors(snakemake.params.plotting["tech_colors"], new_index)

        fig, ax = plt.subplots(figsize=(12, 8))

        df.loc[new_index, new_columns].T.plot(
            kind="bar",
            ax=ax,
            stacked=True,
            color=[snakemake.params.plotting["tech_colors"][i] for i in new_index],
        )

        handles, labels = ax.get_legend_handles_labels()

        handles.reverse()
        labels.reverse()

        if bus_carrier in co2_carriers:
            ax.set_ylabel("CO2 [MtCO2/a]")
        else:
            ax.set_ylabel("Energy [TWh/a]")

        ax.set_xlabel("")

        ax.grid(axis="x")

        ax.legend(
            handles,
            labels,
            ncol=1,
            loc="upper left",
            bbox_to_anchor=[1, 1],
            frameon=False,
        )

        fig.savefig(
            snakemake.output.balances[:-10] + bus_carrier + ".pdf", bbox_inches="tight"
        )
        plt.close(fig)

    if not os.path.exists(snakemake.output.balances):
        logger.warning("No balance data was plotted. Creating placeholder file.")
        create_placeholder_plot(
            snakemake.output.balances,
            "No balance data available\n(all values below threshold)",
            ylabel="Energy [TWh/a]",
        )


def historical_emissions(
    co2_totals: str, countries: list[str], options: dict
) -> pd.Series:
    """
    Read historical emissions in Gt of the modelled sectors since 1990.
    """
    co2 = pd.read_csv(co2_totals, index_col=[0, 1]).loc[countries]
    sectors = determine_emission_sectors(options)
    return co2[sectors].sum(axis=1).groupby("year").sum().loc[1990:] / 1e3


def plot_carbon_budget_distribution(co2_totals, options):
    """
    Plot historical carbon emissions in the EU and decarbonization path.
    """
    import seaborn as sns

    sns.set()
    sns.set_style("ticks")
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.labelsize"] = 20
    plt.rcParams["ytick.labelsize"] = 20

    # historic emissions
    countries = snakemake.params.countries
    emissions = historical_emissions(co2_totals, countries, options)
    e_1990 = emissions[1990]

    if snakemake.config["foresight"] == "myopic":
        path_cb = "results/" + snakemake.params.RDIR + "/csvs/"
        co2_cap = pd.read_csv(path_cb + "carbon_budget_distribution.csv", index_col=0)[
            ["cb"]
        ]
        co2_cap *= e_1990
    else:
        supply_energy = pd.read_csv(
            snakemake.input.balances, index_col=[0, 1, 2], header=list(range(n_header))
        )
        if "co2" not in supply_energy.index.get_level_values("bus_carrier"):
            logger.warning(
                "No CO2 balance in energy balances; skipping carbon budget plot."
            )
            return
        co2_emissions = supply_energy.xs("co2", level="bus_carrier").droplevel(
            "component"
        )
        co2_cap = (
            co2_emissions.drop("co2").sum().div(1e9).to_frame(name="co2 emissions")
        )
        co2_cap.index = co2_cap.index.astype(int)

    plt.figure(figsize=(10, 7))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = plt.subplot(gs1[0, 0])
    ax1.set_ylabel("CO$_2$ emissions \n [Gt per year]", fontsize=22)
    # ax1.set_ylim([0, 5])
    ax1.set_xlim([1990, snakemake.params.planning_horizons[-1] + 1])

    ax1.plot(emissions, color="black", linewidth=3, label=None)

    # plot committed and under-discussion targets
    # (notice that historical emissions include all countries in the
    # network, but targets refer to EU)
    ax1.plot(
        [2020],
        [0.8 * emissions[1990]],
        marker="*",
        markersize=12,
        markerfacecolor="black",
        markeredgecolor="black",
    )

    ax1.plot(
        [2030],
        [0.45 * emissions[1990]],
        marker="*",
        markersize=12,
        markerfacecolor="black",
        markeredgecolor="black",
    )

    ax1.plot(
        [2030],
        [0.6 * emissions[1990]],
        marker="*",
        markersize=12,
        markerfacecolor="black",
        markeredgecolor="black",
    )

    ax1.plot(
        [2050, 2050],
        [x * emissions[1990] for x in [0.2, 0.05]],
        color="gray",
        linewidth=2,
        marker="_",
        alpha=0.5,
    )

    ax1.plot(
        [2050],
        [0.0 * emissions[1990]],
        marker="*",
        markersize=12,
        markerfacecolor="black",
        markeredgecolor="black",
        label="EU committed target",
    )

    for col in co2_cap.columns:
        ax1.plot(co2_cap[col], linewidth=3, label=col)

    ax1.legend(
        fancybox=True, fontsize=18, loc=(0.01, 0.01), facecolor="white", frameon=True
    )

    plt.grid(axis="y")
    path = snakemake.output.balances.split("balances")[0] + "carbon_budget.pdf"
    plt.savefig(path, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("plot_summary")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n_header = 1

    plot_costs()

    plot_energy()

    plot_balances()

    if snakemake.input.get("co2_totals"):
        options = snakemake.params.sector
        plot_carbon_budget_distribution(snakemake.input.co2_totals, options)
