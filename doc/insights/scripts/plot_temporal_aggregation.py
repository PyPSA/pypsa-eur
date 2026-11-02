# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Plot the effect of temporal aggregation on system costs, capacities and solver
effort for the insights page ``doc/insights/temporal-aggregation.md``.

Usage::

    pixi run python doc/insights/scripts/plot_temporal_aggregation.py \
        --results results/time-aggregation --output doc/img/temporal-aggregation.svg

Pass ``--csv <path>`` to also export the underlying data.
"""

import argparse
import logging
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

logger = logging.getLogger(__name__)

HOURLY = "full"
EXCLUDED_STORES = ["oil", "methanol"]
REFERENCE_CPU = "AMD EPYC 9475F"

METHODS = {
    "segmentation": ("segmentation", "#2993b5"),
    "averaging": ("averaging", "#D10A49"),
}

GROUPS = {
    "solar": ("#f9d002", r"^solar"),
    "onshore wind": ("#235ebc", r"^onwind$"),
    "offshore wind": ("#6895dd", r"^offwind"),
    "nuclear": ("#ff8c00", r"^nuclear$"),
    "electricity grids": ("#6c9459", r"^(AC|DC|electricity distribution grid)$"),
    "battery": ("#ace37f", r"battery"),
    "hydrogen": ("#bf13a0", r"^H2 "),
    "heating": ("#d15959", r"heat pump|resistive heater|water pits|water tanks"),
    "gas-to-power/heat": ("#a85522", r"OCGT|CCGT|gas CHP|gas boiler"),
    "e-fuels": (
        "#25c49a",
        r"Fischer-Tropsch|methanolisation|electrobiofuels|Haber-Bosch|Sabatier"
        r"|^ammonia store$",
    ),
    "fuel supply": (
        "#baa741",
        r"^(gas|biogas|oil primary|solid biomass|biogas to gas)$",
    ),
    "other": ("#b4b4b4", r"."),
}

CAPACITY_GROUPS = [
    "solar",
    "onshore wind",
    "offshore wind",
    "battery",
    "hydrogen",
    "heating",
    "gas-to-power/heat",
]
STORAGE_GROUPS = ["battery", "hydrogen", "heating", "e-fuels"]
EXCLUDED_CAPACITIES = r"pipeline|(?<!dis)charger|water (pits|tanks) discharger"
MAX_HOURS_BARS = 12


def to_group(carrier: str) -> str:
    """Map a carrier name to its technology group."""
    return next(g for g, (_, pattern) in GROUPS.items() if re.search(pattern, carrier))


def parse_scenario(name: str) -> tuple[str, int]:
    """Split a scenario name like ``averaging-3h`` into method and hours."""
    if name == HOURLY:
        return "hourly", 1
    method, hours = name.rsplit("-", 1)
    return method, int(hours.removesuffix("h"))


def read_scenario(path: Path) -> dict[str, pd.Series]:
    """
    Read cost and capacity summaries and solver statistics of one scenario.

    Oil and methanol stores are excluded because the barrier solution of the
    hourly model leaves them oversized without crossover.
    """
    costs = pd.read_csv(path / "csvs/costs.csv", index_col=[0, 1, 2]).iloc[:, 0]
    excluded = costs.index.get_level_values("component") == "Store"
    excluded &= costs.index.get_level_values("carrier").isin(EXCLUDED_STORES)
    costs = costs[~excluded]
    costs = costs.groupby(costs.index.get_level_values("carrier").map(to_group)).sum()

    capacities = pd.read_csv(path / "csvs/capacities.csv", index_col=[0, 1]).iloc[:, 0]
    component = capacities.index.get_level_values("component")
    carrier = capacities.index.get_level_values("carrier")
    groups = carrier.map(to_group)
    is_store = (component == "Store") & ~carrier.isin(EXCLUDED_STORES)
    is_power = component.isin(["Generator", "Link"]) & ~carrier.str.contains(
        EXCLUDED_CAPACITIES
    )
    power = capacities[is_power].groupby(groups[is_power]).sum()
    storage = capacities[is_store].groupby(groups[is_store]).sum()

    log = (path / "logs/solve_network/solver_2050.log").read_text()
    solved = re.search(
        r"Barrier solved .* ([\d.]+) seconds \(([\d.]+) work units\)", log
    )
    benchmark = pd.read_csv(path / "benchmarks/solve_network_2050.log", sep="\t")

    return {
        "costs": costs / 1e9,
        "power": power.reindex(CAPACITY_GROUPS, fill_value=0) / 1e3,
        "storage": storage.reindex(STORAGE_GROUPS, fill_value=0) / 1e6,
        "solver": pd.Series(
            {
                "seconds": float(solved.group(1)),
                "work_units": float(solved.group(2)),
                "reference_cpu": float(f"CPU model: {REFERENCE_CPU}" in log),
                "memory": benchmark["max_rss"].iloc[0] / 1e3,
            }
        ),
    }


def load_results(results: Path) -> dict[str, pd.DataFrame]:
    """
    Collect all scenarios into frames indexed by method and hours.

    Solving times are normalised to one CPU type with the hardware-independent
    work units reported by the solver.
    """
    paths = [
        p
        for p in results.iterdir()
        if parse_scenario(p.name)[0] in [*METHODS, "hourly"]
    ]
    for p in paths:
        if not (p / "csvs/costs.csv").exists():
            logger.warning("Skipping unfinished scenario %s", p.name)
    scenarios = {
        parse_scenario(p.name): read_scenario(p)
        for p in paths
        if (p / "csvs/costs.csv").exists()
    }
    data = {
        key: pd.DataFrame({s: v[key] for s, v in scenarios.items()}).T.fillna(0)
        for key in ["costs", "power", "storage", "solver"]
    }
    data["costs"] = data["costs"][[g for g in GROUPS if g in data["costs"]]]
    solver = data["solver"]
    on_reference = solver.query("reference_cpu == 1")
    seconds_per_unit = (on_reference.seconds / on_reference.work_units).median()
    solver["time"] = solver.work_units * seconds_per_unit / 3600
    for df in data.values():
        df.index.names = ["method", "hours"]
    return data


def with_hourly(df: pd.DataFrame | pd.Series, method: str):
    """Select one method and prepend the hourly reference."""
    return pd.concat([df.loc[["hourly"]], df.loc[[method]]]).droplevel("method")


def plot_cost_error(axes: dict, data: dict, ylim: tuple):
    """Top row: total cost difference over aggregation, solving time and memory."""
    total = data["costs"].sum(axis=1)
    ref = total["hourly", 1]
    diff = total - ref
    pct = 100 * diff / ref
    solver = data["solver"]

    ax = axes["cost"]
    ax.set_ylim(ylim)
    for method, (label, color) in METHODS.items():
        d = with_hourly(diff, method).sort_index()
        ax.plot(d.index, d, "o-", color=color, clip_on=False, label=label)
    ax.legend(loc="lower left")
    hours = diff.index.get_level_values("hours")
    ax.set_xscale("log")
    ticks = [h for h in [1, 2, 3, 4, 6, 8, 12, 24, 48, 96] if h <= hours.max()]
    ax.set_xticks(ticks, [f"{t}" for t in ticks])
    ax.minorticks_off()
    ax.set_xlim(0.9, 1.1 * hours.max())
    ax.set_xlabel("aggregation (hours per snapshot, log scale)")
    ax.set_ylabel("cost difference (bn€/a)")
    ax.set_title("Total system cost")
    secondary = ax.secondary_yaxis(
        "right", functions=(lambda y: 100 * y / ref, lambda p: p * ref / 100)
    )
    secondary.set_ylabel("cost difference (%)")
    secondary.spines["right"].set_visible(True)

    ylim_pct = tuple(100 * y / ref for y in ylim)
    for name, column, xlabel, title in [
        ("time", "time", "solving time (h, log scale)", "Solving time"),
        ("memory", "memory", "peak memory (GB)", "Memory"),
    ]:
        ax = axes[name]
        for method, (_, color) in METHODS.items():
            x = with_hourly(solver[column], method).sort_index()
            y = with_hourly(pct, method).sort_index()
            ax.plot(x, y, "o-", color=color, clip_on=False)
        ax.plot(solver.loc[("hourly", 1), column], 0, "o", color="#20242C", zorder=3)
        ax.annotate(
            "hourly",
            (solver.loc[("hourly", 1), column], 0),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            fontsize=6,
        )
        ax.set_ylim(ylim_pct)
        ax.set_xlabel(xlabel)
        ax.set_title(title)
    axes["time"].set_xscale("log")
    ticks = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]
    ticks = [t for t in ticks if solver.time.min() / 2 < t < 2 * solver.time.max()]
    axes["time"].set_xticks(ticks, [f"{t:g}" for t in ticks])
    axes["time"].minorticks_off()
    axes["memory"].set_xlim(0, None)
    axes["time"].set_ylabel("cost difference (%)")
    axes["memory"].tick_params(labelleft=False)


def plot_baseline(ax, costs: pd.Series):
    """Stacked bar of the reference cost composition, doubling as colour legend."""
    left = 0
    for group, value in costs.items():
        ax.barh(0, value, left=left, color=GROUPS[group][0], label=group, height=0.6)
        left += value
    ax.set_xlim(0, left)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(False)
    ax.set_xlabel("system cost (bn€/a)")
    ax.set_title("Cost composition of the hourly model")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.9),
        ncols=6,
        handlelength=1,
        columnspacing=1.2,
    )


def plot_differences(axes: list, df: pd.DataFrame, xlabel: str):
    """Diverging stacked bars of differences to the hourly model, one axis per method."""
    diff = (df - df.loc["hourly", 1]).drop("hourly", level="method")
    diff = diff[diff.index.get_level_values("hours") <= MAX_HOURS_BARS]
    pos = diff.clip(lower=0).sum(axis=1).max()
    neg = diff.clip(upper=0).sum(axis=1).min()
    span = 1.1 * max(pos, -neg)

    for ax, (method, (label, _)) in zip(axes, METHODS.items()):
        d = diff.loc[method].sort_index()
        y = range(len(d))
        for sign in [1, -1]:
            part = (sign * d).clip(lower=0) * sign
            left = part.cumsum(axis=1) - part
            for group in d.columns:
                ax.barh(
                    y, part[group], left=left[group], color=GROUPS[group][0], height=0.7
                )
        net = d.sum(axis=1)
        ax.plot(net.clip(-span, span), y, "o", color="#20242C", ms=2.5, zorder=3)
        for yi, (lo, hi) in enumerate(
            zip(d.clip(upper=0).sum(axis=1), d.clip(lower=0).sum(axis=1))
        ):
            if lo < -span:
                ax.plot(-span, yi, "<", color="#20242C", ms=3, clip_on=False)
            if hi > span:
                ax.plot(span, yi, ">", color="#20242C", ms=3, clip_on=False)
        ax.set_title(label, loc="right", color=METHODS[method][1])
        ax.axvline(0, color="#888888", lw=0.6)
        ax.set_xlim(-span, span)
        ax.set_yticks(list(y), [f"{h} h" for h in d.index])
        ax.invert_yaxis()
        ax.grid(axis="y", visible=False)
    for ax in axes:
        ax.set_xlabel(xlabel)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)


def plot(data: dict, output: Path):
    """Assemble the multi-panel figure."""
    total = data["costs"].sum(axis=1)
    lowest = (total - total["hourly", 1]).min()
    ylim = (1.15 * lowest, -0.1 * lowest)

    mosaic = [
        ["cost"] * 2 + ["time"] * 2 + ["memory"] * 2,
        ["base"] * 6,
        ["costs_s"] * 3 + ["costs_a"] * 3,
        ["power_s"] * 3 + ["power_a"] * 3,
        ["storage_s"] * 3 + ["storage_a"] * 3,
    ]
    style = Path(__file__).with_name("insights.mplstyle")
    with plt.style.context(style):
        fig, axes = plt.subplot_mosaic(
            mosaic,
            figsize=(7.5, 10),
            layout="constrained",
            height_ratios=[1.3, 0.35, 1, 1, 1],
        )
        plot_cost_error(axes, data, ylim)
        plot_baseline(axes["base"], data["costs"].loc["hourly", 1])
        rows = [
            ("costs", "cost difference (bn€/a)", "Cost by technology"),
            ("power", "capacity difference (GW)", "Power capacity"),
            ("storage", "capacity difference (TWh)", "Storage capacity"),
        ]
        for key, xlabel, title in rows:
            row = [axes[f"{key}_{m[0]}"] for m in METHODS]
            plot_differences(row, data[key], xlabel)
            row[0].annotate(
                title,
                (0, 1),
                xycoords="axes fraction",
                xytext=(0, 16),
                textcoords="offset points",
                fontweight="bold",
                annotation_clip=False,
            )
        fig.savefig(output)
        plt.close(fig)
    logger.info("Saved figure to %s", output)


def export_csv(data: dict, path: Path):
    """Write all plotted data to one long-format CSV file."""
    frames = {k: df.stack().rename("value") for k, df in data.items()}
    pd.concat(frames, names=["quantity"]).to_csv(path)
    logger.info("Saved data to %s", path)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results", type=Path, default=Path("results/time-aggregation")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("doc/img/temporal-aggregation.svg")
    )
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()

    data = load_results(args.results)
    plot(data, args.output)
    if args.csv:
        export_csv(data, args.csv)
