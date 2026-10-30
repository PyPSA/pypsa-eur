# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Plot computational scaling from a collected solve benchmark CSV.

Examples
--------
python scripts/plot_solve_benchmarks.py \
    results/dispatch-power-IT-noUC/benchmarks.csv

python scripts/plot_solve_benchmarks.py benchmarks.csv \
    --output-dir plots --formats png pdf --log-scale
"""

import argparse
import logging
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

logger = logging.getLogger(__name__)

TIME_RESOLUTION_RE = re.compile(
    r"(?:^|/)th_(?P<periods>\d+)c(?P<hours>\d+)(?P<variant>_[^/]*)?(?:\Z|/)"
)
REQUIRED_COLUMNS = {
    "case",
    "clusters",
    "solver_options",
    "s",
    "objective",
    "objective_constant",
}
OUTCOME_STYLES = {
    "infeasible": {"marker": "P", "label": "Infeasible", "size": 100},
    "time_limit": {"marker": "X", "label": "Time limit", "size": 130},
}
SOLVED_POINT_SIZE = 75
SMSPP_TIMING_STYLES = {
    "smspp_computational_s": {
        "label": "SMS++ computational",
        "linestyle": "--",
    },
    "smspp_optimization_s": {
        "label": "SMS++ subprocess",
        "linestyle": ":",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Collected benchmark CSV.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory (default: <input directory>/plots).",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf"],
        choices=["png", "pdf", "svg"],
        help="Output formats (default: png pdf).",
    )
    parser.add_argument(
        "--log-scale",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use a logarithmic runtime axis (default: enabled).",
    )
    parser.add_argument(
        "--title",
        default="Solver computational performance",
        help="Figure title.",
    )
    parser.add_argument(
        "--dpi", type=int, default=220, help="Raster output resolution (default: 220)."
    )
    return parser.parse_args()


def parse_time_resolution(case: object) -> tuple[int, int] | None:
    """Extract ``(number of periods, hours per period)`` from a case path."""
    match = TIME_RESOLUTION_RE.search(str(case))
    if match is None:
        return None
    return int(match["periods"]), int(match["hours"])


def parse_case_variant(case: object) -> str:
    """Extract the non-numeric suffix after a time-resolution component."""
    match = TIME_RESOLUTION_RE.search(str(case))
    if match is None or match["variant"] is None:
        return ""
    variant = match["variant"].lstrip("_")
    return re.sub(r"_\d+\Z", "", variant)


def parse_scenario_count(case: object) -> int:
    """Infer a trailing scenario count for CSVs without a scenarios column."""
    match = TIME_RESOLUTION_RE.search(str(case))
    if match is None or match["variant"] is None:
        return 1
    scenario_match = re.search(r"_(\d+)\Z", match["variant"])
    return int(scenario_match[1]) if scenario_match else 1


def prepare_data(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(sorted(missing))}"
        )

    parsed = data["case"].map(parse_time_resolution)
    unparsed_cases = sorted(data.loc[parsed.isna(), "case"].astype(str).unique())
    if unparsed_cases:
        logger.warning(
            "Ignoring %d row(s) whose case does not contain th_<periods>c<hours>: %s",
            int(parsed.isna().sum()),
            ", ".join(unparsed_cases),
        )
    data = data.loc[parsed.notna()].copy()
    if data.empty:
        raise ValueError("No time resolutions like 'th_24c24' were found in case.")

    data[["periods", "hours_per_period"]] = pd.DataFrame(
        parsed.dropna().tolist(), index=data.index
    )
    data["case_variant"] = data["case"].map(parse_case_variant)
    inferred_scenarios = data["case"].map(parse_scenario_count)
    data["scenarios"] = pd.to_numeric(
        data.get("scenarios", inferred_scenarios), errors="coerce"
    ).fillna(inferred_scenarios)
    data["modeled_hours"] = data["periods"] * data["hours_per_period"]
    data["clusters"] = pd.to_numeric(data["clusters"], errors="coerce")
    data["runtime_s"] = pd.to_numeric(data["s"], errors="coerce")
    # Snakemake reports max_rss in MiB; display peak resident memory in GiB.
    data["memory_gib"] = pd.to_numeric(
        data.get("max_rss", pd.Series(index=data.index, dtype=float)),
        errors="coerce",
    ) / 1024
    for column in SMSPP_TIMING_STYLES:
        data[column] = pd.to_numeric(
            data.get(column, pd.Series(index=data.index, dtype=float)),
            errors="coerce",
        )
    data["objective"] = pd.to_numeric(data["objective"], errors="coerce")
    data["objective_constant"] = pd.to_numeric(
        data["objective_constant"], errors="coerce"
    )
    data["objective"] += data["objective_constant"]
    data["time_limit_s"] = pd.to_numeric(
        data.get("time_limit_s", pd.Series(index=data.index, dtype=float)),
        errors="coerce",
    )
    termination = data.get(
        "termination_condition", pd.Series("", index=data.index, dtype=str)
    )
    status = data.get("solve_status", pd.Series("", index=data.index, dtype=str))
    outcome_text = (
        termination.fillna("")
        .astype(str)
        .str.cat(status.fillna("").astype(str), sep=" ")
    )
    data["outcome"] = "solved"
    data.loc[outcome_text.str.contains("infeasible", case=False), "outcome"] = (
        "infeasible"
    )
    time_limited = outcome_text.str.contains(
        r"time[_ -]?limit|stopped for time", case=False, regex=True
    )
    failed = (
        data.get("benchmark_file", pd.Series("", index=data.index))
        .fillna("")
        .str.endswith(".failed")
    )
    data.loc[failed | time_limited, "outcome"] = "time_limit"
    missing_limits = time_limited & data["time_limit_s"].isna()
    if missing_limits.any():
        logger.warning(
            "Ignoring %d time-limit row(s) with no time_limit_s; measured runtime is not reliable.",
            int(missing_limits.sum()),
        )
    data.loc[time_limited, "runtime_s"] = data.loc[time_limited, "time_limit_s"]
    data["solver_options"] = data["solver_options"].fillna(data.get("solver"))
    data["solver_options"] = data["solver_options"].fillna("unspecified").astype(str)

    valid = (
        data["clusters"].notna() & data["runtime_s"].notna() & (data["runtime_s"] > 0)
    )
    if not valid.all():
        logger.warning(
            "Ignoring %d row(s) with invalid clusters/runtime.", (~valid).sum()
        )
    data = data.loc[valid].copy()
    if data.empty:
        raise ValueError(
            "No rows have valid clusters and positive runtime in column 's'."
        )

    # Repeated runs of the same configuration are summarized robustly.
    group_columns = [
        "case_variant",
        "scenarios",
        "solver_options",
        "clusters",
        "periods",
        "hours_per_period",
        "modeled_hours",
        "outcome",
    ]
    return (
        data.groupby(group_columns, as_index=False, dropna=False)[
            ["runtime_s", *SMSPP_TIMING_STYLES, "objective", "memory_gib"]
        ]
        .median()
        .sort_values(["periods", "hours_per_period", "clusters", "solver_options"])
    )


def variant_slug(variant: str) -> str:
    """Return a filesystem-safe suffix for a case variant."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", variant).strip("_.-") or "unlabelled"


def resolution_label(periods: int, hours: int) -> str:
    unit = "day" if periods == 1 else "days"
    return f"{periods} × {hours} h ({periods} {unit})"


def seconds_label(seconds: float, _position: float = 0) -> str:
    return f"{seconds:g}"


def needs_log_scale(values: pd.Series, ratio_threshold: float = 20.0) -> bool:
    """Use logarithmic scaling when positive values span a wide range."""
    finite = pd.to_numeric(values, errors="coerce").dropna()
    positive = finite.loc[finite > 0]
    return len(positive) > 1 and positive.max() / positive.min() >= ratio_threshold


def solver_colors(solvers: list[str]) -> dict[str, object]:
    palette = plt.get_cmap("tab10")
    return {solver: palette(index % 10) for index, solver in enumerate(solvers)}


def style_axis(ax: plt.Axes, log_scale: bool) -> None:
    ax.grid(axis="y", color="#D9DEE7", linewidth=0.8)
    ax.grid(axis="x", color="#EDF0F4", linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#AAB2BF")
    if log_scale:
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=8))
        ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(seconds_label))


def draw_series(
    ax: plt.Axes,
    data: pd.DataFrame,
    x: str,
    solvers: list[str],
    colors: dict[str, object],
    y: str,
) -> None:
    for solver in solvers:
        subset = data.loc[data["solver_options"] == solver].sort_values(x)
        if subset.empty:
            continue
        # A time-limit point is an isolated observation, not part of a trend line.
        connected_y = subset[y].where(subset["outcome"] != "time_limit")
        ax.plot(
            subset[x],
            connected_y,
            color=colors[solver],
            marker=None,
            linewidth=2,
            label=solver,
        )
        if y == "runtime_s":
            for timing_column, timing_style in SMSPP_TIMING_STYLES.items():
                timing_subset = subset.dropna(subset=[timing_column])
                if timing_subset.empty:
                    continue
                ax.plot(
                    timing_subset[x],
                    timing_subset[timing_column],
                    color=colors[solver],
                    linestyle=timing_style["linestyle"],
                    marker=None,
                    linewidth=2,
                )
        for outcome, outcome_subset in subset.groupby("outcome"):
            style = OUTCOME_STYLES.get(
                outcome, {"marker": "o", "size": SOLVED_POINT_SIZE}
            )
            ax.scatter(
                outcome_subset[x],
                outcome_subset[y],
                color=colors[solver],
                marker=style["marker"],
                s=style["size"],
                edgecolor="white",
                linewidth=1.1,
                zorder=4,
            )


def figure_legend(
    figure: plt.Figure,
    data: pd.DataFrame,
    solvers: list[str],
    colors: dict[str, object],
) -> None:
    handles = [
        Line2D([], [], color=colors[s], marker="o", linewidth=2, label=s)
        for s in solvers
    ]
    if any(data[column].notna().any() for column in SMSPP_TIMING_STYLES):
        handles.append(
            Line2D([], [], color="#444444", linestyle="-", linewidth=2, label="Total")
        )
        for column, style in SMSPP_TIMING_STYLES.items():
            if data[column].notna().any():
                handles.append(
                    Line2D(
                        [],
                        [],
                        color="#444444",
                        linestyle=style["linestyle"],
                        linewidth=2,
                        label=style["label"],
                    )
                )
    for outcome, style in OUTCOME_STYLES.items():
        if outcome in data["outcome"].values:
            handles.append(
                Line2D(
                    [],
                    [],
                    color="#444444",
                    marker=style["marker"],
                    linestyle="none",
                    markersize=11 if outcome == "time_limit" else 9,
                    label=style["label"],
                )
            )
    figure.legend(
        handles=handles,
        loc="outside lower center",
        ncols=min(len(handles), 4),
        frameon=False,
        title="Solver options / timing / outcome",
    )


def make_node_scaling_figure(
    data: pd.DataFrame,
    solvers: list[str],
    colors: dict[str, object],
    title: str,
    log_scale: bool,
    y: str = "runtime_s",
    y_label: str = "Wall-clock runtime [s]",
    heading: str | None = None,
) -> plt.Figure:
    resolutions = (
        data[["periods", "hours_per_period"]]
        .drop_duplicates()
        .sort_values(["periods", "hours_per_period"])
        .itertuples(index=False, name=None)
    )
    resolutions = list(resolutions)
    columns = min(3, len(resolutions))
    rows = math.ceil(len(resolutions) / columns)
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(5.1 * columns, 3.8 * rows),
        sharey=y == "runtime_s",
        squeeze=False,
        constrained_layout=True,
    )
    for ax, (periods, hours) in zip(axes.flat, resolutions):
        subset = data.loc[
            (data["periods"] == periods) & (data["hours_per_period"] == hours)
        ]
        draw_series(ax, subset, "clusters", solvers, colors, y)
        panel_log_scale = log_scale or (y == "objective" and needs_log_scale(subset[y]))
        style_axis(ax, panel_log_scale)
        ax.set_title(
            resolution_label(periods, hours), loc="left", fontweight="semibold"
        )
        ax.set_xlabel("Network nodes (clusters)")
        ax.set_xticks(sorted(subset["clusters"].unique()))
    for ax in axes[:, 0]:
        ax.set_ylabel(y_label)
    for ax in axes.flat[len(resolutions) :]:
        ax.set_visible(False)
    figure.suptitle(f"{heading or title}\nNetwork size", fontweight="bold", fontsize=15)
    figure_legend(figure, data, solvers, colors)
    return figure


def make_time_scaling_figure(
    data: pd.DataFrame,
    solvers: list[str],
    colors: dict[str, object],
    title: str,
    log_scale: bool,
    y: str = "runtime_s",
    y_label: str = "Wall-clock runtime [s]",
    heading: str | None = None,
) -> plt.Figure:
    clusters = sorted(data["clusters"].unique())
    columns = min(3, len(clusters))
    rows = math.ceil(len(clusters) / columns)
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(5.1 * columns, 3.8 * rows),
        sharey=y == "runtime_s",
        squeeze=False,
        constrained_layout=True,
    )
    ordered_resolutions = (
        data[["modeled_hours", "periods", "hours_per_period"]]
        .drop_duplicates()
        .sort_values("modeled_hours")
    )
    ticks = ordered_resolutions["modeled_hours"].tolist()
    tick_labels = [
        f"{periods}c{hours}"
        for _, periods, hours in ordered_resolutions.itertuples(index=False, name=None)
    ]
    for ax, cluster in zip(axes.flat, clusters):
        subset = data.loc[data["clusters"] == cluster]
        draw_series(ax, subset, "modeled_hours", solvers, colors, y)
        panel_log_scale = log_scale or (y == "objective" and needs_log_scale(subset[y]))
        style_axis(ax, panel_log_scale)
        ax.set_title(f"{cluster:g} nodes", loc="left", fontweight="semibold")
        ax.set_xlabel("Time resolution (periods × hours)")
        ax.set_xticks(ticks, tick_labels, rotation=35, ha="right")
    for ax in axes[:, 0]:
        ax.set_ylabel(y_label)
    for ax in axes.flat[len(clusters) :]:
        ax.set_visible(False)
    figure.suptitle(f"{heading or title}\nModeled time", fontweight="bold", fontsize=15)
    figure_legend(figure, data, solvers, colors)
    return figure


def save_figure(
    figure: plt.Figure, output_stem: Path, formats: list[str], dpi: int
) -> None:
    for extension in formats:
        path = output_stem.with_suffix(f".{extension}")
        figure.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        logger.info("Wrote %s", path)
    plt.close(figure)


def plot_case(
    data: pd.DataFrame,
    args: argparse.Namespace,
    output_dir: Path,
    colors: dict[str, object],
    stem_suffix: str,
    title_suffix: str,
) -> None:
    """Create all benchmark figures for one case variant."""
    solvers = sorted(data["solver_options"].unique())
    title = f"{args.title} — {title_suffix}" if title_suffix else args.title
    figures = (
        ("node_scaling", make_node_scaling_figure),
        ("time_scaling", make_time_scaling_figure),
    )
    for dimension, make_figure in figures:
        figure = make_figure(data, solvers, colors, title, args.log_scale)
        save_figure(
            figure,
            output_dir / f"benchmark_{dimension}{stem_suffix}",
            args.formats,
            args.dpi,
        )

    memory_data = data.loc[data["memory_gib"].notna() & (data["memory_gib"] > 0)].copy()
    if memory_data.empty:
        logger.warning("No positive max_rss values found; skipping memory plots.")
    else:
        # Timing lines are only relevant to the runtime figures.
        memory_data[list(SMSPP_TIMING_STYLES)] = float("nan")
        memory_solvers = sorted(memory_data["solver_options"].unique())
        heading = f"Peak resident memory — {title_suffix}" if title_suffix else "Peak resident memory"
        for dimension, make_figure in figures:
            figure = make_figure(
                memory_data, memory_solvers, colors, title, False,
                y="memory_gib", y_label="Peak resident memory [GiB]", heading=heading,
            )
            save_figure(
                figure, output_dir / f"benchmark_memory_{dimension}{stem_suffix}",
                args.formats, args.dpi,
            )

    valid_objective = data["objective"].notna() & (data["objective"] > 0)
    if not valid_objective.all():
        logger.warning(
            "Ignoring %d row(s) without a positive objective value%s.",
            int((~valid_objective).sum()),
            f" for case variant {title_suffix}" if title_suffix else "",
        )
    objective_data = data.loc[valid_objective].copy()
    if objective_data.empty:
        logger.warning("No numeric objective values found; skipping objective plots.")
        return
    objective_solvers = sorted(objective_data["solver_options"].unique())
    heading = (
        f"Objective function — {title_suffix}" if title_suffix else "Objective function"
    )
    for dimension, make_figure in figures:
        figure = make_figure(
            objective_data,
            objective_solvers,
            colors,
            title,
            False,
            y="objective",
            y_label="Objective value",
            heading=heading,
        )
        save_figure(
            figure,
            output_dir / f"benchmark_objective_{dimension}{stem_suffix}",
            args.formats,
            args.dpi,
        )


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir or args.input.parent / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    data = prepare_data(args.input)
    all_solvers = sorted(data["solver_options"].unique())
    colors = solver_colors(all_solvers)
    case_dimensions = (
        data[["case_variant", "scenarios"]]
        .drop_duplicates()
        .sort_values(["case_variant", "scenarios"])
    )
    qualify_stems = (
        len(case_dimensions) > 1
        or any(case_dimensions["case_variant"])
        or any(case_dimensions["scenarios"] != 1)
    )
    for variant, scenarios in case_dimensions.itertuples(index=False, name=None):
        variant_data = data.loc[
            (data["case_variant"] == variant) & (data["scenarios"] == scenarios)
        ].copy()
        scenario_label = f"{scenarios:g} scenario" + ("" if scenarios == 1 else "s")
        title_parts = [part for part in (variant, scenario_label) if part]
        title_suffix = " — ".join(title_parts)
        stem_parts = [variant_slug(variant)] if variant else []
        if scenarios != 1:
            stem_parts.append(f"scenarios_{scenarios:g}")
        stem_suffix = "_" + "_".join(stem_parts) if qualify_stems else ""
        plot_case(variant_data, args, output_dir, colors, stem_suffix, title_suffix)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
