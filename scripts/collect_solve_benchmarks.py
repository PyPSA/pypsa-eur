# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Collect solve benchmarks and objectives below results prefixes or selected case folders.

python scripts/collect_solve_benchmarks.py \
    --prefix dispatch-power-IT \
    --output results/solve_benchmarks.csv

python scripts/collect_solve_benchmarks.py \
    --folder dispatch-power-IT/th_24c24 dispatch-power-IT/th_48c24 \
    --output results/selected_solve_benchmarks.csv
"""

import argparse
import json
import logging
import re
from pathlib import Path

import pandas as pd
import xarray as xr
import yaml

logger = logging.getLogger(__name__)
SMSPP_BENCHMARK_SUFFIX = ".smspp"

STATUS_RE = re.compile(
    r"Solving status ['\"](?P<status>[^'\"]+)['\"] "
    r"with termination condition ['\"](?P<condition>[^'\"]+)['\"]",
    re.IGNORECASE,
)
TERMINATION_RE = re.compile(
    r"Termination condition:\s*(?P<condition>[^\r\n]+)", re.IGNORECASE
)
TIME_LIMIT_RE = re.compile(
    r"(?:TimeLimit|time[_ ]limit)\s*[:=]\s*(?P<seconds>\d+(?:\.\d+)?)",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument(
        "--prefix",
        dest="prefixes",
        nargs="+",
        help="Prefixes below --results-root to scan recursively.",
    )
    selection.add_argument(
        "--folder",
        dest="folders",
        nargs="+",
        type=Path,
        help=(
            "Specific case folders to scan recursively. Relative paths are resolved "
            "below --results-root; paths beginning with the results-root name and "
            "absolute paths are also accepted."
        ),
    )
    parser.add_argument("--results-root", default=Path("results"), type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def parse_benchmark_name(stem: str, solver: str | None = None) -> dict:
    case_name = stem.removesuffix(f"_{solver}") if solver else stem
    metadata = {
        "benchmark": stem,
        "case_result": case_name,
        "solver": solver,
        "network_type": pd.NA,
        "clusters": pd.NA,
        "opts": pd.NA,
        "sector_opts": pd.NA,
        "planning_horizons": pd.NA,
    }
    if not case_name.startswith("base_s_"):
        return metadata

    body = case_name.removeprefix("base_s_")
    if "_elec_" in body:
        clusters, opts = body.split("_elec_", 1)
        metadata.update(
            {"network_type": "electricity", "clusters": clusters, "opts": opts}
        )
        return metadata

    parts = body.split("_")
    if len(parts) >= 4 and parts[-3:] == ["brownfield", "all", "years"]:
        metadata.update(
            {
                "network_type": "sector",
                "clusters": parts[0],
                "opts": parts[1] if len(parts) > 1 else "",
                "sector_opts": "_".join(parts[2:-3]),
                "planning_horizons": "brownfield_all_years",
            }
        )
    elif len(parts) >= 4:
        metadata.update(
            {
                "network_type": "sector",
                "clusters": parts[0],
                "opts": parts[1],
                "sector_opts": "_".join(parts[2:-1]),
                "planning_horizons": parts[-1],
            }
        )
    return metadata


def benchmark_files(root: Path) -> list[Path]:
    if not root.exists():
        logger.warning("Results prefix does not exist: %s", root)
        return []
    paths: list[Path] = []
    for rule_dir in ["solve_network", "solve_sector_network"]:
        paths.extend(root.glob(f"**/benchmarks/{rule_dir}/*"))
    # Snakemake creates empty placeholders for unused log files on success.
    failures = {
        path
        for path in paths
        if path.is_file() and path.name.endswith(".failed") and path.stat().st_size
    }
    return sorted(
        path
        for path in paths
        if path.is_file()
        and not path.name.endswith(SMSPP_BENCHMARK_SUFFIX)
        and (
            path in failures
            if path.name.endswith(".failed")
            else path.with_name(path.name + ".failed") not in failures
        )
    )


def read_benchmark(path: Path) -> dict:
    benchmark = pd.read_csv(path, sep="\t").iloc[-1].to_dict()
    return {key: pd.NA if value == "NA" else value for key, value in benchmark.items()}


def read_solve_outcome(path: Path, configured_time_limit: object = pd.NA) -> dict:
    """Read the final solve outcome and configured time limit from its Python log."""
    log_path = (
        result_directory(path) / "logs" / path.parent.name / f"{path.name}_python.log"
    )
    if not log_path.exists():
        log_path = result_directory(path) / "logs" / f"{path.name}_python.log"
    metadata = {
        "solve_status": pd.NA,
        "termination_condition": pd.NA,
        "time_limit_s": configured_time_limit,
    }
    if not log_path.exists():
        logger.warning("No Python solve log found for %s", path)
        return metadata

    text = log_path.read_text(errors="replace")
    status_matches = list(STATUS_RE.finditer(text))
    termination_matches = list(TERMINATION_RE.finditer(text))
    limit_matches = list(TIME_LIMIT_RE.finditer(text))
    if status_matches:
        match = status_matches[-1]
        metadata["solve_status"] = match["status"].strip()
        metadata["termination_condition"] = match["condition"].strip()
    elif termination_matches:
        metadata["termination_condition"] = termination_matches[-1]["condition"].strip()
    if limit_matches:
        metadata["time_limit_s"] = float(limit_matches[-1]["seconds"])
    return metadata


def result_directory(path: Path) -> Path:
    """Return the directory containing the benchmarks and networks folders."""
    benchmark_index = path.parts.index("benchmarks")
    return Path(*path.parts[:benchmark_index])


def read_solver_metadata(config: dict) -> dict:
    """Read solver identity and limits from a saved run configuration."""
    metadata = {}
    solving = config.get("solving", {})
    solver = solving.get("solver", {})
    option_name = solver.get("options")
    option_values = solving.get("solver_options", {}).get(option_name, {})
    time_limit = option_values.get("TimeLimit", pd.NA)
    configfile = option_values.get("configfile")
    if pd.isna(time_limit) and configfile and Path(configfile).exists():
        config_text = Path(configfile).read_text(errors="replace")
        match = re.search(r"\bdblMaxTime\s+(\d+(?:\.\d+)?)", config_text)
        if match:
            time_limit = float(match[1])
    metadata.update(
        {
            "solver": config.get("wildcards", {}).get("solver"),
            "solver_name": solver.get("name", pd.NA),
            "solver_options": option_name,
            "time_limit_s": time_limit,
        }
    )
    return metadata


def read_network_metadata(path: Path) -> dict:
    """Read scalar result metadata without loading network component tables."""
    if not path.exists():
        logger.warning("No solved network found for %s", path)
        return {}

    with xr.open_dataset(path) as ds:
        metadata = {
            "network_file": path.as_posix(),
            "objective": ds.attrs.get("network__objective", pd.NA),
            "objective_constant": ds.attrs.get("network__objective_constant", pd.NA),
            "scenarios": int(ds.sizes.get("scenario", 1)),
        }
        meta = ds.attrs.get("meta")

    if not meta:
        return metadata
    try:
        config = json.loads(meta)
        metadata.update(read_solver_metadata(config))
    except (TypeError, json.JSONDecodeError):
        logger.warning("Could not parse network metadata in %s", path)
    return metadata


def collect_benchmarks(roots: list[Path], case_folders: bool = False) -> pd.DataFrame:
    rows = []
    for root in roots:
        for benchmark_path in benchmark_files(root):
            failed = benchmark_path.name.endswith(".failed")
            path = benchmark_path.with_name(benchmark_path.name.removesuffix(".failed"))
            result_dir = result_directory(path)
            network_path = result_dir / "networks" / f"{path.name}.nc"
            config_path = result_dir / "configs" / f"config.{path.name}.yaml"
            network_metadata = {} if failed else read_network_metadata(network_path)
            if failed:
                # Failed runs may not have exported their own config. Reuse the
                # same solver variant from another run within this results prefix.
                solver_match = re.search(r"__(?P<solver>[^/]+)$", path.name)
                solver_id = solver_match["solver"] if solver_match else None
                if not config_path.exists() and solver_id:
                    search_root = root.parent if case_folders else root
                    candidates = sorted(
                        search_root.glob(f"**/configs/config.*__{solver_id}.yaml")
                    )
                    if candidates:
                        config_path = candidates[0]
                if config_path.exists():
                    network_metadata = read_solver_metadata(
                        yaml.safe_load(config_path.read_text())
                    )
                if solver_id:
                    network_metadata["solver"] = solver_id
            metadata_root = root.parent if case_folders else root
            row = {
                "prefix": metadata_root.as_posix(),
                "case": result_dir.relative_to(metadata_root).as_posix(),
                "benchmark_rule": path.parent.name,
                "benchmark_file": benchmark_path.as_posix(),
                "failed": failed,
            }
            row.update(network_metadata)
            row.update(parse_benchmark_name(path.name, row.get("solver")))
            row.update(read_benchmark(benchmark_path))
            smspp_benchmark = path.with_name(path.name + SMSPP_BENCHMARK_SUFFIX)
            if not failed and smspp_benchmark.exists():
                smspp_timings = read_benchmark(smspp_benchmark)
                row["smspp_optimization_s"] = smspp_timings["s"]
                row["smspp_computational_s"] = smspp_timings.get(
                    "computational_time", pd.NA
                )
            row.update(read_solve_outcome(path, row.get("time_limit_s", pd.NA)))
            rows.append(row)

    leading_columns = [
        "prefix",
        "case",
        "benchmark_rule",
        "case_result",
        "solver",
        "solver_name",
        "solver_options",
        "network_type",
        "clusters",
        "opts",
        "sector_opts",
        "planning_horizons",
        "scenarios",
        "benchmark_file",
        "network_file",
        "objective",
        "objective_constant",
        "solve_status",
        "termination_condition",
        "time_limit_s",
    ]
    if not rows:
        return pd.DataFrame(columns=leading_columns)

    df = pd.DataFrame(rows)
    remaining_columns = [
        column for column in df.columns if column not in leading_columns
    ]
    return df.reindex(columns=leading_columns + remaining_columns).sort_values(
        ["prefix", "case", "benchmark_rule", "case_result", "solver"],
        na_position="last",
    )


def resolve_folder(folder: Path, results_root: Path) -> Path:
    """Resolve a selected folder while accepting results/... paths."""
    if folder.is_absolute():
        return folder
    if folder.parts and folder.parts[0] == results_root.name:
        return folder
    return results_root / folder


def main() -> None:
    args = parse_args()
    if args.prefixes:
        roots = [args.results_root / prefix.strip("/") for prefix in args.prefixes]
        case_folders = False
    else:
        roots = [resolve_folder(folder, args.results_root) for folder in args.folders]
        case_folders = True
    df = collect_benchmarks(roots, case_folders=case_folders)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    logger.info("Collected %s solve benchmark rows into %s", len(df), args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
