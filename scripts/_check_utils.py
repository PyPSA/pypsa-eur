"""
Shared utilities for carbon dioxide removal pipeline check scripts.

Loads run parameters from the saved pypsa-eur config so check scripts need
only a single input: the run name.

Usage in any check script:

    from _check_utils import parse_check_args, load_check_params

    args = parse_check_args()        # positional run_name, optional overrides
    p    = load_check_params(args)   # returns dict with RDIR, WC, paths, cfg, …

Run from the pypsa-eur root:

    python scripts/check_EW_pipeline.py EW_2050
    python scripts/check_EW_pipeline.py --config results/EW_2050/EW_2050/configs/config.base_s_90__168h_2050.yaml
"""

import argparse
import sys
from pathlib import Path

import yaml

BASE_DIR = Path(".")   # check scripts are run from the pypsa-eur root


# ── Saved-config discovery ─────────────────────────────────────────────────────

def find_saved_config(run_name: str, base_dir: Path) -> Path:
    """Find the config saved by pypsa-eur under results/<run_name>/**/configs/."""
    pattern = f"results/{run_name}/**/configs/config.*.yaml"
    hits = sorted(base_dir.glob(pattern))
    if not hits:
        sys.exit(
            f"ERROR: no saved config found for run '{run_name}'.\n"
            f"  Looked for: {base_dir / pattern}\n"
            f"  Use --config <path> to point directly at the config file."
        )
    if len(hits) > 1:
        print(
            "WARNING: multiple saved configs found; using first:\n  "
            + "\n  ".join(str(h) for h in hits)
        )
    return hits[0]


def _load_yaml(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f) or {}


def load_config(args: argparse.Namespace, base_dir: Path = BASE_DIR) -> tuple:
    """
    Load the saved run config. Returns (cfg_dict, config_path).

    Priority:
      1. --config <path>    → load that file directly
      2. run_name (positional) → auto-discover under results/<run_name>/**/configs/
    """
    if getattr(args, "config", None):
        config_path = Path(args.config)
        if not config_path.exists():
            sys.exit(f"ERROR: config file not found: {config_path}")
    elif getattr(args, "run_name", None):
        config_path = find_saved_config(args.run_name, base_dir)
    else:
        sys.exit("ERROR: provide a run_name or --config <path>.")

    cfg = _load_yaml(config_path)
    print(f"Config: {config_path}")
    return cfg, config_path


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_check_args(extra_args=None) -> argparse.Namespace:
    """
    Common CLI parser for all carbon dioxide removal check scripts.

    Positional:
      run_name     Run directory name (e.g. rock_weathering_2050). Script finds the saved
                   config under results/<run_name>/**/configs/ automatically.

    Optional overrides (take precedence over config values):
      --config     Direct path to a saved config YAML (skips auto-discovery).
      --base-dir   pypsa-eur root directory (default: current directory).
      --run-name   Override run name from config.
      --clusters   Override cluster count.
      --opts       Override opts wildcard.
      --sector-opts Override sector opts wildcard.
      --horizon    Override planning horizon year.

    extra_args: list of ([flags], kwargs) for script-specific arguments.
    """
    p = argparse.ArgumentParser(
        description="Carbon dioxide removal pipeline check — reads wildcards from the saved run config.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "run_name", nargs="?",
        help="Run name (e.g. EW_2050). Auto-discovers saved config under results/<run_name>/.",
    )
    p.add_argument(
        "--config", "-c", default=None, metavar="PATH",
        help="Direct path to a saved config YAML (overrides run_name auto-discovery).",
    )
    p.add_argument(
        "--base-dir", default=".", metavar="DIR",
        help="pypsa-eur root directory (default: current directory).",
    )
    p.add_argument("--run-name", default=None, metavar="NAME",
                   help="Override run name (= results/resources sub-directory).")
    p.add_argument("--clusters", default=None, metavar="N",
                   help="Override cluster count (e.g. 90).")
    p.add_argument("--opts", default=None, metavar="OPTS",
                   help="Override opts wildcard.")
    p.add_argument("--sector-opts", default=None, metavar="OPTS",
                   help="Override sector opts wildcard (e.g. 168h).")
    p.add_argument("--horizon", default=None, metavar="YEAR",
                   help="Override planning horizon year (e.g. 2050).")

    if extra_args:
        for flags, kwargs in extra_args:
            p.add_argument(*flags, **kwargs)

    args = p.parse_args()
    if not args.run_name and not args.config:
        p.print_help()
        sys.exit(1)
    return args


# ── Parameter extraction ───────────────────────────────────────────────────────

def load_check_params(args: argparse.Namespace) -> dict:
    """
    Load the saved config and return a dict with all derived parameters:

      RDIR             – run directory name
      CLUSTERS         – cluster count string (e.g. "90")
      OPTS             – opts wildcard (e.g. "")
      SECTOR_OPTS      – sector opts wildcard (e.g. "168h")
      PLANNING_HORIZON – planning horizon string (e.g. "2050")
      WC               – full wildcard stem (e.g. "base_s_90__168h_2050")
      SHARED_RES       – shared-resources Path (= RES_RUN if no shared policy)
      BASE_DIR         – Path to pypsa-eur root
      RES              – Path to resources/
      RES_RUN          – Path to resources/<RDIR>/
      RESULTS          – Path to results/<RDIR>/
      cfg              – full config dict (for script-specific lookups)
    """
    base_dir = Path(getattr(args, "base_dir", "."))
    cfg, _ = load_config(args, base_dir)

    run_cfg      = cfg.get("run", {})
    scenario_cfg = cfg.get("scenario", {})

    # CLI --run-name > config run.name > positional run_name
    RDIR = (
        getattr(args, "run_name_override", None)   # --run-name flag (argparse stores as run_name_override below)
        or run_cfg.get("name")
        or getattr(args, "run_name", None)
        or "run"
    )
    # Note: argparse stores --run-name as args.run_name which clashes with the
    # positional. We use dest="run_name_cli" to separate them.
    # In practice, parse_check_args stores --run-name in args.run_name (the flag)
    # and the positional in args.run_name too — last one wins in argparse.
    # Simpler: just use config value as primary, CLI flags as overrides.
    RDIR = run_cfg.get("name") or getattr(args, "run_name", None) or "run"
    if getattr(args, "run_name", None) and not run_cfg.get("name"):
        RDIR = args.run_name

    CLUSTERS = str(
        args.clusters
        or (scenario_cfg.get("clusters") or [90])[0]
    )
    OPTS = (
        args.opts if args.opts is not None
        else str((scenario_cfg.get("opts") or [""])[0])
    )
    SECTOR_OPTS = str(
        args.sector_opts
        or (scenario_cfg.get("sector_opts") or ["168h"])[0]
    )
    PLANNING_HORIZON = str(
        args.horizon
        or (scenario_cfg.get("planning_horizons") or [2050])[-1]
    )

    WC = f"base_s_{CLUSTERS}_{OPTS}_{SECTOR_OPTS}_{PLANNING_HORIZON}"

    RES     = base_dir / "resources"
    RES_RUN = base_dir / "resources" / RDIR
    RESULTS = base_dir / "results"   / RDIR

    # shared resources: if run.shared_resources.policy is a string, resources live there
    shared_policy = run_cfg.get("shared_resources", {}).get("policy", False)
    SHARED_RES = (base_dir / "resources" / shared_policy) if shared_policy else RES_RUN

    print(f"Run:    {RDIR}  |  WC: {WC}")

    return dict(
        RDIR=RDIR,
        CLUSTERS=CLUSTERS,
        OPTS=OPTS,
        SECTOR_OPTS=SECTOR_OPTS,
        PLANNING_HORIZON=PLANNING_HORIZON,
        WC=WC,
        BASE_DIR=base_dir,
        RES=RES,
        RES_RUN=RES_RUN,
        RESULTS=RESULTS,
        SHARED_RES=SHARED_RES,
        cfg=cfg,
    )
