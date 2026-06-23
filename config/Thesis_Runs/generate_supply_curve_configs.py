"""
Generate supply curve config folders for low / medium / high CDR technology cost scenarios.

Usage:
    python generate_supply_curve_configs.py              # generates all 3 scenarios
    python generate_supply_curve_configs.py low high     # generates only named scenarios

Reads configs from supply_curve_test/ as the medium template.
Writes to supply_curve_low/ and supply_curve_high/ (medium = supply_curve_test, unchanged).

Cost CSVs expected at:
    data/costs/custom_costs_medium.csv  ← medium
    data/costs/custom_costs_low.csv     ← low scenario
    data/costs/custom_costs_high.csv    ← high scenario
"""

import re
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent  # pypsa-eur-thesis/

SCENARIOS = {
    "medium": {
        "folder": "supply_curve_test",
        "cost_fn": "data/costs/custom_costs_medium.csv",
        "name_infix": "t",          # e.g. S03t-cdr-150eur-168seg (existing naming)
        "skip_generate": True,      # medium already exists — don't overwrite
    },
    "low": {
        "folder": "supply_curve_low",
        "cost_fn": "data/costs/custom_costs_low.csv",
        "name_infix": "low",        # e.g. S03low-cdr-150eur-168seg
        "skip_generate": False,
    },
    "high": {
        "folder": "supply_curve_high",
        "cost_fn": "data/costs/custom_costs_high.csv",
        "name_infix": "hi",         # e.g. S03hi-cdr-150eur-168seg
        "skip_generate": False,
    },
}

TEMPLATE_FOLDER = Path(__file__).parent / "supply_curve_test"


def generate_scenario(scenario_key: str) -> None:
    spec = SCENARIOS[scenario_key]

    if spec["skip_generate"]:
        print(f"  Skipping {scenario_key} (supply_curve_test is the medium reference).")
        return

    # Validate cost file exists
    cost_fn = REPO_ROOT / spec["cost_fn"]
    if not cost_fn.exists():
        print(f"  WARNING: {cost_fn} not found — create it before running. Skipping.")
        return

    target_folder = Path(__file__).parent / spec["folder"]
    target_folder.mkdir(exist_ok=True)

    template_files = sorted(TEMPLATE_FOLDER.glob("config.S*.yaml"))
    if not template_files:
        print(f"  ERROR: No config.S*.yaml found in {TEMPLATE_FOLDER}")
        return

    infix = spec["name_infix"]
    cost_fn_rel = spec["cost_fn"]  # relative path from repo root, used in YAML

    for src in template_files:
        text = src.read_text()

        # Replace run.name: "S03t-cdr-..." → "S03{infix}-cdr-..."
        # The medium uses infix "t" so we replace that pattern
        text = re.sub(
            r'(name:\s*"S\d+)t(-cdr-[^"]+)',
            rf'\g<1>{infix}\2',
            text,
        )

        # Inject or replace costs.custom_cost_fn
        # The costs: block exists but may not have custom_cost_fn yet
        if "custom_cost_fn:" in text:
            text = re.sub(
                r"(custom_cost_fn:\s*).*",
                rf"\g<1>{cost_fn_rel}",
                text,
            )
        else:
            # Insert custom_cost_fn after the "costs:" line (before emission_prices or next key)
            text = re.sub(
                r"(^costs:\s*\n)",
                rf"\g<1>  custom_cost_fn: {cost_fn_rel}\n",
                text,
                flags=re.MULTILINE,
            )

        # Update header comment
        scenario_label = scenario_key.upper()
        text = re.sub(
            r"# Supply curve test scenario",
            f"# Supply curve {scenario_label} cost scenario",
            text,
        )

        dst = target_folder / src.name
        dst.write_text(text)
        print(f"  wrote {dst.relative_to(REPO_ROOT)}")

    print(f"  {scenario_key}: {len(template_files)} configs written to {target_folder.relative_to(REPO_ROOT)}/")


def main() -> None:
    requested = sys.argv[1:] if len(sys.argv) > 1 else list(SCENARIOS.keys())

    unknown = [s for s in requested if s not in SCENARIOS]
    if unknown:
        print(f"Unknown scenarios: {unknown}. Valid: {list(SCENARIOS.keys())}")
        sys.exit(1)

    for scenario in requested:
        print(f"\nGenerating '{scenario}'...")
        generate_scenario(scenario)


if __name__ == "__main__":
    main()
