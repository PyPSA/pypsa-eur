"""
Patch existing cdr_credit_accounting CSVs that show method=solver_unavailable
by recomputing credited values from capture_proxy (same fallback logic as the
fixed solve_network.py).

Usage:
    pixi run python scripts/repatch_cdr_accounting.py results/supply_curve_v8_medium/
    pixi run python scripts/repatch_cdr_accounting.py results/supply_curve_v7_medium/
    pixi run python scripts/repatch_cdr_accounting.py results/  # all subdirs
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd


def patch_csv(csv_path: Path, dry_run: bool = False) -> bool:
    df = pd.read_csv(csv_path)
    changed = False

    for i, row in df.iterrows():
        if row.get("method") != "solver_unavailable":
            continue

        proxy_dac = row.get("capture_proxy_dac_mtco2_per_yr", 0.0)
        proxy_biogenic = row.get("capture_proxy_biogenic_mtco2_per_yr", 0.0)
        proxy_total = proxy_dac + proxy_biogenic
        credit_limit = row.get("credit_limit_mtco2_per_yr", np.nan)
        tolerance = 1e-3

        if proxy_total <= 0:
            continue

        if not np.isnan(credit_limit) and proxy_total > credit_limit:
            scale = credit_limit / proxy_total
            credited_dac = proxy_dac * scale
            credited_biogenic = proxy_biogenic * scale
        else:
            credited_dac = proxy_dac
            credited_biogenic = proxy_biogenic

        credited_total = credited_dac + credited_biogenic
        capture_total = row.get("capture_proxy_total_mtco2_per_yr", proxy_total)

        df.at[i, "method"] = "capture_proxy_fallback"
        df.at[i, "accounting_error"] = ""
        df.at[i, "credited_dac_mtco2_per_yr"] = credited_dac
        df.at[i, "credited_biogenic_mtco2_per_yr"] = credited_biogenic
        df.at[i, "credited_total_mtco2_per_yr"] = credited_total
        df.at[i, "valid_credited_within_limit"] = (
            True if np.isnan(credit_limit) else credited_total <= credit_limit + tolerance
        )
        df.at[i, "valid_credited_within_capture_proxy"] = (
            credited_total <= capture_total + tolerance
        )
        changed = True

    if changed and not dry_run:
        df.to_csv(csv_path, index=False)

    return changed


def main():
    roots = sys.argv[1:] if len(sys.argv) > 1 else ["."]
    dry_run = "--dry-run" in roots
    roots = [r for r in roots if r != "--dry-run"]

    patched = 0
    for root in roots:
        for csv_path in Path(root).rglob("cdr_credit_accounting_*.csv"):
            if patch_csv(csv_path, dry_run=dry_run):
                print(f"{'[DRY] ' if dry_run else ''}Patched {csv_path}")
                patched += 1

    print(f"\n{'Would patch' if dry_run else 'Patched'} {patched} CSV(s).")


if __name__ == "__main__":
    main()
