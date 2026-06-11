"""Regenerate cdr_credit_accounting_*.csv from solved networks — NO re-solving.

The CDR credit accounting written at solve time always fell back to a capture
proxy that could exceed physical geological sequestration (e.g. credited 626 Mt
vs a 600 Mt storage cap), because the solver-side attribution variable
(CO2AnnualCDRSeq) is underdetermined over the fungible "co2 stored" pool and its
post-solve solution is not recoverable.  `export_cdr_credit_accounting` in
solve_network.py has since been changed to compute credited CDR by deterministic
*sequestration proration* (allocate physical sequestration to capture origins
pro-rata by captured volume; credit only the DAC + biogenic shares).

Every input that method needs — capture-by-origin and physical sequestration —
is present in the solved .nc, so we can simply reload each network and re-export
the accounting CSV in place.  This script does exactly that for the v12 results.

Run:  pixi run python scripts/regenerate_cdr_accounting.py
      pixi run python scripts/regenerate_cdr_accounting.py --dry-run
      pixi run python scripts/regenerate_cdr_accounting.py --root results/supply_curve_v10
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import pypsa
import yaml

REPO = Path(__file__).resolve().parents[1]
# solve_network.py imports `scripts._benchmark`, so the repo ROOT must be on the
# path (so that `scripts` resolves as a package), not just the scripts/ dir.
sys.path.insert(0, str(REPO))
from scripts.solve_network import export_cdr_credit_accounting  # noqa: E402

YEARS = (2030, 2040, 2050)

# n.params keys read by export_cdr_credit_accounting, mapped from the saved
# config's `sector` section (mirrors rule solve_sector_network_myopic.params).
def _hydrate_params(sector: dict) -> dict:
    return {
        # NB: the solve rule binds n.params["cdr_credit_limit"] to the *by_year*
        # dict; export reads cdr_credit_limit_by_year first, then cdr_credit_limit.
        "cdr_credit_limit": sector.get("cdr_credit_limit_by_year"),
        "cdr_credit_scope": sector.get("cdr_credit_scope", []),
        "cdr_credit_timing": sector.get("cdr_credit_timing", "capture"),
        "cdr_credit_price": sector.get("cdr_credit_price", 0.0),
        "cdr_credit_prices_by_scope": sector.get("cdr_credit_prices_by_scope", {}),
        "cdr_credit_standalone": sector.get("cdr_credit_standalone", False),
    }


def _old_credited(csv_path: Path) -> tuple[float, float, str]:
    if not csv_path.exists():
        return float("nan"), float("nan"), "—"
    df = pd.read_csv(csv_path)
    if df.empty:
        return float("nan"), float("nan"), "—"
    row = df.iloc[0]
    return (
        float(row.get("credited_total_mtco2_per_yr", float("nan"))),
        float(row.get("physical_sequestration_mtco2_per_yr", float("nan"))),
        str(row.get("method", "—")),
    )


def regenerate_one(net_path: Path, csv_path: Path, cfg_path: Path, year: int,
                   dry_run: bool) -> dict:
    n = pypsa.Network(str(net_path))
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)
    n.params = _hydrate_params(cfg.get("sector", {}))

    old_cred, old_seq, old_method = _old_credited(csv_path)
    if not dry_run:
        export_cdr_credit_accounting(n, csv_path, str(year))
    new_cred, new_seq, new_method = _old_credited(csv_path)  # re-read after write

    return {
        "scenario": net_path.parts[-3],
        "year": year,
        "old_method": old_method,
        "old_credited": old_cred,
        "new_credited": (old_cred if dry_run else new_cred),
        "physical_seq": new_seq if new_seq == new_seq else old_seq,
        "new_method": (old_method if dry_run else new_method),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="results/supply_curve_v12",
                    help="results root to walk (default: v12)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report before/after without writing CSVs")
    args = ap.parse_args()

    root = (REPO / args.root).resolve()
    net_paths = sorted(root.glob("**/networks/base_s_*_*seg_*.nc"))
    if not net_paths:
        raise SystemExit(f"No solved networks found under {root}")

    records = []
    for net_path in net_paths:
        scen_dir = net_path.parents[1]  # .../<scenario>/
        name = net_path.name  # base_s_96__168seg_<year>.nc
        try:
            year = int(name.rsplit("_", 1)[1].split(".")[0])
        except (IndexError, ValueError):
            print(f"skip (cannot parse year): {net_path}")
            continue
        if year not in YEARS:
            continue
        stem = name[: -len(".nc")]
        csv_path = scen_dir / "csvs" / "individual" / f"cdr_credit_accounting_{stem.replace('base_', '')}.csv"
        cfg_path = scen_dir / "configs" / f"config.{stem}.yaml"
        if not cfg_path.exists():
            print(f"skip (no config): {cfg_path}")
            continue
        try:
            rec = regenerate_one(net_path, csv_path, cfg_path, year, args.dry_run)
            records.append(rec)
            flag = "" if (rec["new_credited"] <= (rec["physical_seq"] or 0) + 1e-3) else "  <-- STILL OVER SEQ!"
            print(f"{rec['scenario']:<40} {year}  "
                  f"{rec['old_credited']:7.1f} -> {rec['new_credited']:7.1f} Mt  "
                  f"(seq {rec['physical_seq']:6.1f})  [{rec['new_method']}]{flag}")
        except Exception as exc:  # noqa: BLE001
            print(f"FAILED {net_path}: {exc}")

    if records:
        df = pd.DataFrame(records)
        over = df[df["new_credited"] > df["physical_seq"] + 1e-3]
        print(f"\n{len(df)} networks processed"
              + (" (dry-run, nothing written)" if args.dry_run else ""))
        print(f"credited > physical_seq violations remaining: {len(over)}")
        if len(over):
            print(over.to_string(index=False))


if __name__ == "__main__":
    main()
