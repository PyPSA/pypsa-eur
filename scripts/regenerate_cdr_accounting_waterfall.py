"""Write a COMPLEMENTARY waterfall CDR-credit accounting CSV — NO re-solving.

Companion to ``regenerate_cdr_accounting.py`` (which produces the *pro-rata*
``cdr_credit_accounting_*.csv``).  This script leaves that file untouched and
writes a second file next to it, ``cdr_credit_accounting_waterfall_*.csv``, using
the waterfall attribution rule instead.

Two conventions, both post-processing overlays on the *same* solved networks:

  pro-rata (physical-pool)   physical sequestration is split across capture
                             origins in proportion to captured volume; clean
                             (DAC + biogenic) origins get their captured share of
                             storage.  Conservative lower bound on credited CDR.

  waterfall (optimiser-      eligible DAC + biogenic tonnes get FIRST claim on
  consistent)                geological storage; fossil CO2 is stored only with
                             the residual capacity and earns no credit.  Matches
                             the solver's economic incentive (eligible origins
                             out-bid fossil for scarce storage slots).

Neither is physically observable once the "co2 stored" pool is commingled; both
are accounting conventions.  Waterfall is the policy-/optimiser-consistent rule,
pro-rata is the conservative sensitivity.  Both satisfy, by construction,
    credited_total <= clean_capture      and
    credited_total <= physical_sequestration.

The waterfall file uses the SAME column schema as the pro-rata file (so the same
``credited_cdr(acc, price)`` reader works on either), plus a few waterfall
diagnostics and the pro-rata credited value for side-by-side comparison.

Run:  pixi run python scripts/regenerate_cdr_accounting_waterfall.py
      pixi run python scripts/regenerate_cdr_accounting_waterfall.py --dry-run
      pixi run python scripts/regenerate_cdr_accounting_waterfall.py --root results/supply_curve_v12
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import yaml

REPO = Path(__file__).resolve().parents[1]
# solve_network.py imports `scripts._benchmark`, so the repo ROOT must be on the
# path (so `scripts` resolves as a package), not just the scripts/ dir.
sys.path.insert(0, str(REPO))
from scripts._helpers import get  # noqa: E402
from scripts.solve_network import (  # noqa: E402
    _annual_capture_proxy_by_origin,
    _annual_physical_sequestration,
    _eligible_cdr_scopes,
    _get_cdr_credit_prices_for_period,
)

YEARS = (2030, 2040, 2050)
TOL_MT = 1e-3


def _hydrate_params(sector: dict) -> dict:
    """Mirror rule solve_sector_network_myopic.params (see regenerate_cdr_accounting.py)."""
    return {
        "cdr_credit_limit": sector.get("cdr_credit_limit_by_year"),
        "cdr_credit_limit_by_year": sector.get("cdr_credit_limit_by_year"),
        "cdr_credit_scope": sector.get("cdr_credit_scope", []),
        "cdr_credit_timing": sector.get("cdr_credit_timing", "capture"),
        "cdr_credit_price": sector.get("cdr_credit_price", 0.0),
        "cdr_credit_prices_by_scope": sector.get("cdr_credit_prices_by_scope", {}),
        "cdr_credit_standalone": sector.get("cdr_credit_standalone", False),
    }


def _waterfall_credited(cap_dac, cap_biogenic, cap_fossil, physical_seq_t, credit_limit_mt):
    """Waterfall attribution, all inputs in tonnes/yr, returns Mt/yr.

    Eligible (clean = DAC + biogenic) origins get first claim on physical
    sequestration; fossil takes only the residual and is never credited.  When
    clean capture itself exceeds storage, the available storage is split within
    the clean bucket in proportion to captured volume (both are eligible).
    """
    clean_capture = cap_dac + cap_biogenic
    clean_stored = min(clean_capture, physical_seq_t)  # clean fills storage first
    fossil_residual_stored = max(0.0, min(cap_fossil, physical_seq_t - clean_stored))

    if clean_capture > 0.0:
        credited_dac_mt = (clean_stored * cap_dac / clean_capture) / 1e6
        credited_biogenic_mt = (clean_stored * cap_biogenic / clean_capture) / 1e6
    else:
        credited_dac_mt = 0.0
        credited_biogenic_mt = 0.0

    # Respect an explicit (binding) credit limit if one is configured.
    credited_total = credited_dac_mt + credited_biogenic_mt
    if not np.isnan(credit_limit_mt) and credited_total > credit_limit_mt and credited_total > 0:
        scale = credit_limit_mt / credited_total
        credited_dac_mt *= scale
        credited_biogenic_mt *= scale

    return credited_dac_mt, credited_biogenic_mt, fossil_residual_stored / 1e6


def _prorata_credited_total(cap_dac, cap_biogenic, cap_fossil, physical_seq_t, credit_limit_mt):
    """The pro-rata credited total (Mt/yr), for the side-by-side column only."""
    cap_total = cap_dac + cap_biogenic + cap_fossil
    if cap_total <= 0 or physical_seq_t <= 0:
        return 0.0
    total = (physical_seq_t * (cap_dac + cap_biogenic) / cap_total) / 1e6
    if not np.isnan(credit_limit_mt) and total > credit_limit_mt:
        total = credit_limit_mt
    return total


def build_waterfall_row(n: pypsa.Network, period: int) -> dict:
    prices = _get_cdr_credit_prices_for_period(
        cdr_credit_price=n.params.get("cdr_credit_price", 0.0),
        cdr_credit_scope=n.params.get("cdr_credit_scope") or [],
        cdr_credit_prices_by_scope=n.params.get("cdr_credit_prices_by_scope", {}),
        planning_horizons=str(period),
    )
    eligible_scopes = sorted(_eligible_cdr_scopes(n.params.get("cdr_credit_scope") or []))
    limit_dict = n.params.get("cdr_credit_limit_by_year") or n.params.get("cdr_credit_limit")
    credit_limit_mt = float(get(limit_dict, period)) if (limit_dict and period is not None) else np.nan

    capture = _annual_capture_proxy_by_origin(n)
    cap_dac = capture.get("dac", 0.0)
    cap_biogenic = capture.get("biogenic", 0.0)
    cap_fossil = capture.get("fossil", 0.0)
    physical_seq_t = _annual_physical_sequestration(n)

    credited_dac_mt, credited_biogenic_mt, fossil_residual_mt = _waterfall_credited(
        cap_dac, cap_biogenic, cap_fossil, physical_seq_t, credit_limit_mt
    )
    credited_total = credited_dac_mt + credited_biogenic_mt
    clean_capture_mt = (cap_dac + cap_biogenic) / 1e6
    physical_seq_mt = physical_seq_t / 1e6
    prorata_total = _prorata_credited_total(
        cap_dac, cap_biogenic, cap_fossil, physical_seq_t, credit_limit_mt
    )

    row = {
        "planning_horizon": period,
        "cdr_credit_timing": n.params.get("cdr_credit_timing", "capture"),
        "cdr_credit_standalone": bool(n.params.get("cdr_credit_standalone", False)),
        "method": "sequestration_waterfall",
        "accounting_error": "",
        "eligible_scopes": ",".join(eligible_scopes),
        "price_dac_eur_per_tco2": float(prices.get("dac", 0.0)),
        "price_biogenic_eur_per_tco2": float(prices.get("biogenic", 0.0)),
        "credit_limit_mtco2_per_yr": credit_limit_mt,
        "credited_dac_mtco2_per_yr": credited_dac_mt,
        "credited_biogenic_mtco2_per_yr": credited_biogenic_mt,
        "credited_total_mtco2_per_yr": credited_total,
        # Diagnostics specific to the waterfall convention.
        "fossil_residual_stored_mtco2_per_yr": fossil_residual_mt,
        "credited_prorata_total_mtco2_per_yr": prorata_total,
        # Same capture / sequestration context columns as the pro-rata file.
        "capture_proxy_dac_mtco2_per_yr": cap_dac / 1e6,
        "capture_proxy_biogenic_mtco2_per_yr": cap_biogenic / 1e6,
        "capture_proxy_total_mtco2_per_yr": clean_capture_mt,
        "capture_proxy_fossil_mtco2_per_yr": cap_fossil / 1e6,
        "physical_sequestration_mtco2_per_yr": physical_seq_mt,
    }
    # Invariants: credited must not exceed clean capture or physical sequestration.
    row["valid_credited_within_limit"] = (
        True if np.isnan(credit_limit_mt) else credited_total <= credit_limit_mt + TOL_MT
    )
    row["valid_credited_within_capture_proxy"] = credited_total <= clean_capture_mt + TOL_MT
    row["valid_credited_within_physical_sequestration"] = (
        credited_total <= physical_seq_mt + TOL_MT
    )
    return row


def process_one(net_path: Path, out_csv: Path, prorata_csv: Path, cfg_path: Path,
                year: int, dry_run: bool) -> dict:
    n = pypsa.Network(str(net_path))
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)
    n.params = _hydrate_params(cfg.get("sector", {}))

    row = build_waterfall_row(n, year)
    if not dry_run:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([row]).to_csv(out_csv, index=False)

    # Report the pro-rata value straight from its existing CSV when present.
    prorata_from_file = np.nan
    if prorata_csv.exists():
        pdf = pd.read_csv(prorata_csv)
        if not pdf.empty:
            prorata_from_file = float(pdf.iloc[0].get("credited_total_mtco2_per_yr", np.nan))

    return {
        "scenario": net_path.parts[-3],
        "year": year,
        "prorata_credited": prorata_from_file,
        "waterfall_credited": row["credited_total_mtco2_per_yr"],
        "physical_seq": row["physical_sequestration_mtco2_per_yr"],
        "clean_capture": row["capture_proxy_total_mtco2_per_yr"],
        "within_seq": row["valid_credited_within_physical_sequestration"],
        "within_capture": row["valid_credited_within_capture_proxy"],
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
        tag = stem.replace("base_", "")
        indiv = scen_dir / "csvs" / "individual"
        prorata_csv = indiv / f"cdr_credit_accounting_{tag}.csv"
        out_csv = indiv / f"cdr_credit_accounting_waterfall_{tag}.csv"
        cfg_path = scen_dir / "configs" / f"config.{stem}.yaml"
        if not cfg_path.exists():
            print(f"skip (no config): {cfg_path}")
            continue
        try:
            rec = process_one(net_path, out_csv, prorata_csv, cfg_path, year, args.dry_run)
            records.append(rec)
            flag = "" if rec["within_seq"] and rec["within_capture"] else "  <-- INVARIANT VIOLATED!"
            print(f"{rec['scenario']:<40} {year}  "
                  f"prorata {rec['prorata_credited']:7.1f} -> waterfall {rec['waterfall_credited']:7.1f} Mt  "
                  f"(seq {rec['physical_seq']:6.1f}, clean cap {rec['clean_capture']:6.1f}){flag}")
        except Exception as exc:  # noqa: BLE001
            print(f"FAILED {net_path}: {exc}")

    if records:
        df = pd.DataFrame(records)
        bad = df[~(df["within_seq"] & df["within_capture"])]
        print(f"\n{len(df)} networks processed"
              + (" (dry-run, nothing written)" if args.dry_run else ""))
        print(f"invariant violations: {len(bad)}")
        if len(bad):
            print(bad.to_string(index=False))


if __name__ == "__main__":
    main()
