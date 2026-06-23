"""
Check script for the biochar pipeline in pypsa-eur.
Run from the pypsa-eur root directory:

    python scripts/check_biochar_pipeline.py biochar_2050

Wildcards are read from the saved run config automatically.
Override any value on the CLI:
    --run-name NAME  --clusters N  --horizon YEAR  --sector-opts OPTS
    --config PATH    (direct path to a saved config YAML)
"""

import sys
from pathlib import Path

import pandas as pd

from _check_utils import parse_check_args, load_check_params

# ── Configuration (from saved run config) ────────────────────────────────────
_args = parse_check_args()
_p    = load_check_params(_args)

BASE_DIR         = _p["BASE_DIR"]
RDIR             = _p["RDIR"]
CLUSTERS         = _p["CLUSTERS"]
OPTS             = _p["OPTS"]
SECTOR_OPTS      = _p["SECTOR_OPTS"]
PLANNING_HORIZON = _p["PLANNING_HORIZON"]
WC               = _p["WC"]
RES              = _p["RES"]
RES_RUN          = _p["RES_RUN"]
RESULTS          = _p["RESULTS"]

_assign_capacity_duals = _p["cfg"].get("solving", {}).get("options", {}).get("assign_capacity_duals", False)

# ── Helpers ───────────────────────────────────────────────────────────────────
OK   = "  [OK]"
FAIL = "  [MISSING]"
WARN = "  [WARN]"

passed = []
failed = []


def check_file(path: Path, label: str) -> bool:
    exists = path.exists()
    status = OK if exists else FAIL
    size   = f"  ({path.stat().st_size / 1e6:.1f} MB)" if exists else ""
    print(f"{status}  {label}{size}")
    print(f"        {path}")
    if exists:
        passed.append(label)
    else:
        failed.append(label)
    return exists


def section(title: str):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print("=" * 60)


# ── 1. BUILD: biochar available land (determine_carbon_dioxide_removal_availability_matrix + build_available_land) ──
section("1. BUILD: biochar available land")

biochar_csv = RES_RUN / f"biochar_available_land_s_{CLUSTERS}.csv"

csv_ok = check_file(biochar_csv, f"Biochar available land CSV  s_{CLUSTERS}")

if csv_ok:
    df = pd.read_csv(biochar_csv)
    print(f"        Rows: {len(df)}  |  Columns: {list(df.columns)}")
    if "potential [sqkm]" in df.columns:
        total = df["potential [sqkm]"].sum()
        print(f"        Total potential area: {total:,.0f} km²")
        print(f"        Per-node [km²] — min: {df['potential [sqkm]'].min():.1f}, "
              f"mean: {df['potential [sqkm]'].mean():.1f}, "
              f"max: {df['potential [sqkm]'].max():.1f}")
    if df.isnull().any().any():
        print(f"{WARN}  NaN values detected in biochar potentials CSV.")

# ── 2. Pre-network (prepare_sector_network) ───────────────────────────────────
section("2. PRE-NETWORK: sector-coupled (prepare_sector_network)")

prenet_path = RES_RUN / "networks" / f"{WC}.nc"
prenet_ok   = check_file(prenet_path, f"Pre-network  {WC}.nc")

if prenet_ok:
    try:
        import pypsa

        n = pypsa.Network(str(prenet_path))

        # --- Carriers ---
        biochar_carriers = [c for c in n.carriers.index if "biochar" in c.lower()]
        print(f"\n  Carriers with 'biochar': {biochar_carriers}")

        if "co2 biochar" in n.carriers.index:
            co2_em = n.carriers.at["co2 biochar", "co2_emissions"]
            # co2_emissions = 0.0 is correct: sequestration is modelled via the
            # biochar link drawing from co2 atmosphere (bus0), not via carrier attribute.
            print(f"    co2 biochar co2_emissions = {co2_em}  (expected 0.0 — sequestration via network flow)")
            if abs(co2_em) > 1e-6:
                print(f"{WARN}  co2_emissions != 0.0 — check add_biochar()!")

        # --- Buses ---
        bc_buses = n.buses[n.buses.carrier.str.contains("biochar", case=False, na=False)]
        print(f"\n  Buses (carrier contains 'biochar'): {len(bc_buses)}")
        if not bc_buses.empty:
            print(f"    Carriers present: {bc_buses['carrier'].unique().tolist()}")
            print(f"    Sample (first 5):\n{bc_buses[['location','carrier','unit']].head().to_string()}")

        # --- Links (co2 biochar only) ---
        co2_bc_links = n.links[n.links.carrier == "co2 biochar"]
        print(f"\n  Links (carrier == 'co2 biochar'): {len(co2_bc_links)}")
        if not co2_bc_links.empty:
            cols = [c for c in ["bus0", "bus1", "carrier", "p_nom_extendable", "capital_cost"] if c in co2_bc_links.columns]
            print(co2_bc_links[cols].head(10).to_string())

        # --- Stores (co2 biochar only) ---
        co2_bc_stores = n.stores[n.stores.carrier == "co2 biochar"]
        n_nodes = len(n.buses[n.buses.carrier == "AC"])
        print(f"\n  Stores (carrier == 'co2 biochar'): {len(co2_bc_stores)}  (nodes: {n_nodes})")
        # Note: stores with carrier 'biochar heat' (heat waste sinks) are separate
        # and counted independently below.
        if not co2_bc_stores.empty:
            cols = [c for c in ["bus", "carrier", "e_nom_max", "capital_cost", "e_nom_extendable"] if c in co2_bc_stores.columns]
            print(co2_bc_stores[cols].head(10).to_string())
            if "e_nom_max" in co2_bc_stores.columns:
                finite = co2_bc_stores["e_nom_max"][co2_bc_stores["e_nom_max"] < 1e18]
                print(f"\n    e_nom_max (finite) [tCO2] — "
                      f"min: {finite.min():.0f}, mean: {finite.mean():.0f}, max: {finite.max():.0f}")
                total_max = finite.sum()
                print(f"    Total e_nom_max (max potential): {total_max:,.0f} tCO2  "
                      f"({total_max / 1e6:.3f} MtCO2)")

        # --- Stores (biochar heat waste sinks) ---
        heat_bc_stores = n.stores[n.stores.carrier == "biochar heat"]
        if not heat_bc_stores.empty:
            print(f"\n  Stores (carrier == 'biochar heat', heat waste sinks): {len(heat_bc_stores)}")

        if not biochar_carriers:
            print(f"\n{WARN}  No biochar carriers found — add_biochar may NOT have run!")
        elif co2_bc_links.empty or co2_bc_stores.empty:
            print(f"\n{WARN}  Missing co2 biochar links or stores — check add_biochar() execution!")
        else:
            print(f"\n{OK}  Biochar components present in pre-network.")

    except Exception as e:
        print(f"\n{WARN}  Could not load pre-network: {e}")
else:
    print(f"\n{FAIL}  Pre-network missing — prepare_sector_network has not run yet.")

# ── 3. Solved network (solve_sector_network) ──────────────────────────────────
section("3. OPTIMAL SOLUTION: solved network (solve_sector_network)")

opt_path = RESULTS / "networks" / f"{WC}.nc"
opt_ok   = check_file(opt_path, f"Optimal network  {WC}.nc")

if opt_ok:
    try:
        import pypsa

        _dual_status = "enabled" if _assign_capacity_duals else "disabled (set assign_capacity_duals: true to enable)"
        print(f"\n  assign_capacity_duals: {_dual_status}")

        n_opt = pypsa.Network(str(opt_path))

        co2_bc_links  = n_opt.links[n_opt.links.carrier == "co2 biochar"]
        co2_bc_stores = n_opt.stores[n_opt.stores.carrier == "co2 biochar"]

        # --- Links optimal (co2 biochar) ---
        print(f"\n  Links (carrier == 'co2 biochar'): {len(co2_bc_links)}")
        if not co2_bc_links.empty and "p_nom_opt" in co2_bc_links.columns:
            active = co2_bc_links[co2_bc_links["p_nom_opt"] > 0]
            total  = co2_bc_links["p_nom_opt"].sum()
            print(f"  Links with p_nom_opt > 0: {len(active)}")
            print(f"  Total p_nom_opt (co2 biochar links): {total:,.2f} MW")
            if active.empty:
                print(f"{WARN}  All co2 biochar links have p_nom_opt = 0 (not deployed).")
            else:
                print(f"{OK}  co2 biochar links deployed in optimal solution.")
                print(active[["bus0", "bus1", "carrier", "p_nom_opt"]].to_string())

        # --- Stores optimal (co2 biochar only) ---
        print(f"\n  Stores (carrier == 'co2 biochar'): {len(co2_bc_stores)}")
        if not co2_bc_stores.empty and "e_nom_opt" in co2_bc_stores.columns:
            active = co2_bc_stores[co2_bc_stores["e_nom_opt"] > 0]
            total  = co2_bc_stores["e_nom_opt"].sum()
            print(f"  Stores with e_nom_opt > 0: {len(active)}")
            print(f"  Total e_nom_opt: {total:,.0f} tCO2  ({total/1e6:.3f} MtCO2)")
            if active.empty:
                print(f"{WARN}  All co2 biochar stores have e_nom_opt = 0 (not deployed).")
            else:
                print(f"{OK}  Biochar stores deployed. Total CO2 stored: {total/1e6:.3f} MtCO2")
                print(f"\n  Per-node e_nom_opt [tCO2] stats:")
                print(f"    min:  {co2_bc_stores['e_nom_opt'].min():,.0f}")
                print(f"    mean: {co2_bc_stores['e_nom_opt'].mean():,.0f}")
                print(f"    max:  {co2_bc_stores['e_nom_opt'].max():,.0f}")

        # --- Shadow price of e_nom_max (assign_capacity_duals) ---
        if _assign_capacity_duals:
            if "mu_ext_e_nom_upper" in co2_bc_stores.columns:
                mu = co2_bc_stores["mu_ext_e_nom_upper"].dropna()
                binding = mu[mu.abs() > 1e-6]
                print(f"\n  mu_ext_e_nom_upper (shadow price of e_nom_max) [{len(mu)} stores]:")
                print(f"    non-zero (binding): {len(binding)}")
                if not mu.empty:
                    print(f"    min:  {mu.min():,.4f}")
                    print(f"    mean: {mu.mean():,.4f}")
                    print(f"    max:  {mu.max():,.4f}")
                if not binding.empty:
                    print(f"\n    Binding nodes (top 10):")
                    print(binding.sort_values().tail(10).to_string())
            else:
                print(f"\n{WARN}  mu_ext_e_nom_upper not found in stores — was assign_capacity_duals: true when solving?")

        if co2_bc_links.empty and co2_bc_stores.empty:
            print(f"\n{WARN}  No co2 biochar components in optimal network!")

    except Exception as e:
        print(f"\n{WARN}  Could not load optimal network: {e}")
else:
    print(f"\n{FAIL}  Optimal network missing — solve_sector_network has not run yet.")

# ── Summary ───────────────────────────────────────────────────────────────────
section("SUMMARY")
print(f"  Run:            {RDIR if RDIR else '(none)'}")
print(f"  Wildcard:       {WC}")
print(f"  Resources dir:  {RES_RUN.resolve()}")
print(f"  Results dir:    {RESULTS.resolve()}")
print()
print(f"  Passed: {len(passed)}")
print(f"  Failed: {len(failed)}")
if failed:
    print(f"\n  Missing files:")
    for f in failed:
        print(f"    - {f}")
