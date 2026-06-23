"""
Check script for the perennials pipeline in pypsa-eur.
Run from the pypsa-eur root directory:

    python scripts/check_perennials_pipeline.py peren_2050

Wildcards are read from the saved run config automatically.
Override any value on the CLI:
    --run-name NAME  --clusters N  --horizon YEAR  --sector-opts OPTS
    --shared-resources NAME
    --config PATH    (direct path to a saved config YAML)
"""

import sys
from pathlib import Path

import pandas as pd

from _check_utils import parse_check_args, load_check_params

# ── Configuration (from config file + CLI overrides) ──────────────────────────
_args = parse_check_args(extra_args=[
    (["--shared-resources"],
     {"default": None,
      "metavar": "NAME",
      "help": "Override shared resources directory name. "
              "Reads run.shared_resources.policy from config if not set, "
              "falls back to run name."}),
])
_p = load_check_params(_args)

BASE_DIR         = _p["BASE_DIR"]
RDIR             = _p["RDIR"]
CLUSTERS         = _p["CLUSTERS"]
OPTS             = _p["OPTS"]
SECTOR_OPTS      = _p["SECTOR_OPTS"]
PLANNING_HORIZON = _p["PLANNING_HORIZON"]
WC               = _p["WC"]
RES              = _p["RES"]
RESULTS          = _p["RESULTS"]

_shared_policy = (
    _args.shared_resources
    or _p["cfg"].get("run", {}).get("shared_resources", {}).get("policy")
    or RDIR
)
SHARED_RES_POLICY = str(_shared_policy)
SHARED_RES        = BASE_DIR / "resources" / SHARED_RES_POLICY

# ── Helpers ──────────────────────────────────────────────────────────────────
OK   = "  [OK]"
FAIL = "  [MISSING]"
WARN = "  [WARN]"

def check_file(path: Path, label: str) -> bool:
    exists = path.exists()
    status = OK if exists else FAIL
    size = f"  ({path.stat().st_size / 1e6:.1f} MB)" if exists else ""
    print(f"{status}  {label}{size}")
    print(f"        {path}")
    return exists


def section(title: str):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print('='*60)


# ── 1. Retrieve rule outputs ─────────────────────────────────────────────────
section("1. RETRIEVE: eurostat crops + perennial yields")

check_file(
    SHARED_RES / "eurostat_apro_cpshr_nuts2_raw.csv",
    "Eurostat NUTS2 crops"
)
check_file(
    SHARED_RES / "eurostat_apro_cpshr_nuts0_raw.csv",
    "Eurostat NUTS0 crops"
)
check_file(
    SHARED_RES / "perennials_yields_1G_biofuels.csv",
    "Perennials yields (all NUTS)"
)

# ── 2. Build rule outputs ────────────────────────────────────────────────────
section("2. BUILD: perennial potentials (clustered)")

yields_clustered = SHARED_RES / f"perennials_yields_1G_biofuels_s_{CLUSTERS}.csv"
if check_file(yields_clustered, f"Perennials yields clustered s_{CLUSTERS}"):
    df = pd.read_csv(yields_clustered)
    print(f"        Rows: {len(df)}  |  Columns: {list(df.columns)}")
    if 'perennials' in df.columns:
        print(f"        perennials [tDM/ha] — min: {df['perennials'].min():.2f}, mean: {df['perennials'].mean():.2f}, max: {df['perennials'].max():.2f}")
    else:
        print(f"        {WARN} no 'perennials' column found!")
    biofuels_1G_cols = [c for c in df.columns if 'biofuels_1G' in c]
    if biofuels_1G_cols:
        print(f"        biofuels_1G columns: {biofuels_1G_cols}")
        for col in biofuels_1G_cols:
            print(f"          {col}: mean={df[col].mean():.3f} MWh/ha")
    else:
        print(f"        {WARN} No 'biofuels_1G_*' columns found — biomass.classes in config.default.yaml")
        print(f"        {WARN} must use biofuels_1G_* names, NOT 'not included', for perennials to work!")

# ── 3. Pre-network (prepare_sector_network output) ───────────────────────────
section("3. PRE-NETWORK: sector-coupled (prepare_sector_network)")

prenet_path = SHARED_RES / "networks" / f"{WC}.nc"
prenet_ok = check_file(prenet_path, f"Pre-network {WC}.nc")

if prenet_ok:
    try:
        import pypsa
        n = pypsa.Network(str(prenet_path))

        # Check carriers
        perenn_carriers = [c for c in n.carriers.index if "perennial" in c.lower()]
        print(f"\n  Carriers with 'perennial': {perenn_carriers}")

        # Check links
        perenn_links = n.links[n.links.carrier.str.contains("perennial", case=False, na=False)]
        print(f"  Links (carrier=perennial): {len(perenn_links)}")
        if not perenn_links.empty:
            print(perenn_links[["bus0", "bus1", "carrier", "p_nom"]].to_string(index=True))

        # Check stores
        perenn_stores = n.stores[n.stores.carrier.str.contains("perennial", case=False, na=False)]
        print(f"  Stores (carrier=perennial store): {len(perenn_stores)}")
        if not perenn_stores.empty:
            print(perenn_stores[["bus", "carrier", "e_nom_max"]].head(10).to_string(index=True))
            if "e_nom_max" in perenn_stores.columns:
                finite_max = perenn_stores["e_nom_max"][perenn_stores["e_nom_max"] < 1e18]
                total_cap = finite_max.sum()
                print(f"\n  Total store capacity (e_nom_max): {total_cap:,.0f} tCO2  ({total_cap/1e6:.3f} MtCO2)")
                if total_cap == 0:
                    print(f"  {WARN} ALL stores have e_nom_max=0!")
                    print(f"  {WARN} This usually means biomass.classes in config.default.yaml")
                    print(f"  {WARN} is missing biofuels_1G_* entries — check and re-run prepare_sector_network.")

        if not perenn_carriers:
            print(f"\n{WARN}  No perennial carriers found — add_perennials may NOT have run!")
        else:
            print(f"\n{OK}  Perennial components found in pre-network.")

    except Exception as e:
        print(f"\n{WARN}  Could not load network: {e}")
else:
    print(f"\n{FAIL}  Pre-network missing — prepare_sector_network has not run yet.")

# ── 4. Optimal solution ──────────────────────────────────────────────────────
section("4. OPTIMAL SOLUTION: solved network (solve_sector_network)")

opt_path = RESULTS / "networks" / f"{WC}.nc"
opt_ok = check_file(opt_path, f"Optimal network {WC}.nc")

if opt_ok:
    try:
        import pypsa
        n_opt = pypsa.Network(str(opt_path))

        perenn_links = n_opt.links[n_opt.links.carrier.str.contains("perennial", case=False, na=False)]
        perenn_stores = n_opt.stores[n_opt.stores.carrier.str.contains("perennial", case=False, na=False)]

        print(f"\n  Links (carrier=perennial): {len(perenn_links)}")
        if not perenn_links.empty:
            cols = ["carrier", "p_nom_opt"] if "p_nom_opt" in perenn_links.columns else ["carrier", "p_nom"]
            print(perenn_links[cols].to_string(index=True))
            if "p_nom_opt" in perenn_links.columns:
                active = perenn_links[perenn_links["p_nom_opt"] > 0]
                print(f"\n  Links with p_nom_opt > 0: {len(active)}")
                if active.empty:
                    print(f"{WARN}  Perennial links exist but have zero optimal capacity.")
                else:
                    print(f"{OK}  Perennial links are deployed in the optimal solution.")

        print(f"\n  Stores (carrier=perennial store): {len(perenn_stores)}")
        if not perenn_stores.empty:
            cols = ["carrier", "e_nom_opt"] if "e_nom_opt" in perenn_stores.columns else ["carrier", "e_nom"]
            print(perenn_stores[cols].to_string(index=True))
            if "e_nom_opt" in perenn_stores.columns:
                active = perenn_stores[perenn_stores["e_nom_opt"] > 0]
                total = perenn_stores["e_nom_opt"].sum()
                print(f"\n  Stores with e_nom_opt > 0: {len(active)}")
                print(f"  Total e_nom_opt: {total:,.0f} tCO2  ({total / 1e6:.3f} MtCO2)")
                if active.empty:
                    print(f"{WARN}  All perennial stores have e_nom_opt = 0 (not deployed).")
                else:
                    print(f"\n  Per-node e_nom_opt [tCO2] stats:")
                    print(f"    min:  {perenn_stores['e_nom_opt'].min():,.0f}")
                    print(f"    mean: {perenn_stores['e_nom_opt'].mean():,.0f}")
                    print(f"    max:  {perenn_stores['e_nom_opt'].max():,.0f}")

        if perenn_links.empty and perenn_stores.empty:
            print(f"\n{WARN}  No perennial components in optimal network!")

    except Exception as e:
        print(f"\n{WARN}  Could not load optimal network: {e}")
else:
    print(f"\n{FAIL}  Optimal network missing — solve_sector_network has not run yet.")

# ── Summary ──────────────────────────────────────────────────────────────────
section("SUMMARY")
print(f"  Run:             {RDIR}")
print(f"  Wildcard:        {WC}")
print(f"  Shared resources:{SHARED_RES.resolve()}")
print(f"  Results dir:     {RESULTS.resolve()}")
