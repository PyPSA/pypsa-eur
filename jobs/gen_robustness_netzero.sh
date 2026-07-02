#!/bin/bash
# =============================================================================
# gen_robustness_netzero.sh
# Generate the NET-ZERO sensitivity config set (Robustness_Analysis_NetZero) from
# the Scenario 1 MEDIUM-cost configs.
# Internal tag: v14_netzero.
#
# The main analysis runs with NO economy-wide CO2 cap (co2_budget: null), so CDR
# deploys purely on the standalone credit economics. This sensitivity turns the cap
# back ON: an economy-wide CO2 budget declining to net-zero in 2050, identical in
# every other respect to v12 medium (cost CSVs custom_costs_medium_hoB.csv, central
# ETS 119/279/463, gas-for-industry CC cap 1.0, credit sweep S00-S10, solver opts).
#
# CO2 budget trajectory (fraction of 1990 emissions, per planning horizon):
#   2030: 0.45 , 2040: 0.10 , 2050: 0.0  (PyPSA-Eur default, net-zero at 2050)
# Implemented through add_co2limit() as a global co2_atmosphere GlobalConstraint
# per planning year (scripts/prepare_sector_network.py).
#
# Only three things differ from the v12 medium configs:
#   1. run.name tag : t10hoB -> t10NZhoB
#   2. run.prefix   : supply_curve_v12 -> supply_curve_v14_netzero
#   3. co2_budget   : null -> {2030: 0.45, 2040: 0.10, 2050: 0.0}
#
# Usage (from repo root):
#   bash jobs/gen_robustness_netzero.sh
# =============================================================================

set -euo pipefail
shopt -s nullglob

REPO_ROOT="/work3/s240459/pypsa-eur-thesis"
SRC_DIR="$REPO_ROOT/config/Thesis_Runs/Scenario_1_Deployment_Response/Sensitivity_Medium"
DST_DIR="$REPO_ROOT/config/Thesis_Runs/Robustness_Analysis_NetZero/netzero_default"

echo "Generating v14_netzero configs from: $SRC_DIR"
rm -rf "$DST_DIR"
mkdir -p "$DST_DIR"

n=0
for src in "$SRC_DIR"/config.S*.yaml; do
    [[ -f "$src" && "$src" != *-fixed.yaml ]] || continue
    base="$(basename "$src")"
    out="$DST_DIR/$base"
    # 1+2: retag run name and prefix; 3: replace the `co2_budget: null` line with the
    #      net-zero trajectory block (awk handles the single-line -> multi-line swap).
    sed \
        -e "s/t10hoB/t10NZhoB/g" \
        -e 's/prefix: "supply_curve_v12"/prefix: "supply_curve_v14_netzero"/' \
        "$src" \
    | awk '
        /^co2_budget:/ {
            print "co2_budget:  # Net-zero sensitivity: economy-wide CO2 cap (fraction of 1990) per horizon"
            print "  2030: 0.45"
            print "  2040: 0.10"
            print "  2050: 0.0"
            next
        }
        { print }
    ' > "$out"
    n=$((n+1))
done
echo "  netzero_default: wrote $n config(s) -> $DST_DIR"

echo ""
echo "Verification:"
echo "  run.name / run.prefix (S02):"
grep -hE 'name:|prefix:' "$DST_DIR"/config.S02-*.yaml | head -2 | sed 's/^/    /'
echo "  co2_budget block (S02):"
awk '/^co2_budget:/{f=1} f{print "    "$0} f&&/2050:/{exit}' "$DST_DIR"/config.S02-*.yaml
echo ""
echo "Done. Submit with: bash jobs/run_robustness_netzero.sh"
