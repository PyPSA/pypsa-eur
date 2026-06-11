#!/bin/bash
# =============================================================================
# gen_supply_curve_v13_ets.sh
# Generate the EU ETS sensitivity config set (v13_ets) from the v12 MEDIUM configs.
#
# v13_ets varies ONLY the exogenous EU ETS price path (costs.emission_prices.co2).
# Everything else is identical to v12 medium (cost CSVs, gas-for-industry CC cap,
# credit-price sweep S00-S10, solver options). Central ETS (119/279/463) is already
# solved as v12 medium and is NOT regenerated here.
#
# ETS paths (low / central / high), 2030 / 2040 / 2050, EUR/tCO2:
#   low     :  70 / 130 / 400   (Enerdata POLES core + removals-sensitivity lower bound)
#   central : 119 / 279 / 463   (Gunther et al. 2025, PRIMES/PIK)  <- = v12 medium, reused
#   high    : 160 / 400 / 630   (LSEG 90%-target + Enerdata removals-sensitivity upper bound)
#
# Run tags: medium-cost tag t10 -> t10L (low ETS) / t10H (high ETS).
# Results land in results/supply_curve_v13_ets/<run.name>/ via run.prefix.
#
# Usage (from repo root):
#   bash jobs/gen_supply_curve_v13_ets.sh
# =============================================================================

set -euo pipefail
shopt -s nullglob

REPO_ROOT="/work3/s240459/pypsa-eur-thesis"
SRC_DIR="$REPO_ROOT/config/Myruns/supply_curve_v12/supply_curve_v12_medium"
DST_ROOT="$REPO_ROOT/config/Myruns/supply_curve_v13_ets"

# ETS triples (2030 2040 2050)
LOW=(70 130 400)
HIGH=(160 400 630)

emit_variant() {
    local variant="$1" tag="$2" y30="$3" y40="$4" y50="$5"
    local dst="$DST_ROOT/ets_${variant}"
    rm -rf "$dst"
    mkdir -p "$dst"
    local n=0
    for src in "$SRC_DIR"/config.S*.yaml; do
        [[ -f "$src" && "$src" != *-fixed.yaml ]] || continue
        local base out
        base="$(basename "$src")"
        out="$dst/$base"
        sed \
            -e "s/t10hoB/${tag}hoB/g" \
            -e 's/prefix: "supply_curve_v12"/prefix: "supply_curve_v13_ets"/' \
            -e "s/^\(      2030:\) 119$/\1 ${y30}/" \
            -e "s/^\(      2040:\) 279$/\1 ${y40}/" \
            -e "s/^\(      2050:\) 463$/\1 ${y50}/" \
            "$src" > "$out"
        n=$((n+1))
    done
    echo "  ets_${variant}: wrote $n config(s) -> $dst  (ETS ${y30}/${y40}/${y50}, tag ${tag})"
}

echo "Generating v13_ets configs from: $SRC_DIR"
emit_variant low  t10L "${LOW[@]}"
emit_variant high t10H "${HIGH[@]}"

echo ""
echo "Verification (emission_prices.co2 triples found):"
for d in "$DST_ROOT"/ets_low "$DST_ROOT"/ets_high; do
    echo "  $d:"
    for f in "$d"/config.S*.yaml; do
        awk '/emission_prices/{f=1} f&&/2030:/{a=$2} f&&/2040:/{b=$2} f&&/2050:/{print a","b","$2; f=0}' "$f"
    done | sort | uniq -c | sed 's/^/    /'
done
echo ""
echo "Run-name tags:"
grep -h -m1 'name:' "$DST_ROOT"/ets_low/config.S01-*.yaml "$DST_ROOT"/ets_high/config.S01-*.yaml | sed 's/^/    /'
echo ""
echo "Done. Submit with: bash jobs/run_supply_curve_v13_ets.sh"
