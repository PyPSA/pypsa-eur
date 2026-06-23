#!/bin/bash
# =============================================================================
# gen_robustness_ets.sh
# Generate the EU ETS sensitivity config set (Robustness_Analysis_ETS) from the
# Scenario 1 MEDIUM-cost configs.
# (formerly gen_supply_curve_v13_ets.sh; internal tag: v13_ets.)
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
# Additional 2050-endpoint points (keep 2030/2040, vary only 2050):
#   low430  :  70 / 130 / 430   (extra resolution near the low endpoint)
#   high525 : 160 / 400 / 525   (extra resolution below the high endpoint)
#
# 2050 infill sweep (medium cost) to pin the binding-credit-price flips exactly:
#   low : 405 / 420 / 435 / 450  -> resolves the 100->150 EUR/t flip on the low path
#   high: 500 / 550              -> resolves the flip on the high path
#
# Alternative high anchor (LOWER 2030/2040 than the 160/400 high family):
#   high2   : 150 / 360 / 500   (tag t10H2)  -> a separate high trajectory, not a
#                                               2050-infill of the 160/400 high path
#
# Run tags: medium-cost tag t10 -> t10L<y50> (low path) / t10H<y50> (high path),
#   e.g. t10L (=400) / t10H (=630) / t10L430 / t10H525 / t10L405 / t10H500 ...
# Results land in results/supply_curve_v13_ets/<run.name>/ via run.prefix.
#
# Usage (from repo root):
#   bash jobs/gen_robustness_ets.sh
# =============================================================================

set -euo pipefail
shopt -s nullglob

REPO_ROOT="/work3/s240459/pypsa-eur-thesis"
SRC_DIR="$REPO_ROOT/config/Thesis_Runs/Scenario_1_Deployment_Response/Sensitivity_Medium"
DST_ROOT="$REPO_ROOT/config/Thesis_Runs/Robustness_Analysis_ETS"

# Variant table: "variant_name  tag  2030 2040 2050"
# (low/high 2050-endpoints + the 2050 infill points; central = v12 medium, not regenerated)
VARIANT_TABLE=(
    "low      t10L     70  130 400"
    "high     t10H    160  400 630"
    "low430   t10L430  70  130 430"
    "high525  t10H525 160  400 525"
    "low405   t10L405  70  130 405"
    "low420   t10L420  70  130 420"
    "low435   t10L435  70  130 435"
    "low450   t10L450  70  130 450"
    "high500  t10H500 160  400 500"
    "high550  t10H550 160  400 550"
    "high2    t10H2   150  360 500"
)

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
for row in "${VARIANT_TABLE[@]}"; do
    read -r vname tag y30 y40 y50 <<< "$row"
    emit_variant "$vname" "$tag" "$y30" "$y40" "$y50"
done

echo ""
echo "Verification (emission_prices.co2 triples found):"
for row in "${VARIANT_TABLE[@]}"; do
    read -r vname _ _ _ _ <<< "$row"
    d="$DST_ROOT/ets_${vname}"
    echo "  $d:"
    for f in "$d"/config.S*.yaml; do
        awk '/emission_prices/{f=1} f&&/2030:/{a=$2} f&&/2040:/{b=$2} f&&/2050:/{print a","b","$2; f=0}' "$f"
    done | sort | uniq -c | sed 's/^/    /'
done
echo ""
echo "Run-name tags (S01 of each variant):"
for row in "${VARIANT_TABLE[@]}"; do
    read -r vname _ _ _ _ <<< "$row"
    grep -h -m1 'name:' "$DST_ROOT/ets_${vname}"/config.S01-*.yaml | sed 's/^/    /'
done
echo ""
echo "Done. Submit with: bash jobs/run_robustness_ets.sh"
