#!/bin/bash
# =============================================================================
# run_supply_curve_v13_ets.sh
# Submit one LSF job per scenario for the EU ETS price SENSITIVITY (v13_ets).
# Varies only the exogenous ETS path (costs.emission_prices.co2) at MEDIUM cost.
# Central ETS (119/279/463) is already solved as v12 medium and is NOT re-run here.
#
# ETS paths (2030/2040/2050, EUR/tCO2):  low = 70/130/400 ,  high = 160/400/630
# Results go to: results/supply_curve_v13_ets/<run.name>/
#
# Usage (from repo root):
#   bash jobs/run_supply_curve_v13_ets.sh                 # submit ALL (low + high, S00-S10)
#   bash jobs/run_supply_curve_v13_ets.sh low             # only the low-ETS variant
#   bash jobs/run_supply_curve_v13_ets.sh high S01 S03    # high-ETS, specific credit prices
#
# Flags:
#   --clean   Delete existing results for selected scenarios before submitting
#
# NOTE on HPC throughput: per ticket #41882 settle on <= 6-7 concurrent jobs.
# 22 jobs total here -> submit in waves, e.g. `... low` first, then `... high`.
#
# Monitor:
#   bjobs -w
#   bjobs -J 'v13e-*'
# =============================================================================

set -euo pipefail
shopt -s nullglob

REPO_ROOT="/work3/s240459/pypsa-eur-thesis"
BASE_DIR="$REPO_ROOT/config/Myruns/supply_curve_v13_ets"
LOG_DIR="$REPO_ROOT/logs/supply_curve_v13_ets"
mkdir -p "$LOG_DIR"

WALLTIME="72:00"
CORES=8
MEM_MB=15000
QUEUE="hpc"

# --- parse args: optional variant filter (low/high) + optional scenario prefixes + --clean
CLEAN=0
VARIANTS=()
PREFIXES=()
for arg in "$@"; do
    case "$arg" in
        --clean) CLEAN=1 ;;
        low|high) VARIANTS+=("$arg") ;;
        *) PREFIXES+=("$arg") ;;
    esac
done
[[ ${#VARIANTS[@]} -eq 0 ]] && VARIANTS=(low high)

# --- collect configs
CONFIGS=()
for v in "${VARIANTS[@]}"; do
    DIR="$BASE_DIR/ets_${v}"
    [[ -d "$DIR" ]] || { echo "Warning: missing $DIR — run jobs/gen_supply_curve_v13_ets.sh first"; continue; }
    if [[ ${#PREFIXES[@]} -gt 0 ]]; then
        for prefix in "${PREFIXES[@]}"; do
            for candidate in "$DIR"/config.${prefix}.yaml "$DIR"/config.${prefix}-*.yaml; do
                [[ -f "$candidate" && "$candidate" != *-fixed.yaml ]] && CONFIGS+=("$candidate")
            done
        done
    else
        for candidate in "$DIR"/config.S*.yaml; do
            [[ -f "$candidate" && "$candidate" != *-fixed.yaml ]] && CONFIGS+=("$candidate")
        done
    fi
done

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
    echo "Error: no configs selected."
    exit 1
fi

# --- optional clean
if [[ "$CLEAN" -eq 1 ]]; then
    RESULTS_TO_CLEAN=()
    for cfg in "${CONFIGS[@]}"; do
        run_name=$(grep -m1 'name:' "$cfg" | awk '{print $2}' | tr -d '"')
        run_prefix=$(grep -m1 'prefix:' "$cfg" | awk '{print $2}' | tr -d '"')
        d="$REPO_ROOT/results/$run_prefix/$run_name"
        [[ -d "$d" ]] && RESULTS_TO_CLEAN+=("$d")
    done
    if [[ ${#RESULTS_TO_CLEAN[@]} -eq 0 ]]; then
        echo "No existing result directories found — nothing to clean."
    else
        echo "The following result directories will be deleted:"
        for d in "${RESULTS_TO_CLEAN[@]}"; do echo "  $d  ($(du -sh "$d" 2>/dev/null | cut -f1))"; done
        read -r -p "Confirm deletion? [y/N] " confirm
        if [[ "$confirm" =~ ^[Yy]$ ]]; then
            for d in "${RESULTS_TO_CLEAN[@]}"; do echo "  Deleting $d ..."; rm -rf "$d"; done
        else
            echo "Aborted."; exit 0
        fi
    fi
fi

echo ""
echo "Submitting ${#CONFIGS[@]} ETS-sensitivity scenario(s) [v13_ets, medium cost, ETS low=70/130/400 high=160/400/630]..."
echo ""

for cfg in "${CONFIGS[@]}"; do
    basename_cfg=$(basename "$cfg" .yaml)
    scenario_id="${basename_cfg#config.}"
    variant_dir=$(basename "$(dirname "$cfg")")   # ets_low / ets_high
    variant="${variant_dir#ets_}"
    job_name="v13e-${variant}-${scenario_id}"

    bsub <<BSUB_SCRIPT
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -q ${QUEUE}
#BSUB -n ${CORES}
#BSUB -W ${WALLTIME}
#BSUB -R "rusage[mem=${MEM_MB}] span[hosts=1]"
#BSUB -o ${LOG_DIR}/${variant}_${scenario_id}_%J.out
#BSUB -e ${LOG_DIR}/${variant}_${scenario_id}_%J.err

set -euo pipefail

module purge
module load python3/3.10.18
module load gurobi/12.0.3

cd "${REPO_ROOT}"

export PATH="\$HOME/.pixi/bin:\$HOME/.local/bin:\$PATH"
PIXI_BIN="\$(command -v pixi || true)"
[[ -z "\$PIXI_BIN" && -x "\$HOME/.pixi/bin/pixi" ]] && PIXI_BIN="\$HOME/.pixi/bin/pixi"
[[ -z "\$PIXI_BIN" ]] && { echo "Error: pixi not found"; exit 127; }

WORK_ROOT="/work3/\$USER"
export TMPDIR="\${__LSF_JOB_TMPDIR__:-\$WORK_ROOT/tmp}"
export XDG_CACHE_HOME="\${XDG_CACHE_HOME:-\$WORK_ROOT/.cache}"
export SNAKEMAKE_OUTPUT_CACHE="\${SNAKEMAKE_OUTPUT_CACHE:-\$WORK_ROOT/.snakemake_cache}"
mkdir -p "\$TMPDIR" "\$XDG_CACHE_HOME" "\$SNAKEMAKE_OUTPUT_CACHE"

echo "Host:     \$(hostname)"
echo "Scenario: ${scenario_id} [v13_ets ETS-${variant}, medium cost]"
echo "Config:   ${cfg}"
echo "Started:  \$(date)"

run_snakemake() {
  "\$PIXI_BIN" run snakemake \
    -j ${CORES} \
    --nolock \
    --rerun-incomplete \
    --rerun-triggers params input code \
    --keep-going \
    --latency-wait 120 \
    --printshellcmds \
    --configfile "${cfg}"
}

set +e
run_snakemake 2>&1
status=\$?
set -e

if [[ "\$status" -ne 0 ]]; then
    echo "Snakemake failed with status \$status"
    exit "\$status"
fi

echo "Completed: \$(date)"
BSUB_SCRIPT

    echo "  Submitted: ${variant}/${scenario_id}  (job name: $job_name)"
done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  bjobs -w"
echo "  bjobs -J 'v13e-*'"
echo ""
echo "Results: results/supply_curve_v13_ets/"
echo "Logs:    $LOG_DIR"
