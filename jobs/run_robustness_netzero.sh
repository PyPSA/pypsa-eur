#!/bin/bash
# =============================================================================
# run_robustness_netzero.sh
# Submit one LSF job per scenario for the NET-ZERO sensitivity robustness analysis.
# Internal tag: v14_netzero.
#
# Turns the economy-wide CO2 cap back ON (co2_budget 0.45/0.10/0.0 for 2030/2040/2050,
# net-zero at 2050) at MEDIUM cost and central ETS. Everything else identical to v12
# medium. Comparator is the existing results/supply_curve_v12/S*t10hoB-... runs
# (same configs, co2_budget: null) — no re-run of the comparator needed.
#
# Results go to: results/supply_curve_v14_netzero/<run.name>/
#
# Usage (from repo root):
#   bash jobs/run_robustness_netzero.sh                 # submit S00-S10 (11 jobs)
#   bash jobs/run_robustness_netzero.sh S00 S02 S03 S04 # only specific credit prices
#
# Flags:
#   --clean   Delete existing results for selected scenarios before submitting
#
# NOTE on HPC throughput: per ticket #41882 settle on <= 6-7 concurrent jobs.
# 11 jobs total -> submit in waves if the queue is busy.
#
# Monitor:
#   bjobs -w
#   bjobs -J 'v14nz-*'
# =============================================================================

set -euo pipefail
shopt -s nullglob

REPO_ROOT="/work3/s240459/pypsa-eur-thesis"
CFG_DIR="$REPO_ROOT/config/Thesis_Runs/Robustness_Analysis_NetZero/netzero_default"
LOG_DIR="$REPO_ROOT/logs/supply_curve_v14_netzero"
mkdir -p "$LOG_DIR"

WALLTIME="36:00"
CORES=8
MEM_MB=4000
QUEUE="hpc"

# --- parse args: optional scenario prefixes (e.g. S01) + --clean
CLEAN=0
PREFIXES=()
for arg in "$@"; do
    if [[ "$arg" == "--clean" ]]; then CLEAN=1; continue; fi
    PREFIXES+=("$arg")
done

[[ -d "$CFG_DIR" ]] || { echo "Missing $CFG_DIR — run jobs/gen_robustness_netzero.sh first"; exit 1; }

# --- collect configs
CONFIGS=()
if [[ ${#PREFIXES[@]} -gt 0 ]]; then
    for prefix in "${PREFIXES[@]}"; do
        for candidate in "$CFG_DIR"/config.${prefix}.yaml "$CFG_DIR"/config.${prefix}-*.yaml; do
            [[ -f "$candidate" && "$candidate" != *-fixed.yaml ]] && CONFIGS+=("$candidate")
        done
    done
else
    for candidate in "$CFG_DIR"/config.S*.yaml; do
        [[ -f "$candidate" && "$candidate" != *-fixed.yaml ]] && CONFIGS+=("$candidate")
    done
fi

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
echo "Submitting ${#CONFIGS[@]} net-zero-sensitivity scenario(s) [v14_netzero, medium cost, central ETS]..."
echo ""

for cfg in "${CONFIGS[@]}"; do
    basename_cfg=$(basename "$cfg" .yaml)
    scenario_id="${basename_cfg#config.}"
    job_name="v14nz-${scenario_id}"

    bsub <<BSUB_SCRIPT
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -q ${QUEUE}
#BSUB -n ${CORES}
#BSUB -W ${WALLTIME}
#BSUB -R "rusage[mem=${MEM_MB}] span[hosts=1]"
#BSUB -o ${LOG_DIR}/${scenario_id}_%J.out
#BSUB -e ${LOG_DIR}/${scenario_id}_%J.err

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
echo "Scenario: ${scenario_id} [v14_netzero, medium cost, central ETS]"
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

    echo "  Submitted: ${scenario_id}  (job name: $job_name)"
done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  bjobs -w"
echo "  bjobs -J 'v14nz-*'"
echo ""
echo "Results: results/supply_curve_v14_netzero/"
echo "Logs:    $LOG_DIR"
