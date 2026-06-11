#!/bin/bash
# =============================================================================
# run_supply_curve_v12.sh
# Submit one LSF job per supply curve scenario — v12 = v10-medium + gas-for-industry CC cap (1.0 GW/yr).
# Results go to: results/supply_curve_v12/<scenario>/
#
# Usage (from repo root):
#   bash jobs/run_supply_curve_v12.sh           # submit all scenarios
#   bash jobs/run_supply_curve_v12.sh S01 S03   # submit specific scenarios by prefix
#
# Flags:
#   --clean   Delete existing results for selected scenarios before submitting
#
# Monitor jobs:
#   bjobs -w
#   bjobs -J 'v12-*'
# =============================================================================

set -euo pipefail
shopt -s nullglob

REPO_ROOT="/work3/s240459/pypsa-eur-thesis"
CONFIG_DIR="$REPO_ROOT/config/Myruns/supply_curve_v12/supply_curve_v12_medium"
LOG_DIR="$REPO_ROOT/logs/supply_curve_v12"
mkdir -p "$LOG_DIR"

WALLTIME="72:00"
CORES=8
MEM_MB=15000
QUEUE="hpc"

CLEAN=0
POSITIONAL=()
for arg in "$@"; do
    if [[ "$arg" == "--clean" ]]; then
        CLEAN=1
    else
        POSITIONAL+=("$arg")
    fi
done
set -- "${POSITIONAL[@]+"${POSITIONAL[@]}"}"

if [[ $# -gt 0 ]]; then
    CONFIGS=()
    for prefix in "$@"; do
        before_count=${#CONFIGS[@]}
        match=("$CONFIG_DIR"/config.${prefix}.yaml "$CONFIG_DIR"/config.${prefix}-*.yaml)
        for candidate in "${match[@]}"; do
            [[ -f "$candidate" && "$candidate" != *-fixed.yaml ]] && CONFIGS+=("$candidate")
        done
        if [[ ${#CONFIGS[@]} -eq $before_count ]]; then
            echo "Warning: no config found matching prefix '$prefix' — skipping"
        fi
    done
else
    CONFIGS=()
    for candidate in "$CONFIG_DIR"/config.S*.yaml; do
        [[ -f "$candidate" && "$candidate" != *-fixed.yaml ]] && CONFIGS+=("$candidate")
    done
fi

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
    echo "Error: no configs found in $CONFIG_DIR"
    exit 1
fi

RESULTS_TO_CLEAN=()
for cfg in "${CONFIGS[@]}"; do
    run_name=$(grep -m1 'name:' "$cfg" | awk '{print $2}' | tr -d '"')
    run_prefix=$(grep -m1 'prefix:' "$cfg" | awk '{print $2}' | tr -d '"')
    if [[ -n "$run_prefix" && -n "$run_name" ]]; then
        result_dir="$REPO_ROOT/results/$run_prefix/$run_name"
    elif [[ -n "$run_name" ]]; then
        result_dir="$REPO_ROOT/results/$run_name"
    fi
    if [[ -n "${result_dir:-}" && -d "$result_dir" ]]; then
        RESULTS_TO_CLEAN+=("$result_dir")
    fi
done

if [[ "$CLEAN" -eq 1 ]]; then
    if [[ ${#RESULTS_TO_CLEAN[@]} -eq 0 ]]; then
        echo "No existing result directories found — nothing to clean."
    else
        echo "The following result directories will be deleted:"
        for d in "${RESULTS_TO_CLEAN[@]}"; do
            echo "  $d  ($(du -sh "$d" 2>/dev/null | cut -f1))"
        done
        echo ""
        read -r -p "Confirm deletion? [y/N] " confirm
        if [[ "$confirm" =~ ^[Yy]$ ]]; then
            for d in "${RESULTS_TO_CLEAN[@]}"; do
                echo "  Deleting $d ..."
                rm -rf "$d"
            done
            echo "Done."
        else
            echo "Aborted."
            exit 0
        fi
    fi
fi

echo ""
echo "Submitting ${#CONFIGS[@]} supply curve scenario(s) [v12 = v10-medium + gas-for-industry CC cap (1.0 GW/yr)]..."
echo ""

for cfg in "${CONFIGS[@]}"; do
    basename_cfg=$(basename "$cfg" .yaml)
    scenario_id="${basename_cfg#config.}"
    job_name="v12-${scenario_id}"

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
echo "Scenario: ${scenario_id} [v12 = v10-medium + gas-for-industry CC cap (1.0 GW/yr)]"
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

    echo "  Submitted: $scenario_id  (job name: $job_name)"
done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  bjobs -w"
echo "  bjobs -J 'v12-*'"
echo ""
echo "Results: results/supply_curve_v12/"
echo "Logs:    $LOG_DIR"
