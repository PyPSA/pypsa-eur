#!/bin/bash
# =============================================================================
# run_supply_curve_rerun_2040plus.sh
#
# Re-runs supply-curve scenarios from the 2040 solve onwards.
#
# Default scenarios: S0, S01, ..., S08.
# The script deletes the 2040/2050 solved networks so Snakemake rebuilds the
# downstream summaries and plots from the first unstable planning horizon.
#
# Usage:
#   bash run_supply_curve_rerun_2040plus.sh
#   bash run_supply_curve_rerun_2040plus.sh S0 S03 S07
#
# Notes:
#   - Uses the scenario configs in config/Myruns/supply_curve/.
#   - Skips any *-fixed.yaml helper configs to avoid duplicate submissions.
#   - Assumes 2030 outputs are kept.
# =============================================================================

set -euo pipefail
shopt -s nullglob

REPO_ROOT="/work3/s240459/pypsa-eur-thesis"
CONFIG_DIR="$REPO_ROOT/config/Myruns/supply_curve"
LOG_DIR="$REPO_ROOT/logs/supply_curve"
mkdir -p "$LOG_DIR"

WALLTIME="72:00"
CORES=1
MEM_MB=48000
QUEUE="hpc"

DEFAULT_PREFIXES=(S0 S01 S02 S03 S04 S05 S06 S07 S08 S09 S10)

if [[ $# -gt 0 ]]; then
    PREFIXES=("$@")
else
    PREFIXES=("${DEFAULT_PREFIXES[@]}")
fi

CONFIGS=()
for prefix in "${PREFIXES[@]}"; do
    before_count=${#CONFIGS[@]}
    match=("$CONFIG_DIR"/config.${prefix}.yaml "$CONFIG_DIR"/config.${prefix}-*.yaml)
    for candidate in "${match[@]}"; do
        [[ -f "$candidate" && "$candidate" != *-fixed.yaml ]] && CONFIGS+=("$candidate")
    done
    if [[ ${#CONFIGS[@]} -eq $before_count ]]; then
        echo "Warning: no config found matching prefix '$prefix' — skipping"
    fi
done

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
    echo "Error: no rerun configs found"
    exit 1
fi

echo "Preparing ${#CONFIGS[@]} scenario(s) for 2040+ rerun..."
echo ""

for cfg in "${CONFIGS[@]}"; do
    basename_cfg=$(basename "$cfg" .yaml)
    scenario_id="${basename_cfg#config.}"
    job_name="sc-${scenario_id}-r2040"
    run_name=$(grep -m1 'name:' "$cfg" | awk '{print $2}' | tr -d '"')
    result_dir="$REPO_ROOT/results/$run_name"

    if [[ ! -d "$result_dir" ]]; then
        echo "Warning: result directory missing for $scenario_id ($result_dir)"
        echo "         Submitting anyway; Snakemake will build what it can."
    else
        echo "Cleaning stale 2040/2050 networks for $scenario_id ..."
        rm -f           "$result_dir/networks/base_s_96__72seg_2040.nc"           "$result_dir/networks/base_s_96__72seg_2050.nc"
    fi

    bsub <<BSUB_SCRIPT
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -q ${QUEUE}
#BSUB -n ${CORES}
#BSUB -W ${WALLTIME}
#BSUB -R "rusage[mem=${MEM_MB}]"
#BSUB -o ${LOG_DIR}/${scenario_id}-r2040_%J.out
#BSUB -e ${LOG_DIR}/${scenario_id}-r2040_%J.err

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
echo "Scenario: ${scenario_id}"
echo "Config:   ${cfg}"
echo "Started:  \$(date)"

set +e
"\$PIXI_BIN" run snakemake   -j ${CORES}   --nolock   --rerun-incomplete   --rerun-triggers mtime params   --keep-going   --printshellcmds   --configfile "${cfg}"   2>&1
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
echo "Monitor with:"
echo "  bjobs -J 'sc-*-r2040'"
echo "Logs: $LOG_DIR"
