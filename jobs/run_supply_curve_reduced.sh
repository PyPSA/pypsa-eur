#!/bin/bash
# =============================================================================
# Submit reduced-resolution supply-curve scenarios.
#
# Defaults:
#   - configs: config/Myruns/supply_curve_reduced/config.S*.yaml
#   - target: make_cumulative_costs
#   - resolution encoded in configs: 41 clusters x 24seg
#
# The default target intentionally avoids Snakemake's plotting-heavy `all`
# target. It builds solved networks plus summary CSVs needed for the supply
# curve, but skips maps and interactive plots that caused memory pressure.
#
# Usage:
#   bash run_supply_curve_reduced.sh
#   bash run_supply_curve_reduced.sh S0 S05 S10
#   TARGET=solve_sector_networks bash run_supply_curve_reduced.sh S0
# =============================================================================

set -euo pipefail
shopt -s nullglob

REPO_ROOT="/work3/s240459/pypsa-eur-thesis"
CONFIG_DIR="$REPO_ROOT/config/Myruns/supply_curve_reduced"
LOG_DIR="$REPO_ROOT/logs/supply_curve_reduced"
mkdir -p "$LOG_DIR"

WALLTIME="${WALLTIME:-36:00}"
CORES="${CORES:-1}"
MEM_MB="${MEM_MB:-48000}"  # per core; 1 core -> 48 GB total (rusage[mem] is per-slot)
QUEUE="${QUEUE:-hpc}"
TARGET="${TARGET:-make_cumulative_costs}"
TOTAL_MEM_MB=$((CORES * MEM_MB))

if [[ ! -d "$CONFIG_DIR" ]] || ! compgen -G "$CONFIG_DIR/config.S*.yaml" >/dev/null; then
    echo "Reduced configs missing. Generate them with:"
    echo "  python config/Myruns/supply_curve_reduced/generate_reduced_configs.py"
    exit 2
fi

if [[ $# -gt 0 ]]; then
    PREFIXES=("$@")
else
    PREFIXES=(S0 S01 S02 S03 S04 S05 S06 S07 S08 S09 S10)
fi

CONFIGS=()
for prefix in "${PREFIXES[@]}"; do
    before_count=${#CONFIGS[@]}
    match=("$CONFIG_DIR"/config."$prefix".yaml "$CONFIG_DIR"/config."$prefix"-*.yaml)
    for candidate in "${match[@]}"; do
        [[ -f "$candidate" ]] && CONFIGS+=("$candidate")
    done
    if [[ ${#CONFIGS[@]} -eq $before_count ]]; then
        echo "Warning: no reduced config found for prefix '$prefix'"
    fi
done

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
    echo "Error: no reduced configs selected"
    exit 1
fi

echo "Submitting ${#CONFIGS[@]} reduced supply-curve scenario(s)"
echo "Target:   $TARGET"
echo "Walltime: $WALLTIME"
echo "Memory:   ${MEM_MB} MB/core x ${CORES} cores = ${TOTAL_MEM_MB} MB"
echo ""

for cfg in "${CONFIGS[@]}"; do
    basename_cfg=$(basename "$cfg" .yaml)
    scenario_id="${basename_cfg#config.}"
    job_name="scr-${scenario_id}"

    bsub <<BSUB_SCRIPT
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -q ${QUEUE}
#BSUB -n ${CORES}
#BSUB -W ${WALLTIME}
#BSUB -R "rusage[mem=${MEM_MB}]"
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
echo "Scenario: ${scenario_id}"
echo "Config:   ${cfg}"
echo "Target:   ${TARGET}"
echo "Started:  \$(date)"

set +e
"\$PIXI_BIN" run snakemake \
  "${TARGET}" \
  -j ${CORES} \
  --nolock \
  --rerun-incomplete \
  --rerun-triggers params input code \
  --resources mem_mb=${TOTAL_MEM_MB} \
  --printshellcmds \
  --configfile "${cfg}" \
  2>&1
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
echo "  bjobs -w"
echo "  bjobs -J 'scr-*'"
echo "Logs: $LOG_DIR"
