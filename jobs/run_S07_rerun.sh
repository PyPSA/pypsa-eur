#!/bin/bash
# =============================================================================
# run_S07_rerun.sh
#
# Re-runs S07 (525 EUR/tCO2) from the 2040 solve onwards.
#
# Background:
#   The S07 2040 solve terminated sub-optimally (Gurobi status 13) because
#   BarHomogeneous=1 caused the Homogeneous Barrier to converge to a
#   pathological solution with ~10,000+ GW DAC and ~120,000 GW bioCCS.
#   These corrupt 2040 capacities became infeasible brownfield lower bounds
#   in 2050, causing the 2050 solve to fail immediately.
#
# Fix:
#   config.S07-525eur-fixed.yaml uses BarHomogeneous=0 (standard barrier).
#   The corrupt 2040 network is deleted so snakemake re-runs from there.
#   --rerun-triggers mtime skips the already-good 2030 solve.
#
# Usage:
#   bash run_S07_rerun.sh
# =============================================================================

set -euo pipefail

REPO_ROOT="/work3/s240459/pypsa-eur-thesis"
RESULT_DIR="$REPO_ROOT/results/S07-cdr-525eur"
LOG_DIR="$REPO_ROOT/logs/supply_curve"
FIXED_CFG="$REPO_ROOT/config/Myruns/supply_curve/config.S07-525eur-fixed.yaml"

WALLTIME="60:00"
CORES=1
MEM_MB=48000
QUEUE="hpc"

echo "Deleting corrupt S07 2040 network..."
rm -f "$RESULT_DIR/networks/base_s_96__24h_2040.nc"
echo "  Done."

echo "Submitting S07 re-run job..."

bsub <<BSUB_SCRIPT
#!/bin/bash
#BSUB -J sc-S07-525eur-rerun
#BSUB -q ${QUEUE}
#BSUB -n ${CORES}
#BSUB -W ${WALLTIME}
#BSUB -R "rusage[mem=${MEM_MB}]"
#BSUB -o ${LOG_DIR}/S07-525eur-rerun_%J.out
#BSUB -e ${LOG_DIR}/S07-525eur-rerun_%J.err

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
echo "Scenario: S07-cdr-525eur (re-run, BarHomogeneous=0)"
echo "Config:   ${FIXED_CFG}"
echo "Started:  \$(date)"

set +e
"\$PIXI_BIN" run snakemake \
  -j ${CORES} \
  --nolock \
  --rerun-incomplete \
  --rerun-triggers mtime \
  --keep-going \
  --printshellcmds \
  --configfile "${FIXED_CFG}" \
  2>&1
status=\$?
set -e

if [[ "\$status" -ne 0 ]]; then
    echo "Snakemake failed with status \$status"
    exit "\$status"
fi

echo "Completed: \$(date)"
BSUB_SCRIPT

echo "Done. Monitor with:  bjobs -J 'sc-S07-525eur-rerun'"
echo "Logs: $LOG_DIR/S07-525eur-rerun_<jobid>.out"
