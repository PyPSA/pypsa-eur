#!/bin/bash
#BSUB -J jupyter
#BSUB -q hpc
#BSUB -n 6
#BSUB -R "span[hosts=1]"
#BSUB -W 8:00
#BSUB -R "rusage[mem=8000]"
#BSUB -o logs/jupyter.%J.out
#BSUB -e logs/jupyter.%J.err

set -euo pipefail

cd /work3/s240459/pypsa-eur-thesis

export PATH="$HOME/.pixi/bin:$HOME/.local/bin:$PATH"

PORT=8889
NODE=$(hostname)

echo "============================================================"
echo "Jupyter is running on compute node: $NODE"
echo "Port: $PORT"
echo "------------------------------------------------------------"
echo "On your LOCAL machine, open the SSH tunnel:"
echo "  ssh -N -L ${PORT}:${NODE}:${PORT} s240026@login.hpc.dtu.dk"
echo "Then open the http://127.0.0.1:${PORT}/... URL printed below."
echo "============================================================"

exec pixi run jupyter lab --no-browser --ip=0.0.0.0 --port=${PORT}
