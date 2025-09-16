#!/bin/bash
#SBATCH -o logs/slurm-%x-%j.out

DIR="$1"
CONFIG="$2"
NRUNS="$3"

source "${DIR}/setup.sh"

"${DIR}/scripts/scanBM.py" "$CONFIG" "$NRUNS"
