#!/bin/bash
#SBATCH -o logs/slurm-%x-%j.out

DIR="$1"
CONFIG="$2"
NRUNS="$3"

source "${DIR}/../setup.sh"

"${DIR}/scripts/scan1D.py" "$CONFIG" "$NRUNS"
