#!/bin/bash
#SBATCH -o logs/slurm-%x-%j.out

CONFIG="$1"
NRUNS="$2"

source "../setup.sh"

echo "${INSTALL_DIR}/multithreaded/scripts/scan1D.py" -c "$CONFIG" -n "$NRUNS"
"${INSTALL_DIR}/multithreaded/scripts/scan1D.py" -c "$CONFIG" -n "$NRUNS"
