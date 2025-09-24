#!/bin/bash
#SBATCH -o logs/slurm-%x-%j.out

CONFIG="$1"
NRUNS="$2"

source "../setup.sh"

"${INSTALL_DIR}/multithreaded/scripts/scan1D.py" -c "$CONFIG" -n "$NRUNS"
