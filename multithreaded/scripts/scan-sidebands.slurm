#!/bin/bash
#SBATCH -o logs/slurm-%x-%j.out

RUN_DIR="$1"
CONFIG="$2"
NRUNS="$3"

source "${RUN_DIR}/../setup.sh"

echo "${INSTALL_DIR}/multithreaded/scripts/scan-sidebands.py" -c "$CONFIG" -n "$NRUNS"
"${INSTALL_DIR}/multithreaded/scripts/scan-sidebands.py" -c "$CONFIG" -n "$NRUNS"
