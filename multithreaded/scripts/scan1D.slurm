#!/bin/bash
#SBATCH -o logs/slurm-%x-%j.out

RUN_DIR="$1"
CONFIG="$2"
NRUNS="$3"
LIMITS="$4"
FIT="$5"
SIGMA="$6"

source "${RUN_DIR}/../setup.sh"

echo "${INSTALL_DIR}/multithreaded/scripts/scan1D.py" -c "$CONFIG" -n "$NRUNS" -i "$LIMITS" -f "$FIT" -s "$SIGMA"
"${INSTALL_DIR}/multithreaded/scripts/scan1D.py" -c "$CONFIG" -n "$NRUNS" -i "$LIMITS" -f "$FIT" -s "$SIGMA"
