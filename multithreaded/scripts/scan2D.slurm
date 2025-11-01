#!/bin/bash
#SBATCH -o logs/slurm-%x-%j.out

RUN_DIR="$1"
CONFIG="$2"

source "${RUN_DIR}/../setup.sh"

echo "${INSTALL_DIR}/multithreaded/fit2D_mass" $CONFIG
time "${INSTALL_DIR}/multithreaded/fit2D_mass" $CONFIG >& $CONFIG.log
