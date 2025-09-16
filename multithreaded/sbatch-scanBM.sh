#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CONFIG="$1"
NRUNS="$2"

mkdir -p logs

sbatch -p INTEL_HASWELL,INTEL_IVY -J scan-$(basename -s .json "$CONFIG") "${DIR}/scripts/scanBM.slurm" "$DIR" "$CONFIG" "$NRUNS"
