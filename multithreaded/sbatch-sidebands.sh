#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CONFIG="$1"
NRUNS="$2"

mkdir -p logs

sbatch -p INTEL_HASWELL,INTEL_CASCADE,INTEL_SKYLAKE -J scanSide-$(basename -s .json "$CONFIG") "${DIR}/scripts/scan-sidebands.slurm" "$DIR" "$CONFIG" "$NRUNS"

