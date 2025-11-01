#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CONFIG="$1"
NRUNS="$2"
LIMITS=""
FIT=""
SIGMA="1.0"
if [ "$#" -eq 3 ]; then
    LIMITS="$3"
fi
if [ "$#" -eq 4 ]; then
    LIMITS="$3"
    FIT="$4"
fi
if [ "$#" -eq 5 ]; then
    LIMITS="$3"
    FIT="$4"
    SIGMA="$5"
fi

mkdir -p logs

sbatch -p INTEL_HASWELL,INTEL_CASCADE,INTEL_SKYLAKE -J scan1D-$(basename -s .json "$CONFIG") "${DIR}/scripts/scan1D.slurm" "$DIR" "$CONFIG" "$NRUNS" "$LIMITS" "$FIT" "$SIGMA"
