#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CONFIG="$1"
NRUNS="$2"
START="1"
if [ "$#" -eq 3 ]; then
    START="$3"
fi

mkdir -p logs

"${DIR}/scripts/scan2D.py" -c "$CONFIG" -n "$NRUNS" -s "$START"
