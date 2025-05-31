#!/bin/bash -e
#SBATCH --ntasks=1			# Number of processes
#SBATCH --time=24:00:00			# Time limit hrs:min:sec

#JSON_FILE="config_1D_DM_DCBplusGaus.json"
#JSON_FILE="config_draw.json"
JSON_FILE="$1"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#source $DIR/setup.sh

cd results
cp $DIR/configs/$JSON_FILE .

rootline="${DIR}/macros/runDraw.C(\"$JSON_FILE\")"
#echo root -l -b -q $rootline
root -l -b -q $rootline
