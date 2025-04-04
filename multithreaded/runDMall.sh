#!/bin/bash
#SBATCH --ntasks=1			# Number of processes
#SBATCH --time=24:00:00			# Time limit hrs:min:sec

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source $DIR/setup.sh
mkdir -p results
cd results

for i in {0..5}; do
	for j in {0..1}; do
		rootline="'${DIR}/macros/run_D_M_Every.C($i,$j, \"$1\", false)'"
		echo eval root -l -b -q $rootline
	done
done
