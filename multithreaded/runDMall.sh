#!/bin/bash -e
#SBATCH --ntasks=1			# Number of processes
#SBATCH --time=24:00:00			# Time limit hrs:min:sec

JSON_FILE="config_1D_DM.json"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source $DIR/setup.sh
mkdir -p results/D_M_results
mkdir -p results/B_M_results

cd results
for sign in muplus muminus
do
	cp $DIR/configs/$JSON_FILE D_M_results/$JSON_FILE.$sign
	# Define the key to be replaced and the new value
	KEY="sign"
	NEW_VALUE=$sign
	# Use sed to replace the old value with the new value
	sed -i "s/\"$KEY\":.*\".*\"/\"$KEY\": \"$NEW_VALUE\"/" "D_M_results/$JSON_FILE.$sign"

	rootline="${DIR}/macros/run_D_M_Every.C(\"D_M_results/$JSON_FILE.$sign\")"
	echo root -l -b -q $rootline
	root -l -b -q $rootline
done