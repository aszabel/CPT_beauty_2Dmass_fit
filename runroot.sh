#!/bin/bash
t=$(($1-1))
START=${PWD}
k=1
m=$(($t*3+$k))
m=$1
sbatch -p INTEL_HASWELL --time=10:00:00 --cpus-per-task=40 --nodes=1 --job-name="fit2Dfull"  --output="outputs/fit2D_${m}.log" root_test.sh $1 
