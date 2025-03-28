#!/bin/bash
mkdir -p results/fit2D_$1
JSON_FILE="config_example.json"
cp -r fit2D_mass configs/config_example.json root.sh D_M_results B_M_results results/fit2D_$1
cd results/fit2D_$1


# Define the key to be replaced and the new value
KEY="randSeed"
OLD_VALUE=-1
NEW_VALUE=$1


# Use sed to replace the old value with the new value
sed -i "s/\"$KEY\": $OLD_VALUE/\"$KEY\": $NEW_VALUE/" "$JSON_FILE"

sbatch -p INTEL_CASCADE --time=5:00:00 --cpus-per-task=40 --nodes=1 --job-name="fit2Dfull"  --output="fit2D.log" root.sh $JSON_FILE 
