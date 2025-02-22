#!/bin/bash

for i in `seq 0 4`
do	
	mkdir -p trees
	JSON=config_temp.json
	cp config_inter.json $JSON

	# Define the key to be replaced and the new value
   	KEY="randSeed"
   	OLD_VALUE=0
   	NEW_VALUE=$i


	# Use sed to replace the old value with the new value
	sed -i "s/\"$KEY\": $OLD_VALUE/\"$KEY\": $NEW_VALUE/" "${JSON}"
	./fit2time $JSON
done

