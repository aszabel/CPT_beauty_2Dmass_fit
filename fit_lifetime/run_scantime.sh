/bin/bash

for i in `seq 0 200`
do	
	mkdir -p trees
	JSON=trees/config_temp.json
	cp configs/config_inter.json $JSON

	m=$(echo "-0.01+${i}*0.0001"|bc -l| awk '{printf "%.4f", $1}')
	echo $m
	# Define the key to be replaced and the new value
   	KEY="muPTmin"
   	OLD_VALUE=500.0
   	NEW_VALUE=$m
	# Use sed to replace the old value with the new value
	sed -i "s/\"$KEY\": $OLD_VALUE/\"$KEY\": $NEW_VALUE/" "${JSON}"
	
	# Define the key to be replaced and the new value
   	KEY="randSeed"
   	OLD_VALUE=1
   	NEW_VALUE=$i
	# Use sed to replace the old value with the new value
	sed -i "s/\"$KEY\": $OLD_VALUE/\"$KEY\": $NEW_VALUE/" "${JSON}"
	./fit2time $JSON
done

