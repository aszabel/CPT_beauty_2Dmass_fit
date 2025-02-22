#!/bin/bash
spack load root

t=$(($1-1))
START=${PWD}
for k in `seq 1 1`
do
        m=$(($t*1+$k))

	mkdir -p toy_res/toy_$m
	JSON_FILE="config_inter.json"
	JSON_FILE2="config_unbinned.json"
	JSON_SIDE="config_sidebands.json"
	cp -r root.sh fit2D_mass fit1D_dmass $JSON_FILE $JSON_FILE2 $JSON_SIDE root_test.sh D_M_results B_M_results check_toys.C toy_res/toy_$m
	cd toy_res/toy_$m


	# Define the key to be replaced and the new value
	KEY="input_file"
        OLD_VALUE=ToyMC_BM_DM_xxx.root
        NEW_VALUE=ToyMC_BM_DM_${m}.root



	# Use sed to replace the old value with the new value
	sed -i "s/$OLD_VALUE/$NEW_VALUE/" "$JSON_FILE"
	sed -i "s/$OLD_VALUE/$NEW_VALUE/" "$JSON_FILE2"
	sed -i "s/$OLD_VALUE/$NEW_VALUE/" "$JSON_SIDE"

	
	./fit1D_dmass ${JSON_SIDE}	


	JSON_FILE_UNBINNED="config_unbinned.json"
	JSON_FILE=("config_inter.json" "config_unbinned_best.json")
	cp $JSON_FILE_UNBINNED ${JSON_FILE[1]} 
	HOME=${PWD}
	RESULTS=("results_binned100kMU_2" "results_unbinned100kMU_2")
	Nfits=(30 5)

	for j in `seq 0 1`
	do	
		cd $HOME
		best1=$(eval root -b -q -l 'check_toys.C\(\"${RESULTS[0]}\"\)')
		best=$(echo $best1| awk '{print $NF}')
		OLD_value="fit2D_xxxx"
		NEW_value="fit2D_"$best
		echo TUUUUU
		echo $NEW_value
		sed -i "s/$OLD_value/${NEW_value}/" "${JSON_FILE[$j]}"
		#for i in `seq 1 ${Nfits[$j]}`
		#do
		#	cd $HOME
		#	mkdir -p ${RESULTS[$j]}/fit2D_$i
		#	cp -r root.sh fit2D_mass ${JSON_FILE[$j]} root_test.sh D_M_results B_M_results ${RESULTS[$j]}/fit2D_$i
		#	cd ${RESULTS[$j]}/fit2D_$i


			# Define the key to be replaced and the new value
			KEY="randSeed"
			OLD_VALUE=-1
			NEW_VALUE=1


			# Use sed to replace the old value with the new value
			sed -i "s/\"$KEY\": $OLD_VALUE/\"$KEY\": $NEW_VALUE/" "${JSON_FILE[$j]}"

  			time ./fit2D_mass ${JSON_FILE[$j]}
			#sbatch -p INTEL_HASWELL --time=10:00:00 --cpus-per-task=40 --nodes=1 --job-name="fit2Dfull"  --output="outputs/fit2D_${m}_,${i}.log" root.sh ${JSON_FILE[$j]} 
		#done
	done
	cd $HOME
	best1=$(eval root -b -q -l 'check_toys_final.C\(\"${RESULTS[1]}\"\)')
	best=$(echo $best1| awk '{print $NF}')
	echo "TO TU---->"
	echo $best
	cd $START
done
