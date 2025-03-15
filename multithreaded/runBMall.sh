for i in {0..5}; do
	for j in {0..1}; do
		rootline="'run_B_M_Every.C($i,$j, false)'"
		eval root -l -b -q $rootline
	done
done
