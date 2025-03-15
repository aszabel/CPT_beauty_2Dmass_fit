#!/bin/bash
bad=(7 61 75 131 161 171 202 204 235 246 251 264 266 293 353 373 380 406 416 425 496)



for i in ${bad[@]}
do
	./runroot.sh $i
	#./root_test.sh $i
done
	
