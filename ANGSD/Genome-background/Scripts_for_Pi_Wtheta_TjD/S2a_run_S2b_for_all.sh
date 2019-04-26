#!/bin/bash 

#### This is for running the SFS script for all groups #### 

for i in "NAM" "FL" "PAC" "ATL"
do
	for j in "0f" "4f"
	do
		sbatch S1b_generate_GB_SFS.sh $i $j
	done
done
