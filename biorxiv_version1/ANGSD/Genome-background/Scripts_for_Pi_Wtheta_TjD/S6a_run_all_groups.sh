#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH -t 3:00:00

#### Organizing all output files ####

## Commend line args 
# $1 = population 

## load module and set dir 
module load Biopython
dir=''

## do for all groups 
for i in "0f" "4f"
do
	cp $dir/"S6b_Sort_popgen_stats.py" $dir/Output_$i/All_stats/$1/
	gene_list="gene_list_noZ.txt"
	cp $dir/$gene_list $dir/Output_$i/All_stats/$1/		
	cd $dir/Output_$i/All_stats/$1/
	python S6b_Sort_popgen_stats.py $1 $gene_list $i
	cp *_output_* $dir/Output_$i/
	mv *_output_* $dir/RESULTS
	rm $dir/Output_$i/All_stats/$1/"S6b_Sort_popgen_stats.py"
	rm $dir/Output_$i/All_stats/$1/$gene_list
done
