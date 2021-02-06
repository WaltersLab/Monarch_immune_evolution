#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH -t 3:00:00
#SBATCH --mail-user=wtan4@emory.edu
#SBATCH --mail-type=ALL

# Commend line args 
# $1 = population 

module load Biopython

dir='/scratch/whtan/ANGSD_Nov2018_GB'
cd $dir

for i in "0f" "4f"
do
	cp $dir/"S6b_Extract_stats_new_FST_1203.py" $dir/FST_Output_$i/All_stats/$1/	
	gene_list="gene_list_noZ.txt"
	cp $dir/$gene_list $dir/FST_Output_$i/All_stats/$1/
	cd $dir/FST_Output_$i/All_stats/$1/
	python S6b_Extract_stats_new_FST_1203.py $1 $gene_list $i
	cp *_output_* $dir/FST_Output_$i/
	mv *_output_* $dir/RESULTS
	rm $dir/FST_Output_$i/All_stats/$1/"S6b_Extract_stats_new_FST_1203.py"
	rm $dir/FST_Output_$i/All_stats/$1/$gene_list
done
