#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p bigm
#SBATCH --mem=20G
#SBATCH -t 3:00:00

# Commend line args 
# $1 = population 

# load module 
module load Biopython

dir=''

for i in "full" "0f" "4f"
do
	cp $dir/"S7b_Sort_popgen_stats.py" $dir/Output_$i/All_stats/$1/
	if [ $i == "full" ]
	then 
		gene_list="Gene_List_full.txt"
		cp $dir/"Gene_List_full.txt" $dir/Output_$i/All_stats/$1/
	else
		gene_list="Gene_List_04x.txt"
		cp $dir/"Gene_List_04x.txt" $dir/Output_$i/All_stats/$1/	
	fi
	
	cd $dir/Output_$i/All_stats/$1/
	python S7b_Sort_popgen_stats.py $1 $gene_list $i
	cp *_output_* $dir/Output_$i/
	mv *_output_* $dir/RESULTS
	rm $dir/Output_$i/All_stats/$1/"S7b_Sort_popgen_statspy"
	rm $dir/Output_$i/All_stats/$1/$gene_list
	
done
