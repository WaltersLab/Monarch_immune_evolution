#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=20G
#SBATCH -t 24:00:00

# Set dirs and commend line args
dir=''
cd $dir 
# $1 = input GFF 
# $2 = population 
# $3 = sites 

# run script 
python S4c_generate_input_to_run_angsd_for_all_genes.py $1 $2 $3 