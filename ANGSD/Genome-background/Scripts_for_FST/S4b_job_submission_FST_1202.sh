#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=20G
#SBATCH -t 24:00:00

# Set dirs and commend line args
dir='/scratch/whtan/ANGSD_Nov2018_GB'
cd $dir 
# $1 = input GFF 
# $2 = population 
# $3 = sites 

# run script 
python S4c_run_angsd_for_all_genes_inputted_FST_1202.py $1 $2 $3 