#!/bin/bash

### This is a script for running ANGSD to get popgen stats with several inputs ### 

# Load angsd 
module load angsd/0.919

# Info of command line args
# $1 = population name for bam.filelist 
# $2 = gene name (ID)
# $3 = "region" 
# $4 = sites (0f or 4f) 

## Set directories 
dir='/scratch/whtan/ANGSD_Nov2018_GB'
dir_out='/scratch/whtan/ANGSD_Nov2018_GB/FST_Output_'$4
dir_r='/scratch/whtan/Ref'
cd $dir_out/

### Running ANGSD to calculate popgen stats 
## SFS pre-generated 

## Step 3: Estimate pairwise Fst
# index the sample so that the same sites are analysed for each population
angsd sites index $dir/Sites/$2"_sites_"$4.txt        # index the individual gene site file 
realSFS fst index $dir_out/"NAM_"$4.saf.idx $dir_out/$1"_"$4.saf.idx -sfs $dir_out/"NAM."$1"_"$4.ml -sites $dir/Sites/$2"_sites_"$4.txt -r $3 -fstout $2"_"$4"_NAM."$1

## Step 4 get the global estimate
realSFS fst stats $2"_"$4"_NAM."$1.fst.idx > $2"_"$4"_NAM."$1"_Fst_values.txt"
	
## sort outputs
cp $2"_"$4"_NAM."$1"_Fst_values.txt" $dir_out/All_stats/$1  # copy out the final result 
mkdir $dir_out/$1/out_$2
mv $2"_"$4"_NAM."$1* $dir_out/$1/out_$2  # sort all files into their own folder 
