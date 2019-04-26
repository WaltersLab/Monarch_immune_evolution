#!/bin/bash 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -c 12 
#SBATCH --mem=30G 
#SBATCH -t 24:00:00 

#### This is a script for running ANGSD to get SFS estimated using genome inforamation as prior #### 
# exclude genes /sites on the Z-chr 

# Load angsd 
module load angsd/0.919

## Commend line args 
# $1 = population name for bam.filelist 
# $2 = 0f or 4f sites

## Set directories
dir=''
dir_r=''  # reference

### Running ANGSD to calculate popgen stats 
## Step 1 index the SITE file 
# This is needed for doing the -sites filtering
# the "site" text file is pre-generated for each gene group. Each file has all 0x or 4x sites for the given group. 
angsd sites index $dir/Sites/"All_"$2"_sites_nonZ.txt"

## Step 2 First estimate the folded site allele frequency likelihood (folded sfs)
# -r for a specific region; -rf for input a file, which will be useful for miltiple sites. 
# -s for a list of selected sites. Note that when doing -s, adding also -r or -rf is preferred. 
# Did not do any filtering. 
# fold =1 specifies folded sfs (no ansectory info) (but still need to plug in a anc ref) 
angsd -bam $dir/$1#_bam.filelist -sites $dir/Sites/"All_"$2"_sites_nonZ.txt"  \
-doSaf 1 -anc $dir_r/Dp_genome_v3_masked.fasta -GL 1 -P 12 -out $1"_"$2"_outFold" -fold 1 

## step 3 Obtain the maximum likelihood estimate of the SFS
# realSFS in under angsd 
realSFS $1"_"$2"_outFold".saf.idx -P 12 > $1"_"$2"_outFold".sfs
