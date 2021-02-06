#!/bin/bash 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -c 16
#SBATCH --mem=30G 
#SBATCH -t 48:00:00 

### This is a script for running ANGSD to get SFS estimated using genome inforamation as prior ### 
# exclude genes /sites on the Z-chr 

# Load angsd 
module load angsd/0.919

## Commend line args 
# $1 = 0f or 4f sites

## Set directories
dir=''
dir_r='' # reference

### Running ANGSD to calculate popgen stats 
## Step 1 index the SITE file 
# This is needed for doing the -sites filtering
# the "site" text file is pre-generated for each gene group. Each file has all 0x or 4x sites for the given group. 
angsd sites index $dir/Sites/"All_"$1"_sites_nonZ.txt"

# Step 2: Reconstructing the site frequency spectrum
# calculate per pop saf
# no filtering 
angsd -bam $dir/NAM#_bam.filelist -sites $dir/Sites/"All_"$1"_sites_nonZ.txt" -anc $dir_r/Dp_genome_v3_masked.fasta -out "NAM_"$1 -dosaf 1 -gl 1 
angsd -bam $dir/FL#_bam.filelist -sites $dir/Sites/"All_"$1"_sites_nonZ.txt" -anc $dir_r/Dp_genome_v3_masked.fasta -out "FL_"$1 -dosaf 1 -gl 1 
angsd -bam $dir/PAC#_bam.filelist -sites $dir/Sites/"All_"$1"_sites_nonZ.txt" -anc $dir_r/Dp_genome_v3_masked.fasta -out "PAC_"$1 -dosaf 1 -gl 1 
angsd -bam $dir/ATL#_bam.filelist -sites $dir/Sites/"All_"$1"_sites_nonZ.txt" -anc $dir_r/Dp_genome_v3_masked.fasta -out "ATL_"$1 -dosaf 1 -gl 1 

# Step 3: Calculate the 2-d sfs for prior
realSFS "NAM_"$1.saf.idx "FL_"$1.saf.idx > "NAM.FL_"$1.ml
realSFS "NAM_"$1.saf.idx "PAC_"$1.saf.idx > "NAM.PAC_"$1.ml
realSFS "NAM_"$1.saf.idx "ATL_"$1.saf.idx > "NAM.ATL_"$1.ml
