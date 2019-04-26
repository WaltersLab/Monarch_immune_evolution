#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -p bigm
#SBATCH --mem=20G
#SBATCH -t 24:00:00

### Running ANGSD to get popgen stats with several inputs ### 
# run things per pop per gene group. 

# Load angsd 
module load angsd/0.919

## Commend line args 
# $1 = population name for bam.filelist 
# $2 = gene group name 
# $3 = 0f or 4f sites

## Set directories
dir=''
dir_out=''
dir_r=''

### Running ANGSD to calculate popgen stats 
## Step 1 index the SITE file 
# This is needed for doing the -sites filtering
# the "site" text file is pre-generated for each gene group. Each file has all 0x or 4x sites for the given group. 
angsd sites index $dir/Sites_sorted/$2"_group_sites_"$3.txt

## Step 2 First estimate the folded site allele frequency likelihood (folded sfs)
# -r for a specific region; -rf for input a file, which will be useful for miltiple sites. 
# -s for a list of selected sites. Note that when doing -s, adding also -r or -rf is preferred. 
# Did not do any filtering. 
# fold =1 specifies folded sfs (no ansectory info) (but still need to plug in a anc ref) 
angsd -bam $dir/$1#_bam.filelist -sites $dir/Sites_sorted/$2"_group_sites_"$3.txt  \
-rf $dir/Regions_for_sites_sorted/$2"_group_region".txt  \
-doSaf 1 -anc $dir_r/Dp_genome_v3_masked.fasta -GL 1 -P 12 -out $2"_"$3"_outFold" -fold 1 

## step 3 Obtain the maximum likelihood estimate of the SFS
# realSFS in under angsd 
realSFS $2"_"$3"_outFold".saf.idx -P 12 > $2"_"$3"_outFold".sfs

## step 4 Calculate the thetas (folded)

cat $dir/$2"_list_04x".txt | while read gene 
do 
	angsd sites index $dir/Sites_sorted/$gene"_sites_"$3.txt        # index the individual gene site file 
	angsd -bam $dir/$1#_bam.filelist -sites $dir/Sites_sorted/$gene"_sites_"$3.txt \
	-rf $dir/Regions_for_sites_sorted/$gene"_region".txt \
	-out $gene"_"$3"_outFold" -doThetas 1 -doSaf 1 -pest $2"_"$3"_outFold".sfs \
	-anc $dir_r/Dp_genome_v3_masked.fasta -GL 1 -fold 1
	
## step 5 Estimate Tajimas D
	thetaStat do_stat $gene"_"$3"_outFold".thetas.idx 

## step 6 Print out ThetaStat values 
	thetaStat print $gene"_"$3"_outFold".thetas.idx > $gene"_"$3"_outFold".thetas.idx.txt
done
