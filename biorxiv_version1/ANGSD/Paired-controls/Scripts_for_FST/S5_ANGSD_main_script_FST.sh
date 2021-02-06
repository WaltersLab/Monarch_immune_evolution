#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -p bigm
#SBATCH --mem=20G
#SBATCH -t 24:00:00

### Running ANGSD to get popgen stats with several inputs ### 
# FST (NAM as ref) 
# run things per gene group. 

# Load angsd 
module load angsd/0.919

## Commend line args 
# $1 = gene group name 
# $2 = 0f or 4f sites

## Set directories
dir=''
dir_out=''
dir_r=''

### Running ANGSD to calculate popgen stats 
## Pair-wise FST between populations  -- comparisons in pair, but can do multiple populations in one run (here, doing 4 pops: NAM, FL, PAC, ATL)

## Step 1 index the SITE file 
# This is needed for doing the -sites filtering
# the "site" text file is pre-generated for each gene group. Each file has all 0x or 4x sites for the given group. 
angsd sites index $dir/Sites_sorted/$1"_group_sites_"$2.txt

# Step 2: Reconstructing the site frequency spectrum
# calculate per pop saf
# no filtering 
angsd -bam $dir/NAM#_bam.filelist -sites $dir/Sites_sorted/$1"_group_sites_"$2.txt -rf $dir/Regions_for_sites_sorted_FST/$1"_group_region".txt -anc $dir_r/Dp_genome_v3_masked.fasta -out $1"_NAM" -dosaf 1 -gl 1 
angsd -bam $dir/FL#_bam.filelist -sites $dir/Sites_sorted/$1"_group_sites_"$2.txt -rf $dir/Regions_for_sites_sorted_FST/$1"_group_region".txt -anc $dir_r/Dp_genome_v3_masked.fasta -out $1"_FL" -dosaf 1 -gl 1 
angsd -bam $dir/PAC#_bam.filelist -sites $dir/Sites_sorted/$1"_group_sites_"$2.txt -rf $dir/Regions_for_sites_sorted_FST/$1"_group_region".txt -anc $dir_r/Dp_genome_v3_masked.fasta -out $1"_PAC" -dosaf 1 -gl 1 
angsd -bam $dir/ATL#_bam.filelist -sites $dir/Sites_sorted/$1"_group_sites_"$2.txt -rf $dir/Regions_for_sites_sorted_FST/$1"_group_region".txt -anc $dir_r/Dp_genome_v3_masked.fasta -out $1"_ATL" -dosaf 1 -gl 1 

# Step 3: Calculate the 2-d sfs for prior
realSFS $1"_NAM".saf.idx $1"_FL".saf.idx > $1"_NAM.FL".ml
realSFS $1"_NAM".saf.idx $1"_PAC".saf.idx > $1"_NAM.PAC".ml
realSFS $1"_NAM".saf.idx $1"_ATL".saf.idx > $1"_NAM.ATL".ml

# do Step 4 &5 in a loop 
cat $dir/$1"_list_04x".txt | while read gene 
do 
	angsd sites index $dir/Sites_sorted/$gene"_sites_"$2.txt        # index the individual gene site file 
	reg=$(cat $dir/Regions_for_sites_sorted_FST/$gene"_region".txt)
	
	for POP in "FL" "PAC" "ATL"
	do 
		cd $dir_out/TEMP_$POP
		
		# Step 3: Estimate pairwise Fst
		# index the sample so that the same sites are analysed for each population
		realSFS fst index $dir_out/$1"_NAM".saf.idx $dir_out/$1"_"$POP.saf.idx -sfs $dir_out/$1"_NAM."$POP.ml -sites $dir/Sites_sorted/$gene"_sites_"$2.txt -r $reg -fstout $gene"_"$2"_NAM."$POP

		# Step 4 get the global estimate
		realSFS fst stats $gene"_"$2"_NAM."$POP.fst.idx > $gene"_"$2"_NAM."$POP"_Fst_values.txt"
	
		# Sort outputs
		cp $gene"_"$2"_NAM."$POP"_Fst_values.txt" $dir_out/All_stats/$POP  # copy out the final result 
		mkdir $dir_out/$POP/out_$gene/  # sort all files into their own folder 
		mv $gene"_"$2"_NAM."$POP* $dir_out/$POP/out_$gene/
	done
done 
