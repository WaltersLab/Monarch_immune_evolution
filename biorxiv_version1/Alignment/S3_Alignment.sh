#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem=20G
#SBATCH -t 24:00:00

## Reference Mapping and BAM file Processing##
# step 1: alignment 
# step 2: add read groups
# step 3: sory BAM file
# step 4: summary stats  
# step 5: remove duplicates
# step 6: index the BAM file 
# step 7: Indel realignment 


# commend line args 1 = sample ID

# Load modules
# picard is under java, so cannot be loaded directly 
module load Bowtie2
module load SAMtools

# Set dirs 
dir_r='/scratch/whtan/Ref'
dir_g='/scratch/whtan/Genomes'
dir_out='/scratch/whtan/Outputs'

## Reference Mapping using Bowtie2
# "very sensitive local" is good for SNP calling
# p = no. core; x = ref; -1&-2 = PE fq files; 
bowtie2 --phred33 \
--very-sensitive-local \
-p 6 \
-x $dir_r/Dp_genome_v3_masked \
-1 $dir_g/$1_1.fastq \
-2 $dir_g/$1_2.fastq \
| samtools view -S -b -o $dir_out/$1.bam 

## Add read groups for the BAM file using picard
# alternatively this can be done with Bowtie2 in alignment step 
# Here, the read groups is arbitrary assigned 
java -Xmx20G -jar /tools/cluster/6.2/picard-tools/1.87/AddOrReplaceReadGroups.jar \
I=$dir_out/$1.bam \
O=$dir_out/$1_rg.bam \
RGID=H0164.2 \
RGLB=library1 \
RGPL=illumina \
RGPU=H0164ALXX140820.2 \
RGSM=$1
# remove previous file
rm $dir_out/$1.bam

## Sort BAM file using picard
# sory by corordinate
java -Xmx20G -jar /tools/cluster/6.2/picard-tools/1.87/SortSam.jar \
INPUT=$dir_out/$1_rg.bam \
OUTPUT=$dir_out/$1_rg_sorted.bam SORT_ORDER=coordinate
# remove previous file
rm $dir_out/$1_rg.bam

## Alignment summary statistics by picard 
java -Xmx20G -jar /tools/cluster/6.2/picard-tools/1.87/CollectAlignmentSummaryMetrics.jar \
R=$dir_r/Dp_genome_v3_masked.fasta \
INPUT=$dir_out/$1_rg_sorted.bam \
OUTPUT=$dir_out/$1_rg_sorted_summary_metrics.txt

## Mark and remove PCR duplicates using picard
# remove = T to remove; otherwise just flagged
# M= a output summary table 
java -Xmx20G -jar /tools/cluster/6.2/picard-tools/1.87/MarkDuplicates.jar \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 \
I=$dir_out/$1_rg_sorted.bam \
REMOVE_DUPLICATES=TRUE \
O=$dir_out/$1_rg_sorted_rmdup.bam \
M=$dir_out/$1_marked_dup_metrics.txt
# remove previous file
rm $dir_out/$1_rg_sorted.bam

## Create index (.bai) for the sorted, cleaned BAM file using picard
java -Xmx20G -jar /tools/cluster/6.2/picard-tools/1.87/BuildBamIndex.jar \
INPUT=$dir_out/$1_rg_sorted_rmdup.bam \
OUTPUT=$dir_out/$1_rg_sorted_rmdup.bam.bai

## Indel realignment using GATK 
# First, need to create a target (.bam.list)
java -Xmx20G -jar /tools/cluster/6.2/gatk/3.4-46/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-I $dir_out/$1_rg_sorted_rmdup.bam \
-R $dir_r/Dp_genome_v3_masked.fasta  \
-o $dir_out/$1_rg_sorted_rmdup.bam.list
# Then, run realignment
java -Xmx20G -jar /tools/cluster/6.2/gatk/3.4-46/GenomeAnalysisTK.jar \
-T IndelRealigner \
-I $dir_out/$1_rg_sorted_rmdup.bam \
-R $dir_r/Dp_genome_v3_masked.fasta  \
-targetIntervals $dir_out/$1_rg_sorted_rmdup.bam.list \
-o $dir_out/$1_realigned.bam
# remove previous file
rm $dir_out/$1_rg_sorted_rmdup.bam 
rm $dir_out/$1_rg_sorted_rmdup.bam.bai
