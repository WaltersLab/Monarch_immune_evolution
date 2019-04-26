#!/bin/bash

## Creating Indexes for Reference Genome ##
# First, download reference genome from MonarchBase http://monarchbase.umassmed.edu/home.html
# Used the "repeat masked" v3(newest) genome. 
# In this script all indexes needed for Bowtie2 and GATK will be created. 

# Load modules
module load Bowtie2
module load SAMtools

# Set directory
dir_r=''

#Unzip the fq 
gunzip Dp_genome_v3_masked.fasta.gz

# Create index for Bowtie2
bowtie2-build $dir_r/Dp_genome_v3_masked.fasta $dir_r/Dp_genome_v3_masked

# Create dictionary and index for GATK 
java -jar /tools/cluster/6.2/picard-tools/1.87/CreateSequenceDictionary.jar \
R= $dir_r/Dp_genome_v3_masked.fasta \
O= $dir_r/Dp_genome_v3_masked.dict
samtools faidx $dir_r/Dp_genome_v3_masked.fasta
