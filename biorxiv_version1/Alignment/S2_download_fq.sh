#!/bin/bash 

## Download genome fq files by accession numbers##

# commend line args 1 = sample ID

# Load module
module load SRA-Toolkit

# Set directory 
dir_g=''

# Download genome fq from NCBI 
fastq-dump --split-files $1 
