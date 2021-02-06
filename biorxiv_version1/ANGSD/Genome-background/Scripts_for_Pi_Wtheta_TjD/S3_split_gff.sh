#!/usr/bin/bash

#### This is a script for splitting the annotation file  ####

# Specify input dataset and number of subsets
file=Dp_geneset_OGS2.gff3
num_files=8 

# Get number of lines 
total=$(wc -l <${file})
per_file=$((total/num_files))

# Split the file 
split --lines=${per_file} ${file} --numeric-suffixes=1 --additional-suffix=.gff3 GFF

# Check 
wc -l GFF*
