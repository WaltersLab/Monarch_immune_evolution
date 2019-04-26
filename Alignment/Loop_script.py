#!/usr/bin/env python

## This is a scirpt that loops through all samples specified ## 

# import packages 
import sys, os 

# Inputs 
# commend line arg 1 = sample list (.txt); commend line arg 2 = job script (.sh)

# Read in a txt file with sample IDs/names and sutmit jobs accordingly
with open(sys.argv[1], "r") as list: 
	for i in list:
		ID = i.rstrip("\n")
		os.system("sbatch " + sys.argv[2] + " " + ID)

