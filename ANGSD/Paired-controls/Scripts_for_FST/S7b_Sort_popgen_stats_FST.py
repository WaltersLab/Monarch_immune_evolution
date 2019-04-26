#!/bin/usr/env python 

### Parsing all angsd final output files and extract Tajima'sD, Pi, and Watterson's theta statistics ### 

from __future__ import division 

## import packages 
import os
os.system("module load Biopython")
import argparse
import numpy as np 
import pandas as pd 

## Set command line arguments 
parser = argparse.ArgumentParser(description = "Calculate per-gene stats from per-site raw data")
parser.add_argument("population", type = str, help = "assgin the population to use")
parser.add_argument("gene_list", type = str, help = "input a list of genes of interest (.txt)")
parser.add_argument("site", type = str, help = "specify either full, 0f, or 4f")
args = parser.parse_args()

## Create output files
# unweighted and weighted Fst 
output_U = open(str(args.population) + "_fst_unw_output_" + str(args.site) + ".txt", "w") 
output_U.write("%s\t%s\n" % ("gene.id", "stat_value"))
output_W = open(str(args.population) + "_fst_w_output_" + str(args.site) + ".txt", "w")
output_W.write("%s\t%s\n" % ("gene.id", "stat_value"))

# Extract Fst estimate, and output with gene names 
with open(args.gene_list, "r") as list:
	for i in list:
		gene = i.rstrip("\r\n") # get gene name and use as input later
		if args.site == "full": 
			filename = str(gene)+"__NAM."+str(args.population)+"_Fst_values.txt"
		else: 
			filename = str(gene)+"_"+str(args.site)+"_NAM."+str(args.population)+"_Fst_values.txt"
		with open(filename, "r") as outfile:
			for line in outfile:
				line = line.rstrip("\r\n").split(" ")  # rm /r/n and split
				Fst_unw = line[0] # Fst unweighted
				Fst_w = line[1]  # Fst weighted
				output_U.write("%s\t%s\n" % (gene, Fst_unw))
				output_W.write("%s\t%s\n" % (gene, Fst_w))
list.close()
output_U.close()
output_W.close()
