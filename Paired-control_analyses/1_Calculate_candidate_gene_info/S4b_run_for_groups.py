#!/usr/bin/env python

## This is a scirpt that loops through all samples specified ## 

# import packages 
import sys, os 

# Inputs 
# commend line arg 1 = sample list (.txt); commend line arg 2 = job script (.py); commend line arg 3 = 1 (first run) or 2 (re-run)

# Create output files 
if sys.argv[3] == "1":   # if first run 
	out_file = open("expand_list.txt", "w")  # a list to track with genes need to re-run 
	out_file.close()
	out_file_2 = open("updated_focal_gene_list.txt", "w")  # a list to update usable immune genes 
	out_file_2.close()
elif sys.argv[3] == "2": # if re-run 
	out_file_3 = open("expand_list_2.txt", "w")  # another list 
	out_file_3.close()

# Read in a txt file with gene names and do jobs accordingly
with open(sys.argv[1], "r") as list: 
	for i in list:
		gene = i.rstrip("\r\n")[1:-1]  # get gene name for input in the next line
		os.system("python " + sys.argv[2] + " " + gene + " cand_genes_for_" + gene + ".txt " + "Dp_geneset_OGS2.gff3 " + sys.argv[3]) # the job
