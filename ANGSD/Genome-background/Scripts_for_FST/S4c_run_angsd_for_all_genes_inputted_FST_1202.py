#!bin/usr/env python 

### Run ANGSD for measuring popgen stats for all genes (i.e. genic region) in the genome background.  ### 
## Parse the annotation GFF file (OGS 2.0) to get positions for each gene. 
## modified from the script made for extracting exon positions. 
## Need to have angsd bam file list set up beforehand. Do one population a time. 
## Summary: Input GFF file, extract gene ID and positions, assgin them to angsd, run angsd (per gene), and get all angsd outputs with gene names. 

# Import packages 
import argparse, os, re

# Set command line arguments 
parser = argparse.ArgumentParser(description = "get info for all genes and run angsd")
parser.add_argument("annotation_file", type = str, help = "input the annotation file (.GFF3 file) for OSG2.0")
parser.add_argument("population", type = str, help = "assgin the population to use")
parser.add_argument("sites", type = str, help = "assgin the sites to run. either 0f or 4f")
args = parser.parse_args()

# Read in a list of all non-Z genes in the genome 
nonZ_list = open("All_aut_ID_list.txt", "r").readlines()   # input and read as list 
nonZ_list = [i.rstrip("\r\n") for i in nonZ_list]  # remove \r\n 
nonZ_list = [i.strip('\"') for i in nonZ_list]   # remove quota 

# Parsing the GFF file to get info for all 15130 genes
with open(args.annotation_file, "r") as data:   # open GFF file
# work through lines to get info for each gene
		for line in data.readlines():
			record = line.rstrip("\r\n").split("\t")  # remove "\r\n", split lines by tab
			if record[2] == "gene": # get the line that contains whole gene info
				ID = re.split("=|;", record[8])[3]   # get the gene ID
				# only do genes that are not on Z-chr 
				if ID in nonZ_list:				
					region = str(record[0] + ":" + record[3] + "-" + record[4])  # the "region" input for angsd. Format: scaffold:start-end 
					# Run angsd for each gene, with ID and region as input 
					os.system("sh S4d_ANGSD_stats_FST_1202.sh " + args.population + " " + str(ID) + " " + str(region) + " " + args.sites) 
