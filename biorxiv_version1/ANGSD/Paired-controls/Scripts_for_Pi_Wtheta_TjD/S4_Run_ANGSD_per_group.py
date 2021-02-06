#!/usr/bin/env python

## This is a scirpt that submits jobs for all groups ## 

# import packages 
import argparse, os 

# Set command line arguments 
parser = argparse.ArgumentParser(description = "submit jobs for all GFF subfiles")
parser.add_argument("angsd_script", type = str, help = "input the job submission script (.sh) for submitting to slurm")
parser.add_argument("group_list", type = str, help = "input a list of gene groups of interest")
parser.add_argument("population", type = str, help = "assgin the population to use")
parser.add_argument("site", type = str, help = "assgin 0x or 4x sites")
args = parser.parse_args()

# Read in a txt file with all GFF subfile names and submit jobs accordingly
with open(args.group_list, "r") as list: 
	for i in list:
		group = i.rstrip("\r\n")
		os.system("sbatch " + args.angsd_script + " " + args.population + " " + group + " " + args.site)
