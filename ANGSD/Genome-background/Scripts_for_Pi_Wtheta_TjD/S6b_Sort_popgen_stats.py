#!/bin/usr/env python 

### This is a script for parsing all angsd final output files and extract Tajima'sD, Pi, and Watterson's theta statistics ### 

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
output_d = open(str(args.population) + "_tjd_output_" + str(args.site) + ".txt", "w") # for Tajima's D
output_d.write("%s\t%s\n" % ("gene.id", "stat_value"))
output_t = open(str(args.population) + "_wtheta_output_" + str(args.site) + ".txt", "w")  # for watterson's theta
output_t.write("%s\t%s\n" % ("gene.id", "stat_value"))
output_p = open(str(args.population) + "_pi_output_" + str(args.site) + ".txt", "w")  # for Pi
output_p.write("%s\t%s\n" % ("gene.id", "stat_value"))
output_S = open(str(args.population) + "_nSites_output_" + str(args.site) + ".txt", "w")  # for checking the number of sites used (both total and non-NA only)
output_S.write("%s\t%s\t%s\n" % ("gene.id", "total_sites", "non-NA_sites"))

## Extract Tajima's D, theta and pi statistics, and output with gene names 
with open(args.gene_list, "r") as list: # open a list of samples
	for i in list:
		gene = i.rstrip("\r\n") # get gene name and use as input later
		# Tajima's D 
		filename1 = str(gene)+"_"+str(args.site)+"_outFold.thetas.idx.pestPG"
		with open(filename1, "r") as outfile:   # open each angsd final output file
			next(outfile) # skip header 
			TJ_stat = [line.split("\t")[8] for line in outfile]  # here is the Tajima's D value 
			if len(TJ_stat) > 0:   # if there is no output 
				TJ_stat = TJ_stat[0]
			else:
				TJ_stat = "NA"
			output_d.write("%s\t%s\n" % (gene, TJ_stat))  # output
		
		# Watterson's theta and Pi 		 
		filename2 = str(gene)+"_"+str(args.site)+"_outFold.thetas.idx.txt"		
		dat = pd.read_csv(filename2, sep = "\t") # open each angsd theta output file as pd dataframe
		data = dat.copy()  # deep copy the dataframe 	
		if data.shape[0] == 0:   # if there is no output 
			wtheta_m = "NA"
			Pi_m = "NA"
			Pi_len1 = 0
			Pi_len2 = 0
		else: 	
			data.iloc[:, 2:4] = np.exp(data.iloc[:,2:4])  # raw data for pi and wtheta are in natural log scale. Convert back. 
			wtheta = data.iloc[:, 2]  # watterson's theta per-site 
			wtheta_m = np.nanmean(wtheta)  # watterson's theta 
			Pi = data.iloc[:, 3]   # Pi per-site 
			Pi_m = np.nanmean(Pi)   # Pi 
			Pi_len1 = len(Pi)   # the total length of Pi 
			Pi_len2 = np.count_nonzero(~np.isnan(Pi)) # the length of non-NA sites of Pi 
		output_t.write("%s\t%s\n" % (gene, wtheta_m))  # output the mean of the theta (use nanmean in case of nan present)
		output_p.write("%s\t%s\n" % (gene, Pi_m))  # output the mean of the pi (use nanmean in case of nan present)
		output_S.write("%s\t%s\t%s\n" % (gene, Pi_len1, Pi_len2))  # output the number of sites in the file (both total and non-NA only)
list.close()
output_d.close()
output_t.close()
output_p.close()
output_S.close()
