#!bin/usr/env python 

#### Parsing the annocation GFF file (OGS 2.0) to get lengths and distances to focal immune gene of each candidate control genes ####
# Input one focal (immune) gene, and a dataframe of candidate genes, with the monarch annotation GFF file, get the length and distance and output as a txt file.

# Import packages 
import re, argparse, os
import numpy as np

# Set command line arguments 
parser = argparse.ArgumentParser(description = "get info for each candidate control gene")
parser.add_argument("immune_gene", type = str, help = "input the name of the focal immune gene")
parser.add_argument("cand_gene_list", type = str, help = "input a dataframe (.txt) that contains info for candidate genes of interest")
parser.add_argument("annotation_file", type = str, help = "input the annotation file (.GFF3 file) for OSG2.0")
parser.add_argument("run_number", type = str, help = "1 = first run, 2 = re-run")
args = parser.parse_args()

# Create an output file 
# note that the "region" in the output is for generating range files used for ANGSD later
out_file = open(str(args.immune_gene+"_all_cands_calc.txt"), "w")  #create an output file
out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("name", "scaffold", "start", "end", "region", "length", "distance", "reverse", "length_ok", "dis_ok", "main"))  #write out headers  

# Parsing the GFF file to get info for the focal immune gene 
with open(args.annotation_file, "r") as data:   # open GFF file
		#working through lines 
		for line in data.readlines():
			line = line.rstrip("\r\n")  # remove "\r\n"
			record = line.split("\t")   # split lines by tab 
			ID = re.split("=|;", record[8])[3]   # get the ID
			if ID == args.immune_gene and record[2] == "gene":  # get the line that contains whole gene info
				IG_scaf = record[0] # scaffold info
				IG_start = record[3]  # start loc. for that gene
				IG_end = record[4]  # end loc. for that gene 
				IG_length = abs(int(IG_end) - int(IG_start)) # calc. the length of the gene
				IG_region = str(IG_scaf+":"+IG_start+"-"+IG_end)  # this is for "region" input for ANGSD. format = scaf:start-end, e.g., DPSCF300010:1000-3000
				if record[6] == "-":  # record the direction of that gene. if '-', reverse complement. 
					IG_rev = "rc"
				else:
					IG_rev = "f"
out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (args.immune_gene, IG_scaf, IG_start, IG_end, IG_region, IG_length, 0, IG_rev, "N/A", "N/A", 1))  #write out 
data.close()								

# Parsing the GFF file to get info and calulate length and distances for each candidate gene 
K = 0 # for tracking number of candidates 
with open(args.cand_gene_list, "r") as gene_list: # input cand. gene info 
	next(gene_list)   # skip header 
	for gene in gene_list:
		gene_name = gene.rstrip('\r\n').split("\t")[0][1:-1] # rm "\r\n", split lines by tab ("/t"), remove quotes, and get the ID for the cand. gene
		with open(args.annotation_file, "r") as data:   # open GFF file
		#working through lines 
			for line in data.readlines():
				line = line.rstrip("\r\n")  # remove "\r\n"
				record = line.split("\t")   # split lines by tab 
				ID = re.split("=|;", record[8])[1]   # get the ID
				if ID == gene_name and record[2] == "gene":  # get the line that contains whole gene info
					scaf = record[0]  # scaffold info 
					start = record[3]  # start loc. for that gene
					end = record[4]  # end loc. for that gene 
					length = abs(int(end) - int(start)) # calc. the length of the gene 
					dis = min(abs(int(start) - int(IG_start)), abs(int(start) - int(IG_end)), abs(int(end) - int(IG_end)), abs(int(end) - int(IG_start)))  # calc. distance 
					region = str(scaf+":"+start+"-"+end)  # this is for "region" input for ANGSD. format = scaf:start-end, e.g., DPSCF300010:1000-3000
					if IG_length-1500 <= length <= IG_length+1500 or IG_length*0.5 <= length <= IG_length*2:  # check if length fits the criteria
						length_ok = "Y" 
						K = K + 1   # track num. of usable candidates
					else:
						length_ok = "N" 
					if dis <= 50000: # check if distance fits the criteria
						dis_ok = "Y" 
					else:
						dis_ok = "N" 
					if record[6] == "-":   # record the direction of that gene. if '-', reverse complement. 
						rev = "rc"
					else:
						rev = "f"
					out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_name, scaf, start, end, region, length, dis, rev, length_ok, dis_ok, 0))  #write out 
gene_list.close()
data.close()
out_file.close()

## Track the number of candidates found that meet the size criteria 
# print out the num of cands for each focal gene 
print(str(args.immune_gene) + ": " + str(K))
# Update output lists:  tracking usable immune genes and genes that needs re-run 
out_file_3 = open("updated_focal_gene_list.txt", "a")
if args.run_number == "1": 
	out_file_2 = open("expand_list.txt", "a")  # list for genes that need re-run 
elif args.run_number == "2": 
	out_file_2 = open("expand_list_2.txt", "a")  # list for genes that still fail after re-run (should be empty)
if K >= 8:   # consider 8 as enough candidates
	out_file_3.write(str(args.immune_gene) + "\n")  # update list 
elif K < 8:   # if do not have enough candidates
	os.system("rm " + str(args.immune_gene+"_all_cands_calc.txt"))  # rm the output file created for those genes 
	if 0 < K < 8: 
		out_file_2.write(str(args.immune_gene) + "\n")  # need re-run
out_file_2.close()
out_file_3.close()
