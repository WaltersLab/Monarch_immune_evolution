#!bin/usr/env python 

#### Parse the annotation GFF file (OGS 2.0) to get positions for each gene.  ####
# Output a table of all gene information for further purposes. 

# Import packages 
import argparse, re

# Set command line arguments 
parser = argparse.ArgumentParser(description = "get positions for all genes in the genome")
parser.add_argument("annotation_file", type = str, help = "input the annotation file (.GFF3 file) for OSG2.0")
args = parser.parse_args()

# Create an output file 
# note that the "region" in the output is for generating range files used for ANGSD later
out_file = open("All_gene_positions.txt", "w")  #create an output file
out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("name", "scaffold", "start", "end", "reverse", "region", "length"))  #write out headers  

# Parsing the GFF file to get info for all 15130 genes
with open(args.annotation_file, "r") as data:   # open GFF file
# work through lines to get info for each gene
		for line in data.readlines():
			record = line.rstrip("\r\n").split("\t")  
			if record[2] == "gene": # get the line that contains whole gene info
				ID = re.split("=|;", record[8])[3]   # get the gene ID
				scaf = record[0] # scaffold info
				start = record[3]  # start loc. for that gene
				end = record[4]  # end loc. for that gene 
				length = abs(int(end) - int(start)) + 1  # calc. the length of the gene
				region = str(scaf+":"+start+"-"+end)  # this is for "region" input for ANGSD. format = scaf:start-end, e.g., DPSCF300010:1000-3000
				if record[6] == "-":  # record the direction of that gene. if '-', reverse complement. 
					rev = "rc"
				else:
					rev = "f"
				out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ID, scaf, start, end, rev, region, length))  #write out 
data.close()	
out_file.close()
