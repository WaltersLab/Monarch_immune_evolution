#### Extract 0-fold and 4-fold positions for all genes #### 
# This is a script that takes a data matrix (from Walter's lab previous work) and extract 0x and 4x positions. 
# Creates input (site) files per gene for all genes in the genome, for running ANGSD.  

## Input data ---- 
load("Dplex_0x4x.Rdata") 
data <- my.coords.table
data <- as.data.frame(data)
str(data)
# modify transcript name to gene ID 
data$transcript <- as.character(gsub("-TA", "", data$transcript))
head(data)
# extract all 0x positions 
zero_f <- data[which(data$degeneracy == "0"), ]
dim(zero_f)
# extract all 4x positions
four_f <- data[which(data$degeneracy == "4"), ]
dim(four_f)

## Generate a gene info output table ---- 
# This shows how many 0x or 4x sites each gene has. 
# which can be merged with the one obtained from regions. 
library(plyr)
count_table_0f <- count(zero_f, 'transcript')  # freq. of gene id counts, which equals number of sites for that gene. 
dim(count_table_0f)
count_table_4f <- count(four_f, 'transcript')
dim(count_table_4f)
count_table <- merge(count_table_0f, count_table_4f, by.x = "transcript", by.y = "transcript") # merge the 0x and 4x tables. 
dim(count_table)
head(count_table)
colnames(count_table) <- c("gene_id", "n.0f_sites", "n.4f_sites")  # change col name
setwd("C:/R/Immune_genes/generate_inputs")
write.table(count_table, file = "All_gene_nSites.txt", row.names = F, col.names = T, quote = F, sep = "\t") # output


## Generate ANGSD input files ---- 
Generate_sites_files <- function(gene_group, site){   # gene_group = group name; site = 0f or 4f 
  SITE <- ifelse(site == "0f", "zero_f", "four_f")  # decode name 
  # Read in pre-generated lists of gene groups  
  gene_list <- read.table(paste0(gene_group, "_list_04x.txt"), header = F)
  gene_list <- as.character(gene_list$V1) # convert type 
  
  # generate angsd site file: for the entire group 
  group_site_list <- subset(get(SITE), transcript %in% gene_list)  # subsetting 
  if (nrow(group_site_list) == 0){  # in case there is no sites
    print(paste("An error happened, not sites found: ", gene_group, site))
  }  
  group_site_file <- paste(group_site_list$scaffold, group_site_list$coord.genome, sep = "\t")  # convert to the right format for ANGSD
  write.table(group_site_file, file = paste0(gene_group, "_group_sites_", site, ".txt"), row.names = F, col.names = F, sep = "\t", quote = F)  # output
  
  # generate angsd site file: individual genes 
  for (i in gene_list){
    site_list <- subset(get(SITE), transcript == i)  # find match 
    site_file <- paste(site_list$scaffold, site_list$coord.genome, sep = "\t")  # convert to the right format for ANGSD
    if (nrow(site_list) == 0){  # in case there is no sites
      print(paste("An error happened, not sites found: ", gene_group, site))
    }
    write.table(site_file, file = paste0(i, "_sites_", site, ".txt"), row.names = F, col.names = F, sep = "\t", quote = F)  # output
  }
}

## Apply the function 
combs <- matrix(0, 16, 2)  # create an empty matrix to store the combinations 
combs[1:8, 1] <- rep(c("Rec", "Sig", "Mod", "Eff"), each = 2) # immune groups 
combs[9:16, 1] <- paste0(combs[1:8, 1], "_Con") # control groups
combs[, 2] <- c("0f", "4f") # site
# loop 
for (i in 1:nrow(combs)){
  Generate_sites_files(combs[i, 1], combs[i, 2])
}
