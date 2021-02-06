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

# Input a list of all non-Z genes in the genome 
gene_list_noZ <- read.table("All_aut_ID_list.txt") 

# Subset the data to exclude Z-linked genes 
data2 <- subset(data, transcript %in% gene_list_noZ$V1)
dim(data2)

# extract all 0x positions 
zero_f <- data2[which(data2$degeneracy == "0"), ]
dim(zero_f)
zero_f_out <- paste(zero_f$scaffold, zero_f$coord.genome, sep = "\t")  # convert to the right format for ANGSD
head(zero_f_out)

# extract all 4x positions
four_f <- data2[which(data2$degeneracy == "4"), ]
dim(four_f)
four_f_out <- paste(four_f$scaffold, four_f$coord.genome, sep = "\t")  # convert to the right format for ANGSD
head(four_f_out)

# Output AGNSD files (for whole genome)
write.table(zero_f_out, file = "All_0f_sites_nonZ.txt", row.names = F, col.names = F, sep = "\t", quote = F) 
write.table(four_f_out, file = "All_4f_sites_nonZ.txt", row.names = F, col.names = F, sep = "\t", quote = F) 


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
colnames(count_table) <- c("gene_id", "n.0f_sites", "n.4f_sites")  # change col name
write.table(count_table, file = "All_gene_nSites.txt", row.names = F, col.names = T, quote = F, sep = "\t") # output


## Generate ANGSD input files ---- 
# Put everything in a function
Generate_sites_files <- function(site){  #site = 0f or 4f 
  SITE <- ifelse(site == "0f", "zero_f", "four_f")  # decode name 
  gene_list <- unique(get(SITE)$transcript)  # list of genes 
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
Generate_sites_files("4f")
Generate_sites_files("0f")
