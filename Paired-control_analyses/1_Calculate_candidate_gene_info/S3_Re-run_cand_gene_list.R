#### Find Candidate Control Genes for each Immune Genes #### 

## Modified from S1, for re-running on a small set of focal genes that need to expand to whole chromosome for finding enough candidates.


## Import datasets----
# input a list of monarch immune genes
gene_info <- read.csv("Gene_information_list_0307.csv", header = T)
# input a matrix for all monarch gene locations (modified from the data from neo-Z paper)
loc_list <- read.csv("gene_loc_list.csv", header = T)
rownames(loc_list) <- loc_list[,2]
all_gene_list <- gene_info$New_ID   # extract names as a list
# Input a list of genes needed re-run
focal_list <- read.table("expand_list.txt", header = F)

## Find out all immune gene locations ----
all_genes <- as.character(all_gene_list) # all IG
IG_loc_list <- loc_list[all_genes,] # subset on location data

## Find candidate pairs for each genes ----
# first, remove all immune genes (i.e. pairs cannot be immune genes)
cand_loc_list <- subset(loc_list, !(geneID %in% all_genes))  #rm
dim(cand_loc_list) # check dim 
# focal set (re-runs)
focal_list <- read.table("expand_list.txt", header = F)
focal_list <- as.character(focal_list$V1)
focal_loc_list <- subset(IG_loc_list, (IG_loc_list$geneID %in% focal_list)) # subset
# create a summary of candidate gene numbers 
cand_summary_output <- matrix(0, length(focal_list), 4)
colnames(cand_summary_output) <- c("focal_gene_ID", "num_cands_same_scaff", "use_same_chr_instead", "num_cands_same_chr")
# select cand pair for each, and output
for (N in 1:length(focal_list)){
  cand_summary_output[N, 1] <- focal_list[N]
  # get chromosome info 
  chr <- focal_loc_list[N,"chromosome"]
  # pick those on the same chr
  cand_summary_output[N, 2] <- "N/A"
  cand <- cand_loc_list[which(cand_loc_list$chromosome == chr),]
  cand_summary_output[N, 3] <- "Y"
  cand_summary_output[N, 4] <- nrow(cand)
  write.table(cand, file = paste("cand_genes_for_", focal_loc_list[N,2], ".txt", sep = ""),
              sep = "\t", col.names = T)  # output cand. gene lists
}
# output focal (immune) gene list 
write.table(focal_list, file = "second_focal_gene_list.txt", sep = "\n", row.names = F, col.names = F)
# Output candidate gene summary table 
write.table(cand_summary_output, file = "second_cand_gene_summary.txt", sep = "\t", row.names = F, col.names = T)
