#### Find Candidate Control Genes for each Immune Genes #### 

### First, find out locations on chromosomes for each immune gene 
# Immune gene list obtained from William Palmer.  
# Location info from Walters lab's neo-Z paper (Mongue etal 2017 - suppl. File S2)
# Exclude genes on the sex chromosomes

## Second, with those information, output all potential candidates per focal immune 
# Exclude known immune genes 
# Output all genes on the same scaffold, and do downstream processing in Python 


## Import datasets----
# input a list of monarch immune genes
gene_info <- read.csv("Gene_information_list_0307.csv", header = T)
all_gene_list <- gene_info$New_ID   # extract names as a list

# input a matrix for all monarch gene locations (modified from the data from neo-Z paper)
loc_list <- read.csv("gene_loc_list.csv", header = T)
rownames(loc_list) <- loc_list[,2]

## Find out all immune gene locations ----
all_genes <- as.character(all_gene_list) # all IG
IG_loc_list <- loc_list[all_genes,] # subset on location data
# list of all immune genes with their locations
IG_loc_list
write.table(IG_loc_list, file = "All_IG_loc_list.txt", sep = ",", col.names = T)  # output
# Some summaries  
length(which(IG_loc_list$ZorA =='A'))  # no. of autosomes
length(which(IG_loc_list$ZorA =='Z_neo'))  # no. on neo-Z
length(which(IG_loc_list$ZorA =='Z_anc'))  # no. on anc-Z
IG_aut <- IG_loc_list[which(IG_loc_list$ZorA =='A'),] # subset only autos
aut_out <- table(IG_aut$chromosome, exclude = NULL)  # dist. of no. on autosomes
write.csv(aut_out, file = "dist_on_autosomes.csv")  # output

## Find candidate pairs for each genes ----
# first, remove all immune genes (i.e. pairs cannot be immune genes)
cand_loc_list <- subset(loc_list, !(geneID %in% all_genes))  #rm
dim(cand_loc_list) # check dim 
# Decide focal set - Here, do all immune genes that are on autosomes 
focal_list <- as.character(IG_aut$geneID)   # gene names 
focal_loc_list <- IG_aut  # location into matrix

# create a summary of candidate gene numbers 
cand_summary_output <- matrix(0, length(focal_list), 4)
colnames(cand_summary_output) <- c("focal_gene_ID", "num_cands_same_scaff", "use_same_chr_instead", "num_cands_same_chr")

# select potential candidates for each, and output
for (N in 1:length(focal_list)){   # for each immune gene
  cand_summary_output[N, 1] <- focal_list[N]   # gene id
  # get scaffold and chromosome info 
  chr <- focal_loc_list[N,"chromosome"]
  scaf <- focal_loc_list[N,"scaffold"]
  # First, try pick those on the same scaffold
  # If no. candidates < 20, do same chr. instead 
  cand <- cand_loc_list[which(cand_loc_list$scaffold == scaf),]   # all genes on same  scaff.
  cand_summary_output[N, 2] <- nrow(cand)   # num genes on same scaff.
  if (nrow(cand) <= 20) {  # if not enough on same scaff. 
    cand_summary_output[N, 3] <- "Y"  # do whole chr. 
    cand <- cand_loc_list[which(cand_loc_list$chromosome == chr),]  # all genes on same chr.
    cand_summary_output[N, 4] <- nrow(cand)  # num genes on same chr.
  }
  else {
    cand_summary_output[N, 3] <- "N"
    cand_summary_output[N, 4] <- "N/A"
  }
  write.table(cand, file = paste("cand_genes_for_", focal_loc_list[N,2], ".txt", sep = ""),
              sep = "\t", col.names = T)  # output cand. gene lists
}
# output focal (immune) gene list 
write.table(focal_list, file = "focal_gene_list.txt", sep = "\n", row.names = F, col.names = F)
# Output candidate gene summary table 
write.table(cand_summary_output, file = "cand_gene_summary.txt", sep = "\t", row.names = F, col.names = T)
