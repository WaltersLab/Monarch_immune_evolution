#### Identify control genes for each immune gene and select for each permutation ####
# Select control genes with permuting between immune and control identity


## Import datasets ---- 
# imune a list of genes that have ANGSD point estimates
I_list <- read.table("ANGSD_list_noNA.txt")
I_list <- as.character(I_list$V1)  # convert to chr 
# list of immune genes (102 non-Z genes)
gene_list <- read.table("complete_focal_gene_list.txt", header = F) 
str(gene_list)

## Generate the "control" set for permutations ----
# output the control set for each immune gene 
# output all information needed for permuation tests, and use it consistently for all popgen stats. 
# generate a large matrix, with rows as corresponding immune genes, and cols as permutation runs (test statistic + N permutations)

# set up 
set.seed(1011)  # Set random seed 
N = 10000   # number of permutations 
PERM = matrix(0, nrow(gene_list), N+1)  # create a matrix for storing selected controls 
colnames(PERM) <- c("Test_statistic", paste0("run_", seq(1:N)))
rownames(PERM) <- gene_list$V1
ALL_CONS <- list()  # a list of lists of control genes per immune gene

# loop through all immune genes 
for (i in 1:nrow(gene_list)){
  gene <- as.character(gene_list[i,])  # get focal immune gene name 
  
  # Input candidate gene pool and processing 
  candidate_file <- read.table(paste("all_cands_info/", gene, "_all_cands_calc.txt", sep = ""), header = T)  # input the corresponding candidate table 
  candidate_list <- candidate_file[-which(candidate_file$length_ok == "N"), 1] # only keep candidates that fit the size criteria (i.e., focal & potential controls). Only keep gene IDs.
  candidate_list <- subset(candidate_list, (candidate_list %in% I_list)) # only keep candidate genes that is in the intersection of all output files (i.e., no missing data)
 
  # Double check if any of the genes have too few controls (fewer than 6). 
  if((length(candidate_list) -1) < 6){  # if less than 6 controls 
    print(paste(gene, "(", gene_info[which(gene_info$New_ID == gene), 3], ")",
                "num:", length(candidate_list)))  # print a message 
  }
  
  # Output the control list (note that the first gene is the focal immune gene)
  write.table(candidate_list, file = paste0("control_gene_lists/control_gene_list_for_", gene, ".txt"), sep = "\t", quote = F, col.names = F, row.names  = F)
  ALL_CONS[length(ALL_CONS)+1] <- list(as.character(candidate_list)) # append to the list item
  
  # random selection for each permutation rounds (immune gene can be selected. select with replacement)
  set <- sample(candidate_list, N, replace = T)  # selected set for this focal immune gene 
  PERM[i, ] <- c(gene, as.character(set)) # store it in the PERM matrix 
}


## Output the control list (of lists) and the permutation matrix ---- 
# output as R data for next step 
save(ALL_CONS, file = "ALL_CONTROLS.RData")
save(PERM, file = "PERM_MATRIX.RData")

