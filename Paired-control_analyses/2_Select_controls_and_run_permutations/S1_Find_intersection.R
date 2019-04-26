#### Get intersection of all ANGSD per-gene point estimates to find out NAs #### 
# Some background genes resulted in "NA" in a pre-run using individual gene priors. They also differed across pop/test combinations. 
# Those genes cannot be selected as control. 
# So, need to find the intersection across pops, tests, and sites, to produce a list. 


## Import immune gene informaiton ----
# Input a list of immune genes used for analyses (102 non-Z genes)
gene_list <- read.table("complete_focal_gene_list.txt", header = F) 
str(gene_list)
# input a list of all monarch immune genes with detail information 
gene_info <- read.csv("Gene_information_list_final.csv", header = T)
str(gene_info)
gene_info_used <- subset(gene_info, (New_ID %in% gene_list$V1))  # subset 
dim(gene_info_used)
head(gene_info_used)
rownames(gene_info_used) <- gene_info_used$New_ID  #for later sorting purposes


## Import ANGSD data and find intersections across combinations ----
pops <- c("NAM", "FL", "PAC", "ATL")
tests <- c("tjd", "pi", "wtheta", "fst_w")
for (i in 1:4){
  pop1 <- pops[i]  # get the pop name 
  for (j in 1:4){
    Stat <- tests[j]  # get the test name
    ## Input ANGSD data 
    if (j != 4){  # not FST 
      # input 4-fold data 
      stat_data_4f <- read.table(paste0("pre-run/four_fold/", pop1, "_", Stat, "_output.txt"), sep = "\t", header = T) # input data
      # input 0-fold data 
      stat_data_0f <- read.table(paste0("pre-run/zero_fold/", pop1, "_", Stat, "_output.txt"), sep = "\t", header = T) # input data
      # input full-gene data 
      stat_data_full <- read.table(paste0("pre-run/full_gene/", pop1, "_", Stat, "_output.txt", sep = ""), sep = "\t", header = T) # input full-gene data
    }  
    else if(pop1 != "NAM"){    # if FST, but not NAM 
      POP1 <- paste("NAM.", pop1, sep = "")  # rename
      # input 4-fold data 
      stat_data_4f <- read.table(paste0("pre-run/four_fold/", POP1, "_", Stat, "_output.txt"), sep = "\t", header = T) # input data
      # input 0-fold data 
      stat_data_0f <- read.table(paste0("pre-run/zero_fold/", POP1, "_", Stat, "_output.txt"), sep = "\t", header = T) # input data
      # input full-gene data 
      stat_data_full <- read.table(paste0("pre-run/full_gene/", POP1, "_", Stat, "_output.txt"), sep = "\t", header = T) # input full-gene data
    } 
    else{  # FST & NAM -- not applicable
      next  
    }
    # remove genes with NaN 
    if(sum(is.na(stat_data_4f[, 2] > 0))){
      stat_data_4f <- stat_data_4f[-which(stat_data_4f[, 2] == "NaN"), ]
    }
    if(sum(is.na(stat_data_0f[, 2] > 0))){
      stat_data_0f <- stat_data_0f[-which(stat_data_0f[, 2] == "NaN"), ]
    }
    if(sum(is.na(stat_data_full[, 2] > 0))){
      stat_data_full <- stat_data_full[-which(stat_data_full[, 2] == "NaN"), ]
    }
    # find intersections 
    D <- intersect(stat_data_0f$gene.id, stat_data_4f$gene.id)  # intersection in 0x and 4x data in a given pop-test combination 
    D <- intersect(D, stat_data_full[,1])  # intersection in 0x, 4x, and full-gene data in a given pop-test combination 
    if (i == 1 & j == 1){
      U <- D  # as start 
    }
    else{
      U <- intersect(U, D)  # find intersection and update across loops
      print(U)
    }
  }
}   

## Output ----
length(U)    # this is the list of genes that have no NAs in any combinations  
write.table(U, file = "ANGSD_list_noNA.txt", sep = "\t", col.names = F, row.names = F, quote = F) 
