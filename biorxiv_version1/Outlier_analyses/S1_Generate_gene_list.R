#### Finding intersections (across tests and 0f/4f) and generate a gene list for TjD - FST plotting  ####

## read in all stat output files, and find intersection genes of all files. ----
pops <- c("NAM", "FL", "PAC", "ATL")
tests <- c("tjd", "fst_w")
for (i in 1:4){
  pop1 <- pops[i]
  for (j in 1:2){
    Stat <- tests[j]  # get the name of the test 
    ## Input popgen tests data 
    if (!(j == 2 & pop1 == "NAM")){ 
      # input 4-fold data 
      stat_data_4f <- read.table(paste(pop1, "_", Stat, "_output_4f.txt", sep = ""), sep = "\t", header = T) # input data
      stat_data_4f$gene.id <- as.character(paste(stat_data_4f$gene.id))
      if (sum(is.na(stat_data_4f$stat_value)) > 0){
        stat_data_4f <- stat_data_4f[-which(is.na(stat_data_4f$stat_value) == T), ] # remove genes with NaN   
        }
      # input 0-fold data 
      stat_data_0f <- read.table(paste(pop1, "_", Stat, "_output_0f.txt", sep = ""), sep = "\t", header = T) # input data
      stat_data_0f$gene.id <- as.character(paste(stat_data_0f$gene.id))
      if (sum(is.na(stat_data_0f$stat_value)) > 0){
        stat_data_0f <- stat_data_0f[-which(is.na(stat_data_0f$stat_value) == T), ] # remove genes with NaN 
        }
      } else{  # FST & NAM -- not applicable
      next
    }
    D <- intersect(stat_data_0f$gene.id, stat_data_4f$gene.id)
    if (i == 1 & j == 1){
      U <- D
    } 
    else{
      U <- intersect(U, D)  # find intersection 
    }
  }
}   

## Output the updated gene list (intersection of all above for Tjd an Fst)
tjd_fst_list <- as.data.frame(U)
write.table(tjd_fst_list, file = "tjd_fst_gene_list_04x.txt", quote = F, col.names = F, row.names = F)
