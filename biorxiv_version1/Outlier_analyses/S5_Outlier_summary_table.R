#### Compiling a list of outlier genes based on 0-fold results ####

## Get total list of outlier genes for Tjd and FST seperately ----
# Input outlier list per class 
setwd("C:/R/Immune_genes/angsd_new_results/GB_5percent")
rec <- read.table("Recognition_tjd_fst_w_0f_outliers.txt", header = T)
sig <- read.table("Signalling_tjd_fst_w_0f_outliers.txt", header = T)
mod <- read.table("Modulation_tjd_fst_w_0f_outliers.txt", header = T)
eff <- read.table("Effector_tjd_fst_w_0f_outliers.txt", header = T)
data <- rbind(rec, sig, mod, eff)  # combine 
dim(data)
str(data)
# Get all Tjd outlier names 
data_tjd <- data[, -(5:6)]
data_tjd <- data_tjd[which(data_tjd$T1_outliers == 1), ]
dim(data_tjd)
tjd_outliers <- unique(data_tjd$New_ID)
tjd_outliers
# Get all FST outlier names 
data_fst <- data[, -(3:4)]
data_fst <- data_fst[which(data_fst$T2_outliers == 1), ]
dim(data_fst)
fst_outliers <- unique(data_fst$New_ID)
fst_outliers
# Intersect
both_outliers <- intersect(tjd_outliers, fst_outliers)
# output those lists
write.table(tjd_outliers, file = "tjd_outlier_list.txt", col.names = F, row.names = F)
write.table(fst_outliers, file = "fst_outlier_list.txt", col.names = F, row.names = F)
write.table(both_outliers, file = "both_outlier_list.txt", col.names = F, row.names = F)

## Extracting information needed for those outliers and make output ---- 
# Import full results ---- 
pops <- c("NAM", "FL", "PAC", "ATL") # set pops
tests <- c("tjd", "fst_w") # set tests 
sites <- c("0f", "4f") # set sites 
# import all results 
for (k in 1:2){
  site <- sites[k]
  for (i in 1:2){  # for all 2 tests 
    test_in <- tests[i]
    for (j in 1:4){  # for all 4 pops 
      pop = pops[j]
      if (!((test_in == "fst_w") & (pop == "NAM"))){  # exclude the NAM-FST that didn't exist 
        data <- read.table(paste(pop, "_", test_in, "_", site, "_all_IG_full_results.txt", sep = ""), header = T)
        # assign data
        assign(paste(pop, "_", test_in, "_", site, "_data", sep = ""), data) 
        rm(data)
      }	
    }
  }
}

# Compile Tajima's D table ---- 
# First use the NAM 0f as template for the output 
TJD_DATA <- subset(NAM_tjd_0f_data, New_ID %in% tjd_outliers)  #subset 
TJD_DATA <- TJD_DATA[, c(4, 1, 2, 3, 5, 6, 8)] # get things need
TJD_DATA <- cbind(data.frame(0), TJD_DATA) # add a col for nums (sorting purpose)
TJD_DATA[,1] <- rownames(TJD_DATA)
# Add NAM 4f 
ADD <- subset(NAM_tjd_4f_data, New_ID %in% tjd_outliers)
ADD <- ADD[, c(1, 5, 6, 8)]
head(ADD)
TJD_DATA <- merge(TJD_DATA, ADD, by.x = "New_ID", by.y = "New_ID")
colnames(TJD_DATA) <- c("New_ID", "num", "Gene_name", "Class", "Pathway", "NAM_0f", "NAM_0f_p", "NAM_0f_out", "NAM_4f", "NAM_4f_p", "NAM_4f_out")
rm(ADD)

# Add the other 3 pops through loop 
for (i in c("FL", "PAC", "ATL")){
  ADD_0f <- subset(get(paste0(i, "_tjd_0f_data")), New_ID %in% tjd_outliers)
  ADD_0f <- ADD_0f[, c(1, 5, 6, 8)]
  colnames(ADD_0f) <- c("New_ID", paste0(i, "_0f"), paste0(i, "_0f_p"), paste0(i, "_0f_out"))
  ADD_4f <- subset(get(paste0(i, "_tjd_4f_data")), New_ID %in% tjd_outliers)
  ADD_4f <- ADD_4f[, c(1, 5, 6, 8)]
  colnames(ADD_4f) <- c("New_ID", paste0(i, "_4f"), paste0(i, "_4f_p"), paste0(i, "_4f_out"))
  TJD_DATA <- merge(TJD_DATA, ADD_0f, by.x = "New_ID", by.y = "New_ID")
  TJD_DATA <- merge(TJD_DATA, ADD_4f, by.x = "New_ID", by.y = "New_ID")
  rm(ADD_0f)
  rm(ADD_4f)
}
head(TJD_DATA)
TJD_DATA$num <- as.numeric(TJD_DATA$num)  # convert to number to sort 
TJD_DATA <- TJD_DATA[order(TJD_DATA$num), ]  # sort back to orinigal order 
head(TJD_DATA)
TJD_DATA <- TJD_DATA[, c(3, 1, 4:ncol(TJD_DATA))]  # rearrange cols 
head(TJD_DATA)
# output 
write.table(TJD_DATA, file = "Tjd_all_outlier_info_table.txt", row.names = F, col.names = T, quote = F, sep = "\t")

# Compile FST table ---- 
# First use the FL 0f as template for the output 
FST_DATA <- subset(FL_fst_w_0f_data, New_ID %in% fst_outliers)  #subset 
FST_DATA <- FST_DATA[, c(4, 1, 2, 3, 5, 6, 8)] # get things need
FST_DATA <- cbind(data.frame(0), FST_DATA) # add a col for nums (sorting purpose)
FST_DATA[,1] <- rownames(FST_DATA)
# Add FL 4f 
ADD <- subset(FL_fst_w_4f_data, New_ID %in% fst_outliers)
ADD <- ADD[, c(1, 5, 6, 8)]
head(ADD)
FST_DATA <- merge(FST_DATA, ADD, by.x = "New_ID", by.y = "New_ID")
colnames(FST_DATA) <- c("New_ID", "num", "Gene_name", "Class", "Pathway", "FL_0f", "FL_0f_p", "FL_0f_out", "FL_4f", "FL_4f_p", "FL_4f_out")

# Add the other 2 pops through loop 
for (i in c("PAC", "ATL")){
  ADD_0f <- subset(get(paste0(i, "_fst_w_0f_data")), New_ID %in% fst_outliers)
  ADD_0f <- ADD_0f[, c(1, 5, 6, 8)]
  colnames(ADD_0f) <- c("New_ID", paste0(i, "_0f"), paste0(i, "_0f_p"), paste0(i, "_0f_out"))
  ADD_4f <- subset(get(paste0(i, "_fst_w_4f_data")), New_ID %in% fst_outliers)
  ADD_4f <- ADD_4f[, c(1, 5, 6, 8)]
  colnames(ADD_4f) <- c("New_ID", paste0(i, "_4f"), paste0(i, "_4f_p"), paste0(i, "_4f_out"))
  FST_DATA <- merge(FST_DATA, ADD_0f, by.x = "New_ID", by.y = "New_ID")
  FST_DATA <- merge(FST_DATA, ADD_4f, by.x = "New_ID", by.y = "New_ID")
  rm(ADD_0f)
  rm(ADD_4f)
}
head(FST_DATA)
FST_DATA$num <- as.numeric(FST_DATA$num)  # convert to number to sort 
FST_DATA <- FST_DATA[order(FST_DATA$num), ]  # sort back to orinigal order 
head(FST_DATA)
FST_DATA <- FST_DATA[, c(3, 1, 4:ncol(FST_DATA))]  # rearrange cols 
head(FST_DATA)
# output 
write.table(FST_DATA, file = "Fst_all_outlier_info_table.txt", row.names = F, col.names = T, quote = F, sep = "\t")
