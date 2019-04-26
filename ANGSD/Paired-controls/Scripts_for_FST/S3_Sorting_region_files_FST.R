#!/usr/bin/env R 

#### Sorting ANGSD region files  ####

file_list <- list.files(pattern = "_group_region.txt")  # list all group region files 

for (j in file_list){
  data <- read.table(j)  # input data 
  data$V1 <- as.character(data$V1)  # convert 
  DATA <- matrix(0, nrow(data), 3)  # create a matrix for storing output 
  # split by pattern 
  for (i in 1: nrow(data)){  
    temp <- strsplit(data[i, 1], "\\:|\\-")
    DATA[i, ] <- c(temp[[1]][1], temp[[1]][2], temp[[1]][3])
  }
  DATA <- as.data.frame(DATA) # convert 
  DATA$V2 <- as.numeric(paste(DATA$V2))
  DATA$V3 <- as.numeric(paste(DATA$V3))
  DATA_sorted <- DATA[order(DATA$V1, DATA$V2),]  # sort first by scaff name (chr) then by start position (num)
  list <- 0
  for (k in 2: nrow(DATA_sorted)){
    if (DATA_sorted[(k-1), 3] > DATA_sorted[k, 2] & DATA_sorted[(k-1), 1] == DATA_sorted[k, 1]){   #  if overlapping 
      # merge the two by replacing the end pos. of the 1st entry with the 2nd entry end pos, and add a NA is the second entry for removal later 
      DATA_sorted[(k-1), 3] <- DATA_sorted[k, 3]   
      DATA_sorted[k, 2] <- NA
    }
  }
  # remove entries with "NA" marked 
  if(length(which(is.na(DATA_sorted$V2))) > 0){  
    DATA_sorted <- DATA_sorted[-which(is.na(DATA_sorted$V2)), ]   # remove rows with NA 
  }
  # output 
  DATA_out <- as.data.frame(paste0(DATA_sorted$V1, ":", DATA_sorted$V2, "-", DATA_sorted$V3))
  write.table(DATA_out, file = j, col.names = F, row.names = F, quote = F, sep = "\t")
  rm(temp)
  rm(DATA)
  rm(DATA_out)
}
