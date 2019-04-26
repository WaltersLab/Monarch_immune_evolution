#### Contingency table for Tajima's D - FST outliers  ####

## Input datasets -- these are the outputs from the 2D plot generating script 
classes <- c("Recognition", "Signalling", "Modulation", "Effector")  # set gene classes 
short_name <- c("Rec", "Sig", "Mod", "Eff")
Ns <- c(19, 41, 28, 14)    # the number of total genes respectively in each gene classes. These numbers obtained from hte 2D plot script 
input <- matrix(NA, 8, 4)  # create an input matrix 
input[,1] <- rep(classes, 2)
input[,2] <- rep(short_name, 2)
input[,3] <- c(rep("0f", 4), rep("4f", 4))
input[,4] <- rep(Ns, 2)
colnames(input) <- c("Class","Name","Site", "N")
input <- as.data.frame(input)
input$N <- as.numeric(as.character(input$N))
str(input)

# input via loop 
for (i in 1:nrow(input)){
  # read in the corresponding data and assign 
  data <- read.table(paste(input$Class[i], "_tjd_fst_w_", input$Site[i], "_outliers.txt", sep = ""), 
             header = T)
  # converting coordinates info into 4 quadrates  
  mat <- cbind(data[,1:2], data.frame(0, 0, 0, 0))  #for storing results 
  colnames(mat) <- c("pop", "New_ID", "Q1", "Q2", "Q3", "Q4")
  mat$Q1[which(data$T1 >= 0.5 & data$T2 >= 0.5)] <- 1         # Note that NAMs turned only into either Q1 or Q2
  mat$Q2[which(data$T1 < 0.5 & data$T2 >= 0.5)] <- 1 
  mat$Q3[which(data$T1 < 0.5 & data$T2 < 0.5)] <- 1 
  mat$Q4[which(data$T1 >= 0.5 & data$T2 < 0.5)] <- 1 
  # assign 
  assign(paste(input$Name[i], "_", input$Site[i], sep = ""), mat) 
  # print unique gene names per file 
  print(paste(input$Name[i], "_", input$Site[i], sep = ""))
  print(length(unique(mat$New_ID)))
  # rm 
  rm(data)
  rm(mat)
}

## Contingency Table: table per functional class (pops as rows) ----
# create a matrix for storing all information and for output 
ALL_output <- matrix(NA, 24, 9)
colnames(ALL_output) <- c("Class", "Site", "Pop", "Q1", "Q2","Q3", "Q4", "NS", "Pval")

# create a matrix for stats results 
stats_output <- cbind(input, data.frame(NA))
colnames(stats_output) <- c(colnames(input), "Pval")

for (i in 1:nrow(input)){
  # create a contingency table 
  OUTPUT <- matrix(NA, 3, 5)  
  colnames(OUTPUT) <- c("Q1", "Q2", "Q3", "Q4", "NS")
  rownames(OUTPUT) <- c("FL", "PAC", "ATL")
  OUTPUT[,5] <- input$N[i]   # assign the total num of genes per group   
  # calculate counts in each cell
  for (j in rownames(OUTPUT)){
    sum <- subset(get(paste(input$Name[i], "_", input$Site[i], sep = "")), pop == j)  # subset per pop
    OUTPUT[j, 1:4] <- colSums(sum[3:6])  # use colsums to add up the counts, and fill in the table 
    OUTPUT[j,5] <-   OUTPUT[j,5] - sum(colSums(sum[3:6]))   # update the NS (all - sum(Q1~Q4))
  }
  # do Fisher's exact test  (seems like it deal with zero-sum cols /rows automatically)
  stats_output$Pval[i] <- fisher.test(OUTPUT)$p.value
  # fill out the output matrix 
  start <- i*3 -2
  end <- i*3
  ALL_output[start:end, 1] <- as.character(input$Name[i])
  ALL_output[start:end, 2] <- as.character(input$Site[i])
  ALL_output[start:end, 3] <- as.character(rownames(OUTPUT))
  ALL_output[start:end, 4:8] <- OUTPUT
  ALL_output[start, 9] <- fisher.test(OUTPUT)$p.value
}
# output 
write.table(ALL_output, file = "Contingency_table_per_class.txt", col.names = T, row.names = F, sep = "\t", quote = F)


## Contingency Table: table per populations (classes as rows) ----
# create a matrix for storing all information and for output 
ALL_output2 <- matrix(NA, 32, 9)
colnames(ALL_output2) <- c("Pop", "Site","Class", "Q1", "Q2","Q3", "Q4", "NS", "Pval")

# create a matrix for stats results 
stats_output2 <- cbind(input, data.frame(NA))[, -1]
pops <- c("NAM", "FL", "PAC", "ATL")
stats_output2[,1] <- rep(pops, 2) 
colnames(stats_output2) <- c("Pop", "Site", "N", "Pval")

for (i in 1:nrow(stats_output2)){
  # create a contingency table 
  OUTPUT <- matrix(NA, 4, 5)  
  colnames(OUTPUT) <- c("Q1", "Q2", "Q3", "Q4", "NS")
  rownames(OUTPUT) <- c(as.character(input$Name[1:4]))
  OUTPUT[,5] <- input$N[1:4]   # assign the total num of genes per group 
  # calculate counts in each cell
  for (j in input$Name[1:4]){
    sum <- subset(get(paste(j, "_", stats_output2$Site[i], sep = "")), pop == stats_output2$Pop[i])  # subset per pop
    OUTPUT[j, 1:4] <- colSums(sum[3:6])  # use colsums to add up the counts, and fill in the table 
    OUTPUT[j,5] <-   OUTPUT[j,5] - sum(colSums(sum[3:6]))   # update the NS (all - sum(Q1~Q4))
    }
  # do Fisher's exact test  (seems like it deal with zero-sum cols /rows automatically)
  if (stats_output2$Pop[i] == "NAM"){  # remove Q3 and Q4 
    OUTPUT2 <- OUTPUT[, -c(3:4)]
    stats_output2$Pval[i] <- fisher.test(OUTPUT2)$p.value
  }
  stats_output2$Pval[i] <- fisher.test(OUTPUT)$p.value
  # fill out the output matrix 
  start <- i*4 -3 
  end <- i*4
  ALL_output2[start:end, 1] <- stats_output2$Pop[i]
  ALL_output2[start:end, 2] <- as.character(stats_output2$Site[i])
  ALL_output2[start:end, 3] <- as.character(input$Name[1:4])
  ALL_output2[start:end, 4:8] <- OUTPUT
  ALL_output2[start, 9] <- fisher.test(OUTPUT)$p.value
}
# output 
write.table(ALL_output2, file = "Contingency_table_per_pop.txt", col.names = T, row.names = F, sep = "\t", quote = F)

