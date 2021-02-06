#### Perform permutation tests for each popgen stats  #### 
# performing the paired-control analyses using files generated previously 
# See "method.pdf" scanned doc. for method details 

## Load packages
library(hash)

## Import datasets ---- 
# list of lists of control genes 
load("ALL_CONTROLS.RData")  
# all permutations 
load("PERM_MATRIX.RData")  
# Input a list of immune genes used for analyses (102 non-Z genes)
gene_list <- read.table("complete_focal_gene_list.txt", header = F) 
# input a list of all monarch immune genes with detail information 
gene_info <- read.csv("Gene_information_list_final.csv", header = T)
gene_info_used <- subset(gene_info, (New_ID %in% gene_list$V1))  # subset 
rownames(gene_info_used) <- gene_info_used$New_ID  #for later sorting purposes


### Paired control permutation calculations and output----
N <- ncol(PERM) - 1   # number of permutations 

# put all in a function 
Permutation_Test <- function(INPUT){  # input population and site combination 
pop <- INPUT[1]  # pop  = "NAM", "FL", "PAC", "ATL"
site <- INPUT[2] # site = "0f, "4f", "full 
  
## Create a matrix for the output summary table ----
OUTPUT <- matrix(NA, 12, 5)  # empty matrix
colnames(OUTPUT) <- c("All immune", "Recognition", "Signalling", "Modulation", "Effector")
rownames(OUTPUT) <- c("Tajima's D: test statistic", "Tajima's D: > 0 (%)", "Tajima's D: P-value", 
                      "Pi: test statistic", "Pi: > 0 (%)", "Pi: P-value.", 
                      "Watterson's theta: test statistic", "Watterson's theta: > 0 (%)", "Watterson's theta: P-value", 
                      "FST: test statistic", "FST: > 0 (%)", "FST: P-value"
)


### Loop through each popgen stats ----
Stats <- c("tjd", "pi", "wtheta", "fst_w")  # all four stats 
for (s in 1:4){ 
  Stat <- Stats[s]  # get the name of the test 
  ## Input corresponding ANGSD data  
  if (!(pop == "NAM" & s == 4)){   # not the FST & NAM combination
    stat_data <- read.table(paste0("ANGSD_data/", pop, "_", Stat, "_output_", site, ".txt"), sep = "\t", header = T) # input data
  # Just double checking: remove genes with either "NaN" or "NA"
    if(length(which(is.na(stat_data$stat_value))) > 0){  
      stat_data <- stat_data[-which(is.na(stat_data$stat_value)), ]
    }
  } else {  # FST & NAM -- not applicable
    next  # skip 
  }
  
  ## build a dictionary for gene id: stats value
  stats_dict <- hash(as.character(stat_data$gene.id), stat_data$stat_value) # build a dict.
  
  ## Convert gene names in PERM matrix to to stat values 
  data_mat <- matrix(0, nrow(gene_list), N+1)  # a matrix of stats values for focal & control genes for all permutations
  colnames(data_mat) <- colnames(PERM)
  rownames(data_mat) <- rownames(PERM)
  # convert by the dictionary 
  for (i in 1:nrow(PERM)){
    for (j in 1:(N+1)){
      key <- as.character(PERM[i, j])  # the gene name 
      data_mat[i, j] <- stats_dict[[key]]  # store converted stat value to data mat
    }
  }

  ## Calculating test statistics and permuted distribution 
  # the "A vector"   
  NS <- unlist(lapply(ALL_CONS, length))    # get length for controls and flatten 
  A_vec <- NS/(NS-1)   # see "method.pdf" for detials 
  # the "B vector"
  B_vec <- numeric(nrow(data_mat))  # empty vector
  for (i in 1:nrow(data_mat)){
    Cons = ALL_CONS[[i]]   # get all the control genes for a given immune gene 
    Con_sub <- subset(stat_data, gene.id %in% Cons) # subset the ANGSD data 
    B_value <- sum(Con_sub$stat_value)/(nrow(Con_sub)-1)  # see "method.pdf" for detials
    B_vec[i] <- B_value  # store 
  }
  ## Calculation  (see "method.pdf" for detials)
  Sum_mat <- A_vec*data_mat - B_vec    # this is the raw data of test statistics& permutation
  colnames(Sum_mat) <- colnames(PERM)
  rownames(Sum_mat) <- rownames(PERM)
  
  ## Summary stats per functional group 
  plot_data <- matrix(0, N+1, 5)  # a matrix for storing raw values for plotting purpose 
  class <- c("All immune", "Recognition", "Signalling", "Modulation", "Effector")
  colnames(plot_data) <- class
  rownames(plot_data) <- c("Immune_genes", paste0("run_", seq(1:N)))
  # loop through each class 
  for (k in 1:5){
    group <- class[k]
    if (k == 1){   # no subsetting - for all immune genes 
      Sub <- Sum_mat  # use the full matrix
      Sub_plot <- data_mat
    } 
    else{   # for functional classes 
      sub_names <- gene_info_used[which(gene_info_used$Class == group), ]  # get gene names of the class
      Sub <- Sum_mat[rownames(sub_names), ]  # subset the matrix
      Sub_plot <- data_mat[rownames(sub_names), ]
    }
    # Calculating the p-value (see "method.pdf" for detials )
    M <- colSums(Sub)
    Prop <- sum(M[2:(N+1)] < M[1])/N   # the proportion that the test statistic is greater 
    Pvalue <- sum(abs(M[2:(N+1)]) > abs(M[1]))/N   # the p -value 
    
    # Write to the output summary table 
    OUTPUT[(1 + (s-1)*3), k] <- round(M[1], digits = 4)   # test statistic
    OUTPUT[(2 + (s-1)*3), k] <- round(Prop*100, digits = 4)  #  greater (%)
    OUTPUT[(3 + (s-1)*3), k] <- round(Pvalue, digits = 4)  # P-value  
    
    # Histogram 
    png(paste0("Histograms/Hist_for_", pop, "_", site, "_", Stat, "_", group, ".png"))
    Mh <- colMeans(Sub)
    hist(Mh[2:(N+1)], xlab = Stat, main = paste(pop, site, Stat, group, sep = "  ")) # histogram 
    abline(v = Mh[1], col = "red") # test statistic 
    abline(v = mean(Mh[2:(N+1)]), col = "blue", lty = 2) # mean 
    legend("topright", c("Test statistic", "Mean"), col = c("red", "blue"), lty = c(1, 2))
    dev.off()
    
    # output for plotting 
    plot_data[, k] <- colMeans(Sub_plot)  # calculate the mean across genes and fill in 
    write.table(plot_data, file = paste0("data_for_plotting/Paired_control_rawdata_for_",pop, "_", Stat, "_", site, ".txt"), sep = "\t", col.names = T, row.names = T) # output 
  }
}

## Output the summary table 
write.table(OUTPUT, file = paste0("Paired_control_full_results_for_", pop, "_", site, ".txt"), sep = "\t", col.names = T, row.names = T) # output

return(OUTPUT)
}


### Run the function ---- 
# Make all combinations 
Combs <- matrix(0, 12, 2)
Combs[, 1] <- rep(c("NAM", "FL", "PAC", "ATL"), each = 3)  # pops 
Combs[,2] <- rep(c("0f", "4f", "full"), 4)

# apply to the function 
apply(Combs, 1, Permutation_Test)
