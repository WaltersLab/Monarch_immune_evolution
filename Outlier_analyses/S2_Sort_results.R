### Generating the distributions of test statistic values across the whole genome  ### 

## Set manual input 
pop = "ATL"   # NAM, FL, PAC, or ATL 
test = "tjd"   # tjd, fst_w
site <- "0f"  # either "0f" or "4f"

## put all in a function 
RUNALL <- function(pop, test, site){

# decode names 
# population 
if (pop == "NAM"){
  population = "North America"
}  else if (pop == "FL"){
    population = "South Florida"
  } else if (pop == "PAC"){
      population = "Pacific (no Hawaii)"
    } else if (pop == "ATL"){
        population = "Atlantic (no Portugal)"
        }
# test name 		
if (test == "tjd"){
  testname = "Tajima's D"
}  else if (test == "fst_w"){
    testname = "FST"
  } else if (test == "wtheta"){
      testname = "Watterson's theta"
    } else if (test == "pi"){
        testname = "Pi"
        }
# site
SITE <- ifelse(site == "0f", "zero_fold", "four_fold") 

# make a function for S.E. 
SE <- function(data){
  se = sd(data)/sqrt(length(data))
  return(se)
}

## Import datasets ---- 
# results -- dataset of stats for all individual genes in the genome 
data <- read.table(paste(pop, "_", test, "_output_", site, ".txt", sep = ""), header = T)  
colnames(data) <- c("gene.id", "stat_value")  # change the col name to make it universal 
# background genes 
gene_list <- read.table("All_aut_ID_list.txt", header = F) # a list of all non-Z genes 
# immune genes 
IG_func <- read.table("All_Monarch_immune_genes_non-Z.txt", header = T)  # a list immune genes with functions, non-Z genes only
# intersected gene list (genes present in all pops and tests)
int_gene_list <- read.table("tjd_fst_gene_list_04x.txt", header = F)   # pre-generated 
#  Take intersections - restrict "genome" and "immune genes" to those present in the intersected gene list 
gene_list <- intersect(int_gene_list$V1, gene_list$V1)
IG_list <- intersect(IG_func$New_ID, int_gene_list$V1)
IG_func <- subset(IG_func, (New_ID %in% IG_list))

## data processing ----
# remmove all Z-link genes 
data_noZ <- subset(data, (gene.id %in% gene_list)) 

# calculate percentiles (for graphic purposes)
data_noZ <- data_noZ[order(data_noZ$stat_value),]  # order by Tajima's D value 
percentile <- (seq(1, nrow(data_noZ), 1) - 1)/(nrow(data_noZ) - 1)  # percentile = (index - 1) / (N-1)
data_noZ <- cbind(data_noZ, percentile)  # combine
data_noZ <- data_noZ[order(data_noZ$gene.id), ]  # convert back to the orignial order (by name)

# subset of all non-Z immune genes
IG_noZ <- subset(data_noZ, (gene.id %in% IG_list)) 
# subset of all non-Z immune genes in Recognition 
Rec <- IG_func[which(IG_func$Class == "Recognition"),]
Rec_list <- as.character(Rec$New_ID)
Rec_noZ <- subset(data_noZ, (gene.id %in% Rec_list)) 
# subset of all non-Z immune genes in Signalling 
Sig <- IG_func[which(IG_func$Class == "Signalling"),]
Sig_list <- as.character(Sig$New_ID)
Sig_noZ <- subset(data_noZ, (gene.id %in% Sig_list)) 
# subset of all non-Z immune genes in Modulation
Mod <- IG_func[which(IG_func$Class == "Modulation"),]
Mod_list <- as.character(Mod$New_ID)
Mod_noZ <- subset(data_noZ, (gene.id %in% Mod_list)) 
# subset of all non-Z immune genes in Effector
Eff <- IG_func[which(IG_func$Class == "Effector"),]
Eff_list <- as.character(Eff$New_ID)
Eff_noZ <- subset(data_noZ, (gene.id %in% Eff_list)) 


## Report outliers ---- 
up_lim = 0.975
low_lim = 0.025
global_mean <- mean(data_noZ$stat_value)  # mean for all genes
# Find out "outliers" based on the quantile used 
upper <- IG_noZ$gene.id[which(IG_noZ$percentile >= up_lim)]  # find those > 95% 
lower <- IG_noZ$gene.id[which(IG_noZ$percentile <= low_lim)] # find those < 5%  
upper_genes <- subset(IG_func, (New_ID %in% upper))  #get names 
lower_genes <- subset(IG_func, (New_ID %in% lower)) # get names 
upper_genes2 <- cbind(c(rep("Upper", nrow(upper_genes))), upper_genes)
colnames(upper_genes2) <- c("Side", colnames(upper_genes))
lower_genes2 <- cbind(c(rep("Lower", nrow(lower_genes))), lower_genes)
colnames(lower_genes2) <- c("Side", colnames(lower_genes))
OUT <- rbind(upper_genes2, lower_genes2)

IG_noZ$gene.id <- droplevels(IG_noZ$gene.id)  # drop redundant levels
IG_func$New_ID <- droplevels(IG_func$New_ID)  # drop redundant levels
rownames(IG_noZ) <- IG_noZ$gene.id  # add rownames for sorting purpose 
IG_noZ <- IG_noZ[IG_func$New_ID,]  # sort to the same as the function list
T_std <- (IG_noZ$stat_value - global_mean)/sd(data_noZ$stat_value)  # standardized by background dist. 
IG_noZ_out <- cbind(IG_func, IG_noZ[,2:3], T_std, matrix(0, nrow(IG_func), 1))# combine
outliers <- c(as.character(upper_genes$New_ID), as.character(lower_genes$New_ID)) # those "outliers"
for (i in 1: nrow(IG_noZ_out)){
  ifelse(IG_noZ_out[i, 1] %in% outliers, IG_noZ_out[i, 8] <- 1 , IG_noZ_out[i, 8] <- 0)  # note outliers 
}
colnames(IG_noZ_out) <- c(head(colnames(IG_noZ_out), n = length(colnames(IG_noZ_out)) -1), "outliers") 
write.table(IG_noZ_out, file = paste(pop, "_", test, "_", site, "_all_IG_full_results.txt", sep = ""), sep = "\t")
}


### Run all 
for (i in c("NAM", "FL", "PAC", "ATL")){
  for (j in c("tjd", "fst_w")){
    for (k in c("0f", "4f")){
      if (!(i == "NAM" & j == "fst_w")){
        RUNALL(i, j, k)
      }
    }
  }
}
