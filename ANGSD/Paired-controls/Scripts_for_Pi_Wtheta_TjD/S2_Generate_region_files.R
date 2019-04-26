#### Generating region files (ANGSD input) for 0x/4x site analyses on selected groups of genes  ####
# Generate region files for both groups and individual genes (regions still needed even if running site option).  


## Input data ----
# input a dataframe with position information of all genes 
gene_table <- read.table("All_gene_positions.txt", header = T)   
gene_table$name <- as.character(gene_table$name)
gene_table$region <- as.character(gene_table$region)
# input a dataframe with immune and selected controls for the entire 10,000 runs
control_table <- read.table("CONTROLS.txt", header = T)
# input non-Z immune gene info 
IG_noZ <- read.table("All_Monarch_immune_genes_non-Z.txt", header = T)

## Subsetting ---- 
# Put everything in a function
Generate_region_files <- function(gene_group){
  # decode group name 
  if (gene_group == "Rec"){
    group_name <- "Recognition"
  } else if (gene_group == "Sig"){
    group_name <- "Signalling"
  } else if (gene_group == "Mod"){
    group_name <- "Modulation"
  } else {
    group_name <- "Effector"
  }
  # immune gene group 
  gene_list <- sort(as.character(subset(IG_noZ, Class == group_name)$Gene_ID))  # all genes from the given group
  # conrol gene group 
  gene_con <- subset(control_table, Immune_genes %in% gene_list)[, -1] # all corresponding control genes in all runs
  gene_con <- sapply(gene_con, as.character) # convert all gene names (factors) to characters 
  gene_con_list <- sort(unique(c(gene_con)))  # extract all elements and find unique 
  # output 
  write.table(gene_list, file = paste0(gene_group, "_list_04x.txt"), col = F, row = F, quote = F) # immune gene group list  
  write.table(gene_con_list, file = paste0(gene_group, "_Con_list_04x.txt"), col = F, row = F, quote = F)  # control gene group list 
  
  # Region file for groups 
  gene_regions <- subset(gene_table, name %in% gene_list)$region
  gene_con_regions <- subset(gene_table, name %in% gene_con_list)$region
  #output 
  write.table(gene_regions, file = paste0(gene_group, "_group_region.txt"), col = F, row = F, quote = F)
  write.table(gene_con_regions, file = paste0(gene_group, "_Con_group_region.txt"), col = F, row = F, quote = F)
  # Region files for individual genes 
  for (i in c(gene_list, gene_con_list)){
    write.table(subset(gene_table, name == i)$region, file = paste0(i, "_region.txt"), col = F, row = F, quote = F)
  }
}

# Apply the function 
for (i in c("Rec", "Sig", "Mod", "Eff")){
  Generate_region_files(i)
}
