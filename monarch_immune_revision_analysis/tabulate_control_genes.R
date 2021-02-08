control_dir <- "../biorxiv_version1/Paired-control_analyses/2_Select_controls_and_run_permutations/Outputs/Selected_control_gene_lists/"
PairedFiles <- list.files(path = control_dir)

PairedFilesID <- gsub(x = PairedFiles, pattern = "control_gene_list_for_(DPOGS\\d+).txt", replacement = "\\1", perl = T)

readControls <- function(iGene) {
    myCtls <- readLines(paste0(control_dir,PairedFiles[iGene]))  # read file into char vector
    # The control set in each file apparently includes the target immune gene (!?)
    myCtls <- myCtls[ myCtls != PairedFilesID[iGene]] # remove the target gene from control set
    return(myCtls)
}

# Generate a list of control genes for each target immune gene.
controls <- list() 
for (i in seq_along(PairedFilesID)) {
    controls[[PairedFilesID[i]]] <- readControls(i)
}

imctlfile <- "paired_control_geneIDs.txt"
cat("# Immune-control gene pairings\n", file = imctlfile)
cat("ImmuneGene\tControlGenes\n", file = imctlfile, append = T)
for (i in names(controls)) {
	cat(i,"\t", file = imctlfile, append = T)
	outcsv <- paste(controls[[i]], collapse = ",")	
	cat(outcsv,"\n", file = imctlfile, append = T)
}