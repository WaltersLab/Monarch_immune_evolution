#### Multi-population plots #### 
# Fig. 6 and Fig. S2
# Plotted using absolute values rather than differences between immune and controls 

## Load packages and set dir ---- 
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)
library(ggpubr)
setwd("")
data.dir <- ""

## Import datasets and processing ---- 
# create a list of filenames to read in 
pairctl.files <- list.files(path = data.dir, pattern = "Paired")   # get a list of all relavent data 
pop.stat.labels <- gsub(pattern = "Paired_control_rawdata_for_(\\S+).txt", rep = "\\1", x = pairctl.files, perl = T)  # create name labels 
names(pairctl.files) <- pop.stat.labels  # add to the list 
# input the datasets 
pairctl.data <- lapply(pairctl.files, function(x) {
  read.delim(paste0(data.dir, x), as.is = T)
})
# sort them and remove unneeded ones 
myStats <- c("pi_0f", "pi_4f", "pi_full", "wtheta_0f","wtheta_4f", "wtheta_full","tjd_0f", "tjd_4f","tjd_full")  # a list of combinations 
stats.list <- list()  # create an empty list 
for (i in myStats) {   # loop through the list of combs
  stats.list[[i]] <- pairctl.data[grep(i, x = names(pairctl.data) )]   # extract dataframes that belongs to the given group
  names(stats.list[[i]]) <- gsub(pattern = "([A-Z]+)_.+", rep = "\\1", x = names(stats.list[[i]]), perl = T)  # add the name as the pop name 
}
str(stats.list) 

# Seperate obs and controls, and melt each of them. 
stats.list.obs <- list() # create empty lists 
stats.obs.melt <- list()
stats.list.ctl <- list()
stats.ctl.melt <- list()

for (i in myStats) {  #loop through all combs 
  ## Obs: the immune data 
  stats.list.obs[[i]] <- lapply(stats.list[[i]], function(x) {return(x[3,])})  # extract the "immune" data of the group, across all pops 
  # set the means, uppers, lowers, and then merged together 
  means <- stats.list.obs[[i]]   # means 
  means.melt <- suppressMessages(melt(means))  # melt it
  names(means.melt) <- c("class", "statistic", "population") 
  uppers <- lapply(stats.list[[i]], function(x) {return(x[2,])})  # extract the "immune" data (mean + SEM)  of the group, across all pops 
  uppers.melt <- suppressMessages(melt(uppers))  # melt it
  names(uppers.melt) <- c("class", "upper", "population") 
  lowers <- lapply(stats.list[[i]], function(x) {return(x[1,])})  # extract the "immune" data (mean - SEM)  of the group, across all pops 
  lowers.melt <- suppressMessages(melt(lowers))  # melt it
  names(lowers.melt) <- c("class", "lower", "population") 
  # merge the three together 
  MERGED <- merge(means.melt, lowers.melt)
  MERGED <- merge(MERGED,uppers.melt)
  stats.obs.melt[[i]] <- MERGED  # assign into the list of lists 
  
  stats.obs.melt[[i]]$population <- factor(stats.obs.melt[[i]]$population, levels = c("NAM", "FL", "PAC", "ATL")) # change the order of pops (for plotting)
  stats.obs.melt[[i]]$population <- revalue(stats.obs.melt[[i]]$population, c("NAM" = "North America", "FL" = "South Florida", "PAC" = "Pacific", "ATL" = "Atlantic"))  # change the abbreviations to full pop name 
  stats.obs.melt[[i]]$class <- revalue(stats.obs.melt[[i]]$class, c("All.immune" = "A", "Recognition" = "R", "Signalling" = "S", "Modulation" = "M", "Effector" = "E"))  # change the full class name to abbreviations
  
  # Ctl: the control data 
  stats.list.ctl[[i]] <- lapply(stats.list[[i]], function(x) {return(x[-c(1:3),])}) # extract the "control" data of the group, across all pops 
  stats.ctl.melt[[i]] <- suppressMessages(melt(stats.list.ctl[[i]]))  # melt it 
  names(stats.ctl.melt[[i]]) <- c("class", "statistic", "population") 
  stats.ctl.melt[[i]]$population <- factor(stats.ctl.melt[[i]]$population, levels = c("NAM", "FL", "PAC", "ATL"))  # change the order of pops (for plotting)
  stats.ctl.melt[[i]]$population <- revalue(stats.ctl.melt[[i]]$population, c("NAM" = "North America", "FL" = "South Florida", "PAC" = "Pacific", "ATL" = "Atlantic"))  # change the abbreviations to full pop name 
  stats.ctl.melt[[i]]$class <- revalue(stats.ctl.melt[[i]]$class, c("All.immune" = "A", "Recognition" = "R", "Signalling" = "S", "Modulation" = "M", "Effector" = "E"))  # change the full class name to abbreviations
}


## Create a function for making grouped violin plots ----
immune.viol <- function(stat, title=NULL) {
  # split the input name and extract the two 
  testname <- strsplit(stat, "_")[[1]][1]
  site <- strsplit(stat, "_")[[1]][2]
  
  # get plot number 
  if (site != "full"){
    plot_list <- c("pi_0f", "pi_4f", "wtheta_0f", "wtheta_4f", "tjd_0f", "tjd_4f")  # a list of groups   
    } else {
    plot_list <- c("pi_full","wtheta_full", "tjd_full")  # a list of groups 
    }
   plot_num <- LETTERS[which(plot_list == stat)]
  
  # decode the names 
   if (testname == "tjd"){
     TEST_NAME = expression(bold("Tajima's D"))
   } else if (testname == "wtheta"){
     TEST_NAME = expression(bold(paste("Watterson's ", theta)))
   } else if (testname == "pi"){
     TEST_NAME = expression(bold(paste("Nucleotide diversity (", pi, ")")))
   } 
   
  
  if (site == "0f"){
     SITE = "0-fold sites"
   }  else if (site == "4f"){
     SITE = "4-fold sites"
   } else if (site == "full"){
     SITE = ""
   }  
  
  # set plot title 
  if(is.null(title)) {title <- SITE}
  # ggplot2 
  gp <- ggplot(data = stats.ctl.melt[[stat]], mapping= aes(x = class, y=statistic)) + 
    geom_violin(draw_quantiles = 0.5, fill = "azure3", alpha = 0.3, size = 0.3) +
    facet_grid(. ~ population) + 
    geom_point(data = stats.obs.melt[[stat]], color = "darkorange2", size = 1.5) +     # mean values 
    geom_segment(data = stats.obs.melt[[stat]], aes(x = class, y = lower, xend=class, yend = upper), color="darkorange2", size = 1.0, alpha = 0.6) +   # mean +- SE
    labs(title = title, y = TEST_NAME, x = "") +
    theme_bw(base_size = 6) + 
    theme(plot.title = element_text(size = 14, face = "bold"), axis.title.y=element_text(size=10,face="bold"), axis.text.x = element_text(size = 11, face = "bold"), axis.text.y = element_text(size = 9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.3)
          , strip.text.x = element_text(size=11,face="bold"), strip.background = element_rect(size = 0.3), 
          legend.position="none")
  gp <- arrangeGrob(gp, top = textGrob(paste0("(", plot_num, ")"), 
                                       x = unit(0.02, "npc"), y = unit(0.2, "npc"), gp = gpar(fontsize = 16, fontface = "bold")))
  return(gp)
}


## Create the multi-pannel plot ---- 
# Figure 6
plot_list <- c("pi_0f", "pi_4f", "wtheta_0f", "wtheta_4f", "tjd_0f", "tjd_4f")  # a list of groups 
PLOTS <- lapply(plot_list, function(x) {immune.viol(stat = x)})  # run all plots 
tiff(filename = "FIGURE6.tif",  width = 2400, height = 4000, res = 300)  # set output params 
grid.newpage()
plot1 <- ggarrange(PLOTS[[1]], PLOTS[[2]], PLOTS[[3]],  PLOTS[[4]],  PLOTS[[5]],  PLOTS[[6]], ncol = 1, nrow = 6,  align="hv")  # plot in grids
grid.draw(plot1)
dev.off()

# Figure S2 
plot_list <- c("pi_full","wtheta_full", "tjd_full")  # a list of groups 
PLOTS <- lapply(plot_list, function(x) {immune.viol(stat = x)})  # run all plots 
tiff(filename = "FIGURES2.tif",  width = 2400, height = 2400, res = 300)  # set output params 
grid.newpage()
plot2 <- ggarrange(PLOTS[[1]], PLOTS[[2]], PLOTS[[3]], ncol = 1, nrow = 3, align="hv")  # plot in grids
grid.draw(plot2)
dev.off()
