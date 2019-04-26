#### Single-population plots #### 
# Fig. 2, Fig. 7, Fig. S1, Fig. S3. 
# Plotted using absolute values rather than differences between immune and controls 


## load libraries ----
library(ggplot2)
library(reshape2)
library(grid)
library(ggpubr)

## Make a function for plotting ---- 
PLOT_ALL <- function(inputs){
## extract params from the inputs 
pop <- inputs[1]   # population 
site <- inputs[2]  # sites 
test <- inputs[3]   # popgen test 
NUM <- as.numeric(inputs[4])  # for graphic purpose

## decode names 
# population 
if (pop == "NAM"){
  population = "North America"
}  else if (pop == "FL"){
  population = "Florida"
} else if (pop == "PAC"){
  population = "Pacific"
} else if (pop == "ATL"){
  population = "Atlantic"
}
# test name 		
if (test == "tjd"){
  testname = expression(bold("Tajima's D"))
}  else if (test == "fst_w"){
  testname = expression(bold(paste(bolditalic(F)[ST])))
} else if (test == "wtheta"){
  testname = expression(bold(paste("Watterson's ", theta)))
} else if (test == "pi"){
  testname = expression(bold(paste("Nucleotide diversity (", pi, ")")))
}
# site
if (site == "0f"){
  SITE = "0-fold"
}  else if (site == "4f"){
  SITE = "4-fold"
} else if (site == "full"){
  SITE = "All sites"
}  

## Import data 
data <- read.table(paste0("data_for_plotting/Paired_control_rawdata_for_", pop, "_", test, "_", site, ".txt"), header = T)  # the first row is the immune group; the rest of 10,000 are control groups 
head(data)
dim(data)
# Subset the immune group stats and create coordinates for plotting
Immune <- data[1:3, ]   # mean-SE, mean+SE, and mean of the immune group 
Immune_cords <- data.frame(x1 = seq(1:5), y1 = as.numeric(paste(Immune[1,])), 
                           x2 = seq(1:5), y2 = as.numeric(paste(Immune[2,])))

# Subset the control groups and change format 
controls <- data[-(1:3), ]  # remove the first row of immune group 
Controls <- melt(controls)  # convert 
colnames(Controls) <- c("group", "stat") 
  
# Create a df for geom_points
Immune_group <- as.data.frame(cbind(colnames(Immune), as.numeric(Immune[3,])))
Immune_group$V2 <- as.numeric(paste(Immune_group$V2))
colnames(Immune_group) <- c("group", "stat")
str(Immune_group)

## Plotting 
# set title 
if (site != "full"){
  if (test != "fst_w"){
  TITLE <- substitute(underline(X), list(X = paste(population, SITE)))
  } else {
    TITLE <- substitute(underline(X), list(X = paste0(population, "-North America ", SITE)))
  }
} else {
  if (test != "fst_w"){
  TITLE <- substitute(underline(X), list(X = paste(population)))
  } else {
    TITLE <- substitute(underline(X), list(X = paste0(population, "-North America ")))
  }
}
# plot 
Fig <- ggplot(Controls, aes(x = group, y = stat))+
  geom_violin(draw_quantiles = 0.5, fill = "azure3", alpha = 0.3)+
  geom_segment(data = Immune_cords, aes(x = x1, y = y1, xend=x2, yend = y2), color="darkorange2", size = 1.5, alpha = 0.6) +    # mean +- SE
  geom_point(data = Immune_group, color="darkorange2", size = 3.5) +
  labs(title = TITLE)+
  ylab(testname)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.3, size = 16), axis.text.x = element_text(size = 13, angle = 55, hjust = 1.0), axis.text.y = element_text(size = 13, margin = margin(0, 0, 0, 9)), axis.title.x=element_blank(),axis.title.y=element_text(size=15),text = element_text(size=15, colour="black"))

## modify labels
if (test != "fst_w"){
  if (!(NUM %in% 1:2)){
    Fig <- Fig + labs(title = "")  # no sub-title
  }
}
if (NUM %in% c(2,4,6)){
  Fig <- Fig + theme(axis.title.y = element_blank())  # no y-axis
}

return(ggplotGrob(Fig))
}


## Plot the results for each population or FST ----

## make a function to generate the plots 
OUTPUT_PLOT <- function(POP, SITEs, title){

## Create a matrix to set input params
# NAM doesn't have FST 
if (SITEs != "full"){  # main text figures (0X and 4x results)
  if (POP == "NAM"){
    combs <- matrix(0, 6, 4)
    combs[,1] <- POP
    combs[,2] <- c("0f", "4f")
    combs[,3] <- rep(c("pi", "wtheta", "tjd"), each = 2)
    combs[,4] <- 1:6
  } else if(POP == "FST"){
      combs <- matrix(0, 6, 4)
      combs[,1] <- rep(c("FL", "PAC", "ATL"), each = 2)
      combs[,2] <- c("0f", "4f")
      combs[,3] <- "fst_w"
      combs[,4] <- 1:6
  } else{
    combs <- matrix(0, 8, 4)
    combs[,1] <- POP
    combs[,2] <- c("0f", "4f")
    combs[,3] <- rep(c("pi", "wtheta", "tjd", "fst_w"), each = 2)
    combs[,4] <- 1:8
  }
}else {                    # suppl. resutls (full gene region)
  if (POP == "NAM"){
    combs <- matrix(0, 3, 4)
    combs[,1] <- POP
    combs[,2] <- c("full")
    combs[,3] <- c("pi", "wtheta", "tjd")
    combs[,4] <- "N/A"
  } else if(POP == "FST"){
    combs <- matrix(0, 3, 4)
    combs[,1] <- c("FL", "PAC", "ATL")
    combs[,2] <- c("full")
    combs[,3] <- "fst_w"
    combs[,4] <- "N/A"
  } else{
    combs <- matrix(0, 4, 4)
    combs[,1] <- POP
    combs[,2] <- c("full")
    combs[,3] <- c("pi", "wtheta", "tjd", "fst_w")
    combs[,4] <- "N/A"
  } 
}  


## Apply the function to the combs by row 
plots <- apply(combs, 1, PLOT_ALL)

## Create multi-pannel plot
# arrange the plot (NAM only has 6 pannels)
if (length(plots) == 6){
  All <- ggarrange(
    plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], 
    ncol=2, nrow = 3, align="hv", 
    labels = c("(A)", "(B)","(C)","(D)","(E)","(F)"))
} else if (length(plots) == 8){
  All <- ggarrange(
    plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], ncol=2, nrow = 4, align="hv", 
    labels = c("(A)", "(B)","(C)","(D)","(E)","(F)", "(G)", "(H)"))
} else if (length(plots) == 3){
  All <- ggarrange(
    plots[[1]], plots[[2]], plots[[3]],
    ncol=2, nrow = 2, align="hv", 
    labels = c("(A)", "(B)","(C)"))
} else if (length(plots) == 4){
  All <- ggarrange(
    plots[[1]], plots[[2]], plots[[3]], plots[[4]],
    ncol=2, nrow = 2, align="hv", 
    labels = c("(A)", "(B)","(C)","(D)"))
}

# Output the plot 
if (length(plots) > 4){
  tiff(filename = paste0("Plots/", title, ".tif"), width = 2400, height = 3600, res = 300)
} else {
  tiff(filename = paste0("Plots/", title, ".tif"), width = 2400, height = 2400, res = 300)
}
grid.newpage()
grid.draw(All)
if (SITEs != "full"){
  grid.lines(x = 0.51, y = c(0, 1), gp = gpar(col = "darkgray", lty = 2, lwd = 2)) 
}
dev.off()
}


## Run the function to output the plot ----
# the input for POP is either a population (e.g. NAM) or FST
OUTPUT_PLOT("NAM", "04f", "FIGURE2")
OUTPUT_PLOT("FST", "04f", "FIGURE7")
OUTPUT_PLOT("NAM", "full", "FIGURES1")
OUTPUT_PLOT("FST", "full", "FIGURES3")
