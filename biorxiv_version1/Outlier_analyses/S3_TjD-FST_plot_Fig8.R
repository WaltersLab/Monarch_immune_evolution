#### Plotting Tajima's D vs. Fst across all 4 populations #### 

## Import datasets ---- 
# import all results 
for (k in c("0f", "4f")){ # for 0f or 4f 
  site <- k 
  for (i in c("tjd", "fst_w") ){  # for all 2 tests 
    test_in <- i
    for (j in c("NAM", "FL", "PAC", "ATL")){  # for all 4 pops 
      pop = j
      if (!((test_in == "fst_w") & (pop == "NAM"))){  # exclude the NAM-FST that didn't exist 
        data <- read.table(paste(pop, "_", test_in, "_", site, "_all_IG_full_results.txt", sep = ""), header = T)
        # assign data
        assign(paste(pop, "_", test_in, "_", site, "_data", sep = ""), data) 
        rm(data)
      }	
    }
  }
}

## Subset by gene class ----
# make a function for combine datasets.
# input = test_1, test_2, type, class
# test_1 & test_2: the two tests for plotting. Must always keep FST only in test 2.  
# type: 2 for class, 3 for pathway, 4 for gene 
# class: # name of the class or pathway or gene
extact_combine <- function(test_1, test_2, type, class, site){
	for (i in 1:4){
		pop <- c("NAM", "FL", "PAC", "ATL")[i]  # pop 
		if (!((test_2 == "fst_w") & (pop == "NAM"))){  # if not involving the NAM-FST combination 
		  test_1_data <- get(paste(pop, "_", test_1, "_", site, "_data", sep = ""))  # get data for test 1 
		  test_2_data <- get(paste(pop, "_", test_2, "_", site, "_data", sep = "")) # get data for test 2 	
			D <- cbind(matrix(pop, nrow(test_1_data[which(test_1_data[, type] == paste(class)), ]), 1),  # set up matrix 
			test_1_data[which(test_1_data[, type] == paste(class)), ][, c(1,6,7,8)], # test 1 data 
			test_2_data[which(test_2_data[, type] == paste(class)), ][,6:8],  # test 2 data
			numeric(nrow(test_1_data[which(test_1_data[, type] == paste(class)), ]))) # outlier pch 
			colnames(D) <- c("pop", colnames(D)[2],"T1", "T1_std", "T1_outliers", "T2", "T2_std", "T2_outliers", "pch")
			assign(paste("D", i, sep = ""), D)  # assign name 
			rm(D)
		} else {  # for the NAM-FST combination 
		  test_1_data <- get(paste(pop, "_", test_1,"_", site, "_data", sep = ""))  # get data for test 1 
		  D <- cbind(matrix(pop, nrow(test_1_data[which(test_1_data[, type] == paste(class)), ]), 1),  # set up matrix 
			test_1_data[which(test_1_data[, type] == paste(class)), ][, c(1,6,7,8)], # test 1 data 
			matrix(0.5, nrow(test_1_data[which(test_1_data[, type] == paste(class)), ]), 1), # Fst data (set to 0.5 for NAM)
			matrix(0, nrow(test_1_data[which(test_1_data[, type] == paste(class)), ]), 1), # Fst std data (set to 0 for NAM)
			matrix(0, nrow(test_1_data[which(test_1_data[, type] == paste(class)), ]), 2))  # outlier pch 
			colnames(D) <- c("pop", colnames(D)[2],"T1", "T1_std", "T1_outliers", "T2", "T2_std", "T2_outliers", "pch")
			assign(paste("D", i, sep = ""), D)  # assign name 
			rm(D)
		}	
	}	
	DATA <- rbind(D1, D2, D3, D4)  # combine the four
	
	# Assign outlier pch for graphing - outliers defined as being outlier in either D or Fst or both 
	for (i in 1:nrow(DATA)){
		if (DATA[i,1] == "NAM"){
			DATA[i, 9] <- ifelse(DATA[i, 5] + DATA[i, 8] > 0, 16, 1)  # if either outlier in Tjd or Fst, solid point
		}
		else if (DATA[i,1] == "FL"){
			DATA[i, 9] <- ifelse(DATA[i, 5] + DATA[i, 8] > 0, 15, 0) 
		}
		else if (DATA[i,1] == "PAC"){
			DATA[i, 9] <- ifelse(DATA[i, 5] + DATA[i, 8] > 0, 18, 5) 
		}
		else if (DATA[i,1] == "ATL"){
			DATA[i, 9] <- ifelse(DATA[i, 5] + DATA[i, 8] > 0, 17, 2) 
		}
	}    
	return(DATA)	
}	


## Tajima's D & Fst 2D-plot ----  
# Make a functon for 2D plot and report outlier info
Plot_and_outliers <- function(DATA, NUM){
  # decode test name 
  if (test_1 == "tjd"){
    test_1_name = "Tajima's D"
  }  else if (test_1 == "fst_w"){
    test_1_name = "FST"
  } else if (test_1 == "wtheta"){
    test_1_name = "Watterson's theta"
  } else if (test_1 == "pi"){
    test_1_name = "Pi"
  }
  if (test_2 == "tjd"){
    test_2_name = "Tajima's D"
  }  else if (test_2 == "fst_w"){
    test_2_name = "FST"
  } else if (test_2 == "wtheta"){
    test_2_name = "Watterson's theta"
  } else if (test_2 == "pi"){
    test_2_name = "Pi"
  }
  SITE <- ifelse(site == "4f", "4-fold sites", "0-fold sites")
  # set output and params
  par(mar=c(5,6,2.5,0.5))  # set graph par
  # plotting 
  plot(DATA$T1, DATA$T2, 
       xlab = ifelse(NUM %in% 7:8, expression(bold("Tajima's D Percentile in Background")), ""), 
       ylab = ifelse(NUM %in% c(1,3,5,7), expression(bold(bolditalic(F)[ST]~"Percentile in Background")), ""),
       xlim = c(0, 1), ylim = c(0, 1), type = "n",
       cex.main = 1.8, cex.lab = 1.2, cex.axis = 1.1)  # emtpy plot 
  title(paste(class, " genes: ", SITE, sep = ""), adj = 0.02, line = 1, cex.main = 1.4)
  segments(-0.04, 0.5, 1.04, 0.5, lty = 2, col = "black", lwd = 1.3)   # mean = 0 
  segments(0.5, -0.04, 0.5, 1.04, lty = 2, col = "black", lwd = 1.3)   # mean = 0 
  segments(0.025, 0.025, 0.975, 0.025, lty = 3, col = "grey50", lwd = 2)   # the 4 below are for 5% and 95% 
  segments(0.025, 0.025, 0.025, 0.975, lty = 3, col = "grey50", lwd = 2)
  segments(0.975, 0.025, 0.975, 0.975, lty = 3, col = "grey50", lwd = 2)
  segments(0.025, 0.975, 0.975, 0.975, lty = 3, col = "grey50", lwd = 2)
  points(DATA$T1, DATA$T2, pch = DATA$pch, lwd = 1.5, cex = ifelse(DATA$pch > 10, 1.2, 1.0),
         col = c("lightgoldenrod4", "red", "cyan3", "purple" )[DATA$pop])   # add points
  mtext(paste0("(", LETTERS[as.numeric(NUM)], ")"), side = 3, adj = -0.2, line = 1, cex = 1.1) 
  outlier_genes <- DATA[which(DATA$pch > 10),][, c(1, 2, 3, 5, 6, 8)]
  write.table(outlier_genes, file = paste(class, "_", test_1, "_", test_2,"_", site, "_outliers.txt", sep = ""), 
               row.names = F, col.names = T)
}


### run the function ----

tiff("FIGURE8.tif", width = 2500, height = 6200, res = 400)
par(mfrow = c(5,2))  # set frame 

# create a matrix for all combinations
combs <- matrix(0, 8, 3)  
combs[, 1] <- rep(c("Recognition", "Signalling", "Modulation", "Effector"), each = 2)  # class
combs[, 2] <- rep(c("0f", "4f"), 4)  # site 
combs[, 3] <- 1:8 # num (for plotting purpose) 

## set input via loop and run it
for (i in 1:8){
  test_1 = "tjd"
  test_2 = "fst_w"
  type = 2
  class = combs[i, 1]
  site = combs[i, 2]
  NUM = combs[i, 3]
  DATA <- extact_combine(test_1, test_2, type, class, site)  # run function 1
  Plot_and_outliers(DATA, NUM)  # run function 2 
}

## Add legend 
plot.new()  # add two blank plots 
plot.new()
legend("topright", c("North Ameirca", "South Florida", "Pacific", "Atlantic"), 
       col = c("lightgoldenrod4", "red", "cyan3", "purple"),
       pch = c(1, 0, 5, 2), cex = 1.2, ncol = 2)
dev.off()
