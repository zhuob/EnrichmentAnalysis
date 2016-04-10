library(MEQLEA)
library(xtable)

source("simulation-results.R")
setwd("Share/Simulation/Simulation20160318/") # the data set
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/")

FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Manuscript/Figures/"
FigurePath <- getwd()
showcol <- c(1,4:9)
textsize = c(15,15,15, 15)

cases <- c("a","b", "c", "d", "e")
DEProb1 <- "0PCT"
DEProb2 <- "10PCT"
folder <- "Simulation"

for ( i in 1: length(cases)){
  
  case1 <- cases[i]

  d1 <- paste(folder, "/TypeIerror_", case1, "_", DEProb1, ".txt", sep="")
  print(d1)
  d2 <- paste(folder, "/TypeIerror_", case1, "_", DEProb2, ".txt", sep="")
  print(d2)
 
  dat1 <- read.table(d1, header=T)
  dat2 <- read.table(d2, header=T)
  colnames(dat1) <- colnames(dat2) <- c("MEQLEA", "MEQLEA-N", "LM", "MRSGE", "SigPath", "CAMERA-modt",
                                        "CAMERA-rank", "GSEA", "QuSAGE")

  fig.obj <- ArrangeTypeIerror(dat1[, showcol], dat2[,showcol],case = case1, legend=F)
  ggsave(paste(FigurePath,"/parallel_", case1, ".eps", sep =""), fig.obj,
         width = 18, height = 4.5)
  legend_need <- ArrangeTypeIerror(dat1[, showcol], dat2[, showcol], case=case1, legend=T)
  ggsave(paste(FigurePath,"/parallel_legend.eps", sep =""), legend_need,
         width = 20, height = 0.6)
}



### type I error 



## Power Simulation -----------------------------------------------------------
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Share/Simulation/Simulation20160318/")
## power for figure
p1 <- plotPower(folder= folder);
ggsave(paste(FigurePath, "powerA1pct.eps", sep = ""), p1, width = 8, height = 5)
p2 <- plotPower(folder = folder, setup = "A2")
ggsave(paste(FigurePath, "powerA2pct.eps", sep = ""), p2, width = 8, height = 5)


##  power table
library(xtable)
createPowerTable(folder = "Normal1", showcol =  showcol, setup = "A1", case = "a")
createPowerTable(folder = "Normal1", showcol =  showcol, setup = "A2", case = "a")

