library(MEQLEA)
library(xtable)

setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Share/Simulation/Simulation20160318/") # the data set
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/simulation-results.R")

FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Manuscript/Figures/"
FigurePath <- getwd()
showcol <- c(1,4:9)
textsize = c(15,15,15, 15)

cases <- c("a","b", "c", "d", "e")
DEProb1 <- "0PCT"
DEProb2 <- "10PCT"

for ( i in 1: length(cases)){
  
  case1 <- cases[i]
  folder <- "Fixed1"
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


## Power Simulation -----------------------------------------------------------
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/Simulation20160318/")

p1 <- plotPower(folder= "Fixed1");
ggsave(paste(FigurePath, "powerBack0pct.eps", sep = ""), p1, width = 8, height = 5)
p2 <- plotPower(folder = "Fixed1", background = "BACK10")
ggsave(paste(FigurePath, "powerBack10pct.eps", sep = ""), p2, width = 8, height = 5)

##  power table
createPowerTable(folder = "Fixed1", showcol =  showcol, background = "BACK0", case = "a")
createPowerTable(folder = "Fixed1", showcol =  showcol, background = "BACK10", case = "a")


dat1 <- read.table("Fixed1/TypeIerror_a_10PCT.txt", header=T)
dat2 <- read.table("Fixed1/Power_a_15VS10PCT.txt", header=T)
colSums(dat2<0.05)
power2 <-RecalibratePower(dat1, dat2, colnum = showcol)
power2

