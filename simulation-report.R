library(MEQLEA)
library(xtable)
library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Manuscript/Figures/"
FigurePath <- getwd()
showcol <- c(1,4:9)
textsize = c(15,15,15, 15)

cases <- c("a","b", "c", "d", "e")
DEProb1 <- "0PCT"
DEProb2 <- "10PCT"
folder <- "Results/Simulation"

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



### type I error table

alpha <- 0.01
table_typeIerror <- data.frame(matrix(NA, nrow = 1000, ncol = length(showcol)+1))
for ( i in 1: length(cases)){
  
  case1 <- cases[i]
  
  d1 <- read.table(paste(folder, "/TypeIerror_", case1, "_", DEProb1, ".txt", sep=""), header=T)
  ts1 <- typeIerror_quantile(d1, alpha, showcol)
  d2 <- read.table(paste(folder, "/TypeIerror_", case1, "_", DEProb2, ".txt", sep=""), header=T)
  ts2 <- typeIerror_quantile(d2, alpha, showcol)
  table_typeIerror[2*i-1, -1] <- ts1
  table_typeIerror[2*i, -1] <- ts2
  names(table_typeIerror) <- c("case", names(ts1))
  table_typeIerror[c(2*i-1, 2*i), 1] <- paste(case1, c(DEProb1, DEProb2), sep="")
  }

table_typeIerror <- table_typeIerror[complete.cases(table_typeIerror), ]
table_typeIerror

test_sig <- function(prop, alpha, n){
  sd <- sqrt((alpha)*(1-alpha)/n)
  z <- (prop-alpha)/sd
  2*pnorm(abs(z), lower.tail = F)
}

significance <- table_typeIerror
for(j in 1:(ncol(table_typeIerror)-1) ){
  for ( i in 1:nrow(table_typeIerror)){
    significance[i, j+1] <- test_sig(table_typeIerror[i, j+1], alpha, 10000)
    }
  }

## Power Simulation -----------------------------------------------------------
#setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Share/Simulation/Simulation20160318/")
## power for figure
p1 <- plotPower(folder= folder);
ggsave(paste(FigurePath, "/powerA1pct.eps", sep = ""), p1, width = 8, height = 5)
p2 <- plotPower(folder = folder, setup = "A2")
ggsave(paste(FigurePath, "/powerA2pct.eps", sep = ""), p2, width = 8, height = 5)

folder <- "Fixed1"
##  power table
library(xtable)
t1 <- createPowerTable(folder = folder, showcol =  showcol, setup = "A1", case = "a")
t2 <- createPowerTable(folder = folder, showcol =  showcol, setup = "A2", case = "a")
xtable(t1, digits = c(0, 3,3,3, 3, 3))
xtable(t2, digits = c(0, 3,3,3, 3, 3))






