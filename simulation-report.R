library(MEQLEA)
library(xtable)

source("simulation-results.R")
setwd("Share/Simulation/Simulation20160318/") # the data set


FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Manuscript/Figures/"
FigurePath <- getwd()
showcol <- c(1,4:9)
textsize = c(15,15,15, 15)

cases <- c("a","b", "c", "d", "e")
DEProb1 <- "0PCT"
DEProb2 <- "10PCT"

for ( i in 1: length(cases)){
  
  case1 <- cases[i]
  folder <- "Normal1"
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
setwd("Share/Simulation/Simulation20160318/")

p1 <- plotPower(folder= "Normal2");
ggsave(paste(FigurePath, "powerBack0pct.eps", sep = ""), p1, width = 8, height = 5)
p2 <- plotPower(folder = "Normal2", background = "BACK10")
ggsave(paste(FigurePath, "powerBack10pct.eps", sep = ""), p2, width = 8, height = 5)

# what about other methods?
dat1 <- read.table("Fixed1/TypeIerror_a_10PCT.txt", header=T)
names(dat1)
method <- 3
p3 <- plotPower(folder="Fixed1", showcol = method); p3
p4 <- plotPower(folder="Fixed1", showcol = method, background = "BACK10"); p4

##  power table
createPowerTable(folder = "Normal1", showcol =  showcol, background = "BACK0", case = "a")
createPowerTable(folder = "Normal1", showcol =  showcol, background = "BACK10", case = "a")



showcol <- c(1, 4:9)
dat1 <- read.table("Fixed1/TypeIerror_a_0PCT.txt", header=T)
dat2 <- read.table("Fixed1/Power_a_5VS0PCT.txt", header=T)
#colSums(dat2<0.05)
power2 <-RecalibratePower(dat1, dat2, colnum = showcol)
power2

dat3 <- read.table("Fixed1/TypeIerror_c_0PCT.txt", header=T)
dat4 <- read.table("Fixed1/Power_c_5VS0PCT.txt", header=T)
power3 <-RecalibratePower(dat3, dat4, colnum = showcol)
power3
rbind(power2$power, power3$power)
round(rbind(power2$adjusted_alpha, power3$adjusted_alpha),3)


method <- 9
summarizePower("Fixed1", background = "BACK10", coln= method, case = "a")[2, 3]
summarizePower("Fixed1", background = "BACK10", coln= method, case = "c")[2, 3]



