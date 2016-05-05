source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/R/MOMEnrichmentTest.R")                                          # the code 
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/R/SimulateLabDataNew.R")                                          # the code 
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/R/SimulationResults.R") #
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/Simulation20160226/") # the data set
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/GSEA.1.0.R")

## Type I error Simulation ----------------------------------------------------


FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Manuscript/Figures/"

showcol <- c(1,4:9)
textsize = rep(15,15,15, 15)

cases <- c("a0", "a", "c", "e", "f")
updated_case <- c("a", "b", "c", "d", "e") 
DEProb1 <- "0PCT"
DEProb2 <- "10PCT"
i <- 1
for ( i in 1: length(cases)){
  
  case <- cases[i]
  case1 <- updated_case[i]
  folder <- "SizePoint1"
  d1 <- paste(folder, "/TypeIerror_", case, "_", DEProb1, ".txt", sep="")
print(d1)
    d2 <- paste(folder, "/TypeIerror_", case, "_", DEProb2, ".txt", sep="")
  #title <- paste(folder, ", 10PCT, ", case, sep="")
  #title <- ""
 # figname <- paste(case,"_", folder, "_", DEProb, ".eps", sep = "")
  #d1 <- "SizePoint1/Power_a_10VS0PCT.txt"
  dat1 <- read.table(d1)
  dat2 <- read.table(d2)
  colnames(dat1) <- colnames(dat2) <- c("MEQLEA", "LM", "GeneSetTest-modt", "MRSGE", "CAMERA-modt",
                     "CAMERA-rank", "GSEA", "QuSAGE")
 # create.hist2(dat[, showcol], textsize = textsize, figure.num =title, figname = figname, print.figure= F)
#  qqTypeI(dat[, showcol], textsize = textsize, figure.num =title, figname = figname, print.figure= F)
  # dat <- read.table("SizePoint1/TypeIerror_a_0PCT_SizePoint1_2.txt")
  # qqunif(dat$GSEA)
  fig.obj <- ArrangeTypeIerror(dat1[, showcol], dat2[,showcol],case = case1, legend=F)
   ggsave(paste(FigurePath,"/parallel_", case, ".eps", sep =""), fig.obj, 
          width = 16, height = 5)
  legend_need <- ArrangeTypeIerror(dat1[, showcol], dat2[, showcol], case=case1, legend=T)
  ggsave(paste(FigurePath,"/parallel_legend.eps", sep =""), legend_need, 
         width = 16, height = 0.6)
}

 ## Power Simulation -----------------------------------------------------------
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/Simulation20160226/")

p1 <- plotPower(0.1);
ggsave(paste(FigurePath, "powerBack0pct.eps", sep = ""), p1, width = 8, height = 5)
p2 <- plotPower(0.1, background = "BACK10")
ggsave(paste(FigurePath, "powerBack10pct.eps", sep = ""), p2, width = 8, height = 5)

##  power table
createPowerTable(size = 0.1, showcol =  showcol, background = "BACK0", case = "a0")
createPowerTable(size = 0.1, showcol =  showcol, background = "BACK10", case = "a0")


dat1 <- read.table("SizePoint1/TypeIerror_a0_0PCT.txt")
dat2 <- read.table("SizePoint1/Power_a0_10VS0PCT.txt")
power2 <-RecalibratePower(dat1, dat2, colnum = showcol)


