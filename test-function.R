# this is the test funciton 

library(MEQLEA)

######### real data 


a0 <- read_gene_set("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/DataSet/Gender/C1.gmt")
CLS <- GSEA.ReadClsFile(file="/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/DataSet/Gender/Gender.cls")
dataset <- GSEA.Gct2Frame2(filename = "/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/DataSet/Gender/Gender.gct")
geneset <- read_gene_set("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/DataSet/Gender/C1.gmt")

results <- meqlea_multiple(dataset, CLS$class.v, geneset, minSetSize=2)
 
 
## test the data 

library(doParallel)                             # parallel computing packages
library(foreach)

cl <- makeCluster(30)                            # request 3 cores to do the simulation
registerDoParallel(cl)                           # Register cluster

case <- "a"
files <- "TypeIerror_e_0PCT.txt"
destination <-  paste("/home/stats/zhuob/data/computing/", files, sep = "")

nsim <- 100
size <- 50           # number of samples to be simulated

# rho <- c(0.7, 0.5, -0.5)   # extreme case
rho <- c(0.1, 0.05, -0.05) # correlation for case a, e, f, used in simulation
num_gene <- c(500, 100)
prop <- c(0.1, 0.1)
n_gene <- num_gene[1]
delta <- rnorm(n_gene, 2, 1)

fti <- foreach(i = 1:nsim, .combine = rbind, .packages = c("MASS", "qusage")) %dopar% {

source("/home/stats/zhuob/Rcode/Enrichment/MOMEnrichmentTest.R")
source("/home/stats/zhuob/Rcode/Enrichment/SimulateLabDataNew.R")               # it contains group.mean function
source("/home/stats/zhuob/Rcode/Enrichment/GSEA.1.0.R")

   obj <- prepare.simulation(num_gene, prop, delta, case = case, rho)
   #obj <- MEQLEA::prepare_simulation(num_gene, prop, delta, case = case, rho)

  dat <- simulate.microarry(size, obj)
  # dat <- MEQLEA::simulate_expression_data(size, obj)
    
  # cat("\r", i)
  # 
	
  pvals <- compare_test(dat)
#  pvals <- MEQLEA::compare_test(dat)
return(pvals)
}


colnames(fti)[1:9] <- c( "OurTest", "OurTest.OneSided",  "lm",  "ModeratedT",  "MRSGE", "Camera", "CameraRank", "GSEA", "QUSAGE")
colMeans(fti <0.05)

