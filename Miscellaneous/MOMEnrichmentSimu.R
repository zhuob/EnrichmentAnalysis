## this file includes all the test procedure and Simulation of microarray data
source("/home/stats/zhuob/Rcode/Enrichment/MOMEnrichmentTest.R")
source("/home/stats/zhuob/Rcode/Enrichment/SimulateLabData.R")               # it contains group.mean function


library(MASS)
library(doParallel)                             # parallel computing packages
library(foreach)

cl <- makeCluster(30)                            # request 3 cores to do the simulation
registerDoParallel(cl)                           # Register cluster


case <- "c"
destination <-  "/home/stats/zhuob/data/computing/NO_DE_c_50.txt"

nsim <- 1000
size <- 50
n_gene <- 500
go_prop <- 0.0
no_go_prop <- 0.0
delta_mat <- 2



fti <- foreach(i = 1:1000, .combine = rbind) %dopar% {

	dat <- simulate.microarry(size, n_gene, go_prop, no_go_prop, delta_mat, case = case)
 	pvals <- compare_test(dat)
	return( pvals)
}


colnames(fti)[1:6] <- c( "OurTest",  "lm",  "ModeratedT",  "Rank", "Camera", "CameraRank")


write.table(fti, destination)

