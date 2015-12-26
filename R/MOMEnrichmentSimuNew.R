source("/home/stats/zhuob/Rcode/Enrichment/MOMEnrichmentTest.R")
source("/home/stats/zhuob/Rcode/Enrichment/SimulateLabDataNew.R")               # it contains group.mean function

library(MASS)
#library(dplyr)
library(doParallel)                             # parallel computing packages
library(foreach)

cl <- makeCluster(30)                            # request 3 cores to do the simulation
registerDoParallel(cl)                           # Register cluster


case <- "g"
destination <-  "/home/stats/zhuob/data/computing/Power_g_50_New.txt"


	nsim <- 10000
  
	size <- 50           # number of samples to be simulated
	rho <- c(0.1, 0.05, -0.05)   # correlation for case a, e, f
	num_gene <- c(500, 100)
	prop <- c(0.2, 0.0)
	n_gene <- num_gene[1]
	delta <- rnorm(n_gene, 2 , 1)
  
  
  
fti <- foreach(i = 1:nsim, .combine = rbind) %dopar% {

	
        if (case != "g"){
                obj <- prepare.simulation(num_gene, prop, delta, case = case, rho)
        }
	if (case == "g"){
		library(dplyr) 
               	obj <- real.go.term(2, num_gene, prop, delta)
                }

    
	library(MASS)

     	dat <- simulate.microarry(size, obj)
    

	# first standardize the expression data  
	# 1. remove the group mean, and get the sd for each row
	# 2. scale original data by dividing the sd for each row

	dat_sd <- standardize.microarray(dat$data, dat$trt)
    	dat$data <- dat_sd
	
	pvals <- compare_test(dat)
    	return( pvals)
  }
  
  
  colnames(fti)[1:6] <- c( "OurTest",  "lm",  "ModeratedT",  "Rank", "Camera", "CameraRank")
  
  
  write.table(fti, destination)
