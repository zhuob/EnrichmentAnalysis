source("/home/stats/zhuob/Rcode/Enrichment/MOMEnrichmentTest.R")
source("/home/stats/zhuob/Rcode/Enrichment/SimulateLabDataNew.R")               # it contains group.mean function
source("/home/stats/zhuob/Rcode/Enrichment/GSEA.1.0.R")


library(MASS)
#library(dplyr)
library(doParallel)                             # parallel computing packages
library(foreach)

cl <- makeCluster(30)                            # request 3 cores to do the simulation
registerDoParallel(cl)                           # Register cluster


case <- "e"
files <- "Power_e_30VS10PCT.txt"



destination <-  paste("/home/stats/zhuob/data/computing/", files, sep = "")


	nsim <- 10000
  
	size <- 50           # number of samples to be simulated
	
	# rho <- c(0.7, 0.5, -0.5)   # extreme case
	 rho <- c(0.1, 0.05, -0.05) # correlation for case a, e, f, used in simulation	
	num_gene <- c(500, 100)
	prop <- c(0.30, 0.10)	
	n_gene <- num_gene[1]
	delta <- rep(0.1, n_gene)
	# delta <- rnorm(n_gene, 0.2 , 0.1)  # make the delta positive
  
  
  
fti <- foreach(i = 1:nsim, .combine = rbind) %dopar% {

	
        if (case != "g"){
                obj <- prepare.simulation(num_gene, prop, delta, case = case, rho)
        }
	if (case == "g"){
		library(dplyr) 
               	obj <- real.go.term(2, num_gene, prop, delta)
                }

    
	library(MASS)
	library(qusage)
	
     	dat <- simulate.microarry(size, obj)
    

	pvals <- compare_test(dat)
    	return( pvals)
  }
  
  
  colnames(fti)[1:8] <- c( "OurTest",  "lm",  "ModeratedT",  "MRSGE", "Camera", "CameraRank", "GSEA", "QUSAGE")
  
  
  write.table(fti, destination)
