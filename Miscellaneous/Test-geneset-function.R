

case <- "a0"
files <- "TypeIerror_e_0PCT.txt"



destination <-  paste("/home/stats/zhuob/data/computing/", files, sep = "")


nsim <- 10

size <- 50           # number of samples to be simulated

# rho <- c(0.7, 0.5, -0.5)   # extreme case
rho <- c(0.1, 0.05, -0.05) # correlation for case a, e, f, used in simulation	
num_gene <- c(500, 100)
prop <- c(0.0, 0.0)	
n_gene <- num_gene[1]
delta <- rep(0.1, n_gene)
# delta <- rnorm(n_gene, 0.2 , 0.1)  # make the delta positive
r1 <- r2 <- c()

for ( i in 1: nsim){
  
  obj <- prepare.simulation(num_gene, prop, delta, case = case, rho)
  
  dat <- simulate.microarry(size, obj)
  microarray <- dat$data
  trt <- dat$trt
  
  design <- model.matrix(~trt) 
  fit <- lmFit(microarray, design)
  fit <- eBayes(fit)								                # Emperical Bayes t test	
  stat1 <- fit$t[, 2]            							      #	use the moderated t statistics to do enrichment
  stat2 <- fit$coefficients[, 2]
  r1[i] <- cor(stat1, stat2)
  stat3 <- group.mean(microarray, trt)
  stat4 <- stat3[, 50] - stat3[, 1]
  r2[i] <- cor(stat1, stat4)
}




















library(MASS)
#library(dplyr)
library(doParallel)                             # parallel computing packages
library(foreach)

cl <- makeCluster(2)                            # request 3 cores to do the simulation
registerDoParallel(cl)                           # Register cluster


case <- "a"
files <- "TypeIerror_a0_0PCT_either.txt"


pvals <- matrix(NA, 100, 4)
destination <-  paste("/home/stats/zhuob/data/computing/", files, sep = "")


nsim <- 100

size <- 50           # number of samples to be simulated

# rho <- c(0.7, 0.5, -0.5)   # extreme case
rho <- c(0.1, 0.05, -0.05) # correlation for case a, e, f, used in simulation
num_gene <- c(500, 100)
prop <- c(0.2, 0.1)

n_gene <- num_gene[1]
delta <- rep(0.1, n_gene)
# delta <- rnorm(n_gene, 0.2 , 0.1)  # make the delta positive


for ( i in 1: 100){
  if (case != "g"){
    obj <- prepare_simulation(num_gene, prop, delta, case = case, rho)
  }
  if (case == "g"){
    library(dplyr)
    obj <- real_go_term(2, num_gene, prop, delta)
  }
  
  library(MASS)
  library(qusage)
  
  print(i)
  dat <- simulate_expression_data(size, obj)
  
  #pvals <- compare_test(dat)
  
  library(limma)
  
  microarray <- dat$data
  trt <- dat$trt
  go_term <- dat$go_term
  
  #	pval1 <- MOM_test(microarray, trt, go_term, standardize=F)$p             # our test, NO standardization for the si$
  
  ## modified on March 16
  MEQ <- meqlea_single(microarray, trt, go_term, standardize=F)             # our test, NO standardization for the simula$
  #	print(MEQ)
  pval1 <- MEQ$p1                                                      # two sided test
  pval1.2 <- MEQ$p2                                                  # one sided test
  
  design <- model.matrix(~trt)
  fit <- lmFit(microarray, design)
  fit <- eBayes(fit)                                                                              # Emperical Bayes $
  stat <- fit$t[, 2]                                                                    # use the moderated t statis$
  # stat <- fit$coefficients[, 2]
  
  alter <- "either"                                  # This is the default option, see below for explanation.
  index1 <- which(go_term==1)                       # which genes are in GOTERM
  
  
  
  # up: positive statistics
  # down: negative stats
  # either: the set is either up or down-regulated as a whole
  # mixed: tests whether the genes in the set tend to be DE without regard for direction.
  #        In this case, the test will be significant if the set contains mostly large statistics, negative or positive
  
  tes1 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= F)         # moderated t GeneSet
  tes2 <- sig_path(index1, stat)
  pvals <- c(pval1, pval1.2, tes1, tes2)
  
  print(pvals)
  
 # return( pvals)
}



#  colnames(fti)[1:8] <- c( "OurTest",  "lm",  "ModeratedT",  "MRSGE", "Camera", "CameraRank", "GSEA", "QUSAGE")

# modified on March 16
colnames(fti)[1:4] <- c( "OurTest", "OurTest.OneSided", "ModeratedT", "ModerateT2")