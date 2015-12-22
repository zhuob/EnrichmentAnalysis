## evaluate the test statistic for enrichment
# I create two cases: 1. independent t's; 2. correlated t's with true parameters

 
source("/home/stats/zhuob/Rcode/Enrichment/SimulateLabData.R")               # it contains group.mean function
source("/home/stats/zhuob/Rcode/Enrichment/OurEnrichmentTest.R")                # the solve function and test function

library(MASS)


library(doParallel)                             # parallel computing packages
library(foreach)

cl <- makeCluster(24)                            # request 3 cores to do the simulation
registerDoParallel(cl)                          # Register cluster

xi_pre <- seq(0.001, 1, length=1000)
if (INDICATOR == "A") {

system.time(
  fti <- foreach(i = 1:1000, .combine = rbind) %dopar% {

################  Set up the parameters  -----------------
  library(MASS)	
  size <- 50                                     # number of samples to be considered
  n_gene <- 500                                  # number of genes to be considered
#  de_prop <- 0.2                                 # proportion of DE
#  de <- n_gene*de_prop                           # number of DE genes
#  delta <- rexp(n_gene, 0.2)                     # magnitude of DE
#  z <- rbinom(n_gene, 1, de_prop)                # indicator, DE =1, NON_DE=0 
#  delta <- delta*z
#  delta <- rep(0, n_gene)             ## NULL, NO TREATMENT effect


######### prepare go term ---------------------------------
#  n_go <- 0.2*n_gene                                       # number of genes in the go term
#  go_term <- rep(0, n_gene)                                # go term indicator
#  go_prop <- 0.6                                           # proportion of DE genes in the go term
#  go_de <- sample(which(z==1), n_go*go_prop)               # de genes in go term
#  go_not_de <- sample(which(z==0), n_go*(1-go_prop))	 # non DE genes in go term
#  go_in <- c(go_de, go_not_de)
#  go_term[go_in] <- 1     

go_term <- rbinom(n_gene, 1, 0.2)                        # random sample of go term
################# To generate the correlation matrix --------------
  rho <- 0.5                                     # correlation matrix
  sigma0 <- 2                                    # variance part 
  cor.struct <- matrix(rho, n_gene, n_gene);    
  diag(cor.struct) <- 1                          # correlation matrix

    sigma <- 1.5*cor.struct
    dat <- mvrnorm(1000, mu=rep(0, n_gene), sigma)
    
    sigma1 <- cor(dat)
    sigma2 <- 1.5*sigma1
    dat <- mvrnorm(1000, mu=rep(0, n_gene), sigma2)
    
    r1 <- cor(dat)                               # Covariance structure
#   r1 <- diag(1, n_gene)                        #  Make the t_vals independent of each other
   
  xi <- xi_pre[i]   
  sigma <- sigma0*( (1-xi)*diag(1, n_gene) + xi*r1)                     # covariance matrix for both groups
#  sigma <- sigma0*diag(n_gene)       ##  Make the t_vals independent of each other

############# generate expression matrix  ------------------------


  mu_initial <- 0
  mu1 <- rep(mu_initial, n_gene)                 # mean expression value for control
#  mu2 <- delta  + mu_initial                     # mean expression value for treatment
  mu2 <- mu1

  microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix
  n <- size/2
  microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma)) 
  microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))

  trt <- rep(c(0, 1), each= size/2)              # group indicator


############### prepare the t_stat, samp_rho and Go term -------------- 

#########  t_test & sample correlation
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix

  samp_rho <- cor(t(resid_mat))                            # sample correlation matrix
  t_val <- group_mean[, size] - group_mean[, 1]            # control - treatment
  t_val_sim <- microarray[, 1]                             # simulate the true t_val

###*********** looks better if I use the true correlation matrix

######## begin the test procedure ---------------------------------------
  solve_eq <- solve.equations(t_val, samp_rho)             # estimate the xi, beta and sigma
  pval1 <- EnrichTest(t_val, solve_eq, samp_rho, go_term)$p  # use the sample correlation matrix 
  solve_eq2 <- solve.equations(t_val, r1)                     # estimate the xi, beta and sigma
  pval2 <- EnrichTest(t_val, solve_eq2, r1, go_term)$p        # use the true correlation matrix
  pval3 <- summary(lm(t_val ~ go_term))$coe[2, 4]             # simple linear regression
  
return(c(pval1, pval2,pval3, unlist(solve_eq), unlist(solve_eq2)))
 }
)

colnames(fti)[1:3] <- c("p_val_est_cor", "p_val_tru_cor", "p_val_lm")

write.table(fti, "/home/stats/zhuob/data/computing/Enrichment_testStatistic0908_rep3.txt")

}




###############  ###  The BLOCK DIAGNAL CASE  

if (INDICATOR == "B")
{

fti <- foreach(i = 1:1000, .combine = rbind) %dopar% {
  library(MASS)
 size <- 500                                      # number of samples to be considered
  n_gene <- 500                                  # number of genes to be considered
  
  go_term_gene <- 0.2*n_gene                     # the go term 
  go_term <- rep(0, n_gene)
  go_term[1:go_term_gene] <- 1
  
  rho <- 0.5                                     # correlation matrix
  sigma0 <- 2                                    # variance part 
  cor.struct <- matrix(rho, go_term_gene, go_term_gene);    
  diag(cor.struct) <- 1                          # correlation matrix
  sigma1 <- sigma0*cor.struct                     # covariance matrix for both groups
  block_1 <- cor(mvrnorm(n=1000, rep(0, go_term_gene), sigma1))
  
  rho <- 0.2                                     # correlation matrix
  sigma0 <- 2                                    # variance part 
  cor.struct <- matrix(rho, n_gene- go_term_gene, n_gene - go_term_gene);    
  diag(cor.struct) <- 1                          # correlation matrix
  sigma2 <- sigma0*cor.struct                     # covariance matrix for both groups
  
  block_2 <- cor(mvrnorm(n=1000, rep(0, n_gene- go_term_gene), sigma2))
  
  cor_total <- matrix(0, n_gene, n_gene)
  cor_total[1:go_term_gene, 1:go_term_gene] <- block_1
  cor_total[(go_term_gene +1):n_gene, (go_term_gene +1):n_gene] <- block_2
  
  diag(cor_total) <- 1
  
  r1 <- cor_total
#  r1 <- diag(n_gene)
  xi <- xi_pre[i]
  sigma <- sigma0*( (1-xi)*diag(1, n_gene) + xi*r1)                     # covariance matrix for both groups


  delta <- rep(0, n_gene)             ## NULL, NO TREATMENT effect
  mu_initial <- 0
  mu1 <- rep(mu_initial, n_gene)                 # mean expression value for control
  mu2 <- delta  + mu_initial                     # mean expression value for treatment
  
  
  microarray <- matrix(NA, nrow=n_gene, ncol=size)   # expression matrix
  n <- size/2
  microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma)) 
  microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))
  
  trt <- rep(c(0, 1), each= size/2)             # group indicator
  
  
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
  
  samp_rho <- cor(t(resid_mat))                            # sample correlation matrix
  t_val <- group_mean[, size] - group_mean[, 1]            # control - treatment
    
  solve_eq <- solve.equations(t_val, samp_rho)                   # estimate the xi, beta and sigma 
  pval1 <-  EnrichTest(t_val,solve_eq, samp_rho, go_term)$p
  solve_eq2 <- solve.equations(t_val, r1)                  # estimate the xi, beta and sigma
  pval2 <-  EnrichTest(t_val,solve_eq2, r1, go_term)$p
  pval3 <- summary(lm(t_val~ go_term))$coe[2, 4]
 
 return(c(pval1, pval2, pval3, unlist(solve_eq), unlist(solve_eq2))) 
 }

 colnames(fti)[1:3] <- c("p_val_est_cor", "p_val_tru_cor", "p_val_lm")

write.table(fti, "/home/stats/zhuob/data/computing/Enrichment_testStatistic0908_rep4.txt")

}




#################################### Mean of each block ########################################

if (INDICATOR == "C")
{
  fti <- foreach(i = 1:1000, .combine = rbind) %dopar% {

	library(MASS)
	size <- 50                                               # number of samples to be considered
	n_gene <- 500                                            # number of genes to be considered
	go_term <- rbinom(n_gene, 1, 0.2)                        # random sample of go term
	r1 <- generate.cor(2, 0.5)                               # Covariance structure

	xi <- xi_pre[i]                                                # true xi
	sigma0 <- 2    
	sigma <- sigma0*( (1-xi)*diag(1, n_gene) + xi*r1)         # covariance matrix for both groups

############# generate expression matrix  ------------------------

	mu_initial <- 0
	mu1 <- rep(mu_initial, n_gene)                            # mean expression value for control
	mu2 <- mu1
	microarray <- matrix(nrow=n_gene, ncol=size)              # expression matrix
	n <- size/2
	microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma))
	microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))
	trt <- rep(c(0, 1), each= size/2)                         # group indicator

############### prepare the t_stat, samp_rho and Go term --------------

#########  t_test & sample correlation
	group_mean <- as.matrix(group.mean(microarray, trt))      # calculate correlation matrix
	resid_mat <- microarray - group_mean                      # the trt effects are removed from matrix
	samp_rho <- cor(t(resid_mat))                             # sample correlation matrix
	t_val <- group_mean[, size] - group_mean[, 1]             # control - treatment


        r2 <- matrix(NA, n_gene, n_gene)
        r2[go_term==1, go_term==1] <- mean(samp_rho[go_term==1, go_term==1]) 	# use the mean to represent all	the blocks elements
        r2[go_term==1, go_term==0] <- mean(samp_rho[go_term==1, go_term==0])
        r2[go_term==0, go_term==1] <- mean(samp_rho[go_term==0, go_term==1])
        r2[go_term==0, go_term==0] <- mean(samp_rho[go_term==0, go_term==0])
        diag(r2) <- 1


######## begin the test procedure ---------------------------------------
	solve_eq <- solve.equations(t_val, r1)                    # estimate the xi, beta and sigma
	solve_eq_mean <- solve.equations(t_val, r2)               # estimate the xi, beta and sigma

  return( c(unlist(solve_eq), unlist(solve_eq_mean) ) )
  }

write.table(fti, "/home/stats/zhuob/data/computing/Enrichment_testStatistic0908_mean.txt")

}



####################################################################################################
#  ENRICHMENT ANALYSIS BY SIMULATING DE GENES
##-------------------------------------------------------------------------------------------

if ( INDICATOR == "D") {

fti <- foreach(i = 1:1000, .combine = rbind) %dopar% {

################  Set up the parameters  -----------------
  library(MASS)
  size <- 6                                     # number of samples to be considered
  n_gene <- 500                                  # number of genes to be considered
  de_prop <- 0.2                                 # proportion of DE
  de <- n_gene*de_prop                           # number of DE genes
  delta <- rnorm(n_gene, 2, 1)                     # magnitude of DE, mean(DE) = 1/0.5 = 2
  delta <- rep(1,n_gene)
  z <- sample(1:n_gene, de_prop*n_gene)             # indicator, DE =1, NON_DE=0
  delta[-z] <- 0

#### therefore about 100 genes are DE ############
#### Among them 20 are in the GO term, 20/100 = 20% are in GO term,  (100-20)/(500-100) = 20% DE are not in GO term #### 
#### the result is: No enrichment 

######### prepare go term ---------------------------------
  n_go <- 0.2*n_gene                                       # number of genes in the go term
  go_term <- rep(0, n_gene)                                # go term indicator
  go_prop <- 0.2                                           # proportion of DE genes in the go term
  go_de <- sample(which(delta !=0), n_go*go_prop)               # de genes in go term
  go_not_de <- sample(which(delta ==0), n_go*(1-go_prop))	   # non DE genes in go term
  go_in <- c(go_de, go_not_de)
#   go_in <- sample(1:n_gene, n_go)
  go_term[go_in] <- 1

################# To generate the correlation matrix --------------

  set.seed(100)                   # fixed correlation matrix
#   zm <- rnorm(n_gene*size)
#   dim(zm) <- c(n_gene, size)
#   sigma <- cor(t(zm))

  ids <- which(go_term ==1)
  zm <- rnorm(length(ids)*size)
  dim(zm) <- c(length(ids), size)
  sigma.set <- cor((t(zm)))
  sigma <- diag(n_gene)
  sigma[ids, ids] <- sigma.set
#  sigma[ids, ids] <- 0.1
#  diag(sigma) <- 1

# rho <- 0.3
# sigma <- matrix(rho, n_gene, n_gene)
# diag(sigma) <- 1


#    sigma <- diag(n_gene)                    # independent case


set.seed(NULL)                               # nullify the seed to generate random numbers


############# generate expression matrix  ------------------------

  mu_initial <- 0
  mu1 <- rep(mu_initial, n_gene)                 # mean expression value for control
  mu2 <- delta  + mu_initial                     # mean expression value for treatment

  microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix
  n <- size/2
  microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma))
  microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))
  trt <- rep(c(0, 1), each= size/2)              # group indicator

############### prepare the t_stat, samp_rho and Go term --------------

#########  t_test & sample correlation
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix

  samp_rho <- cor(t(resid_mat))                            # sample correlation matrix
  t_val <- group_mean[, size] - group_mean[, 1]            # control - treatment

######## begin the test procedure ---------------------------------------
  solve_eq <- solve.equations(t_val, samp_rho)             # estimate the xi, beta and sigma
  pval1 <- EnrichTest(t_val, solve_eq, samp_rho, go_term)$p  # use the sample correlation matrix
  pval2 <- summary(lm(t_val ~ go_term))$coe[2, 4]             # simple linear regression




### limma package
    ##### LAST MODIFIED OCT 16, 2015
        library(limma)
        design <- model.matrix(~trt)

        fit <- lmFit(microarray, design)
        fit <- eBayes(fit)                                                                    # Emperical Bayes t test
        stat <- fit$t[, 2]                                                                    #
        alter <- "either"                                                                     # undirectional test
        index1 <- which(go_term==1)                                                           # which genes are in GOTERM

# up: positive statistics
# down: negative stats
# either: the set is either up or down-regulated as a whole
# mixed: tests whether the genes in the set tend to be DE without regard for direction.
#        In this case, the test will be significant if the set contains mostly large statistics, negative or positive

        tes1 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= F)         # moderated t GeneSet Rank
        tes2 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= T)         # moderated t GeneSet
        tes3 <- camera(microarray, index1,  design)$PValue                                    # camera proedure
        tes4 <- camera(microarray, index1,  design, use.ranks= T)$PValue                      # rankded version Camera

# pval1 <- EnrichmentTest(microarray, trt, go_term)$p          ## our method

return(c(pval1, pval2, tes1, tes2, tes3, tes4,  unlist(solve_eq)))
 }

colnames(fti)[1:6] <- c("OurTest", "lm",  "ModeratedT", "Rank", "Camera", "CameraRank")

write.table(fti, "/home/stats/zhuob/data/computing/Enrichment_testStatistic1016_temp.txt")


}


############################################# ENRICHMENT TEST WITH FISHER'S EXACT TEST OR LIMMA PACKAGE --------------------------------


if ( INDICATOR == "E") {

fti <- foreach(i = 1:1000, .combine = rbind) %dopar% {

################  Set up the parameters  -----------------
  library(MASS)
  size <- 6                                     # number of samples to be considered
  n_gene <- 500                                  # number of genes to be considered

## LAST MODIFIED on OCT 25, 2015, power simulation
   go_prop <- 0.2                                           # portion of DE
   go_term <- sample(1:n_gene, 100)                         # sample the go genes
   go_de <-go_term[ which(rbinom(100, 1, go_prop)==1)]      # randomly assign DE to genes
 
   no_go_prop = 0.2  
#   no_go_prop <- 0.2                                        # proportion of DE in no_go
    no_go_term <- (1:n_gene)[-go_term]                       # the rest are non go term genes
   no_go_de <- no_go_term[which(rbinom(400, 1, no_go_prop)==1)]; # which genes are DE in non-go term
 
   delta <- rnorm(n_gene, 2, 1)                             # magnitiude of DE
   delta[-c(go_de, no_go_de)] <- 0                          # non de genes have 0 
#   delta <- rep(0, n_gene)                                  # SIMULATE THE CASE IN CAMERA PAPER
  
   go_term2 <- rep(0, n_gene)  
   go_term2[go_term] <- 1

   go_term <- go_term2                                      # prepare the GO term

################# To generate the correlation matrix --------------
  set.seed(100)                   # fixed correlation matrix
 #  zm <- rnorm(n_gene*size)
 #  dim(zm) <- c(n_gene, size)
 #  sigma <- cor(t(zm))

 
 sigma <- diag(n_gene)
 ids <- which(go_term ==1)
#  zm <- rnorm(length(ids)*size)
#  dim(zm) <- c(length(ids), size)
#  sigma.set <- cor((t(zm)))
#  sigma[ids, ids] <- sigma.set
  sigma[ids, ids] <- 0.1

#  sigma[-ids, -ids] <- 0.1
#  sigma[ids, -ids] <- -0.1   ## negative correlation between go term and non-go term 
#  sigma[-ids, ids] <- -0.1   ##
  diag(sigma) <- 1

# rho <- 0.1
# sigma <- matrix(rho, n_gene, n_gene)
# diag(sigma) <- 1


#    sigma <- diag(n_gene)                    # independent case


set.seed(NULL)                               # nullify the seed to generate random numbers

############# generate expression matrix  ------------------------

  mu_initial <- 0
  mu1 <- rep(mu_initial, n_gene)                 # mean expression value for control
  mu2 <- delta  + mu_initial                     # mean expression value for treatment

  microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix
  n <- size/2
  microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma))
  microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))
  trt <- rep(c(0, 1), each= size/2)              # group indicator

############### prepare the t_stat, samp_rho and Go term --------------

  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix

  samp_rho <- cor(t(resid_mat))                            # sample correlation matrix
  t_val <- group_mean[, size] - group_mean[, 1]            # control - treatment


######## begin the test procedure ---------------------------------------
  solve_eq <- solve.equations(t_val, samp_rho)             # estimate the xi, beta and sigma
  pval1 <- EnrichTest(t_val, solve_eq, samp_rho, go_term)$p  # use the sample correlation matrix
  pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]            # linear regression

### limma package
    ##### LAST MODIFIED OCT 16, 2015
	library(limma)
	design <- model.matrix(~trt) 

	fit <- lmFit(microarray, design)
	fit <- eBayes(fit)								      # Emperical Bayes t test	
	stat <- fit$t[, 2]            							      #	
	alter <- "either"                                              			      # undirectional test	 
	index1 <- which(go_term==1)                                                           # which genes are in GOTERM

# up: positive statistics
# down: negative stats
# either: the set is either up or down-regulated as a whole
# mixed: tests whether the genes in the set tend to be DE without regard for direction.
#        In this case, the test will be significant if the set contains mostly large statistics, negative or positive

	tes1 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= F)         # moderated t GeneSet Rank
	tes2 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= T)         # moderated t GeneSet
	tes3 <- camera(microarray, index1,  design)$PValue                                    # camera proedure
	tes4 <- camera(microarray, index1,  design, use.ranks= T)$PValue                      # rankded version Camera

return(c(pval1, pval2, tes1, tes2, tes3, tes4,  unlist(solve_eq)))
 }

colnames(fti)[1:6] <- c("OurTest", "lm",  "ModeratedT", "Rank", "Camera", "CameraRank")

write.table(fti, "/home/stats/zhuob/data/computing/Enrichment_testStatistic1025_temp.txt")

}





if (INDICATOR == "F"){

fti <- foreach(i = 1:1000, .combine = rbind) %dopar% {

library(MASS)
library(limma)
  m <- 500; n <- 50;
#  zm <- rnorm(m*n)
#  dim(zm) <- c(m, n)
# sigma <- 1.5*cor(t(zm))

  rho <- 0.5
  s1 <- diag(m)
  sigma <-1.5* rho^abs(row(s1)-col(s1))
 
  mu <- rep(0, m)
   
 
  y <- t(mvrnorm(n, mu, Sigma = sigma))
  # y <- matrix(rnorm(100*4),100,4)
  design <- cbind(Intercept=1,Group=rep(c(0, 1), each=n/2))
  
  # First set of 5 genes contains 3 that are genuinely differentially expressed
  index1 <- 1:(0.1*m)
  delta <- 1
  y[index1,(n/2+1):n] <- y[index1,(n/2+1):n]+ delta
  
  # Second set of 5 genes contains none that are DE
 # index2 <- 6:10

	index1 <- c(index1, 201:250)
	index0 <- rep(0, m)
	index0[index1] <- 1
  
  fit <- lmFit(y, design )
  fit <- eBayes(fit)
  stat <- fit$t[, 2]
  alter <- "up"
  
  p1 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= F)
  p2 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= T)
  p3 <- camera(y, index1,  design, contrast=2, use.ranks =F )$PValue
  p4 <- camera(y, index1,  design, contrast=2, use.ranks =T )$PValue

  trt <- rep(c(0,1), each=n/2)
  p5 <- EnrichmentTest(y, trt, index0)$p

return(c(p1, p2, p3, p4, p5))

}


colnames(fti)[1:5] <- c( "ModeratedT", "Rank", "Camera", "CameraRank", "OurTest")

write.table(fti, "/home/stats/zhuob/data/computing/Enrichment_testStatistic1013_4.txt")

}

