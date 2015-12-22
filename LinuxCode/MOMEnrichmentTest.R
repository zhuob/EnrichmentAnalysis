## use method of moment to estimate var(delta) and then conduct the
## enrichment analysis...

estimate.sigma <- function(microarray, trt){
  
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
  samp_rho <- cor(t(resid_mat))                 
  
  factors <- unique(trt)
  m1 <- sum(trt == factors[1])
  m2 <- sum(trt == factors[2]) 
  sigma2 <- 1/m1 + 1/m2
  
#  var_diag <- apply(resid_mat, 1, var)
#  sigma0 <- var_diag - sigma2
#  sigma0[sigma0 < 0] <- 0
  
  # what if I don't use samp_rho but some random correlation matrix
	
#  n_gene <- nrow(microarray); size <- length(trt)
#  zm <- rnorm(n_gene*size)
#  dim(zm) <- c(n_gene, size)
#  samp_rho <- cor(t(zm))

  # What if I use the true correlation?
#  samp_rho <- dat$sigma                 				# this has been passed by compare_test()


  ndim <- nrow(microarray)
#  sigma <-  sigma0*diag(ndim) + sigma2*samp_rho

  size <- length(trt)
  idex1 <- which(trt == factors[1])[1]
  idex2 <- which(trt == factors[2])[1]
  t_val <- group_mean[, idex2] - group_mean[, idex1]            # control - treatment

  sigma0 <- var(t_val) - sigma2
  sigma <-  sigma0*diag(ndim) + sigma2*samp_rho
  
  
  return(list(sigma = sigma, t_val = t_val))

}
  

#####  our test
MOM_test <- function(microarray, trt, go_term){
  
   est_sigma <- estimate.sigma(microarray, trt)
   sigma <- est_sigma$sigma
   t_val <- est_sigma$t_val

  
   beta0 <- mean(t_val)
   ones <- rep(1, length(t_val))
   
   ## calculate the test statistics
   inv_sigma <- solve(sigma)
   numer <- (crossprod(go_term, (t_val - beta0*ones)))^2
#   denom <- t(go_term) %*% sigma %*% (go_term) - (crossprod(go_term, ones))^2/(t(ones) %*% inv_sigma %*% ones)
   
   go_bar <- rep(mean(go_term), length(go_term))
   denom2 <- t(go_term-go_bar) %*% sigma %*% (go_term- go_bar)
   #print(c(denom, denom2))
   
   test_stat <- drop(numer/denom2)
   pval <- 1 - pchisq(test_stat, 1)
   
   return(list(stat = test_stat, p = pval))
}
  

compare_test <- function(dat){
  
  microarray <- dat$data
  trt <- dat$trt
  go_term <- dat$go_term
  sigma <- dat$sigma
  
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
  
  samp_rho <- cor(t(resid_mat))                            # sample correlation matrix
  t_val <- group_mean[, size] - group_mean[, 1]            # control - treatment
  
  
  #### our test
  pval1 <- MOM_test(microarray, trt, go_term)$p
  pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]            # linear regression
  
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
  tes4 <- camera(microarray, index1,  design, use.ranks= T)$PValue       
  
  return(c(pval1, pval2, tes1, tes2, tes3, tes4))
  
}

