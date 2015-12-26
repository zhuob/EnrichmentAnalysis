## use method of moment to estimate var(delta) and then conduct the
## enrichment analysis...

library(limma)         # needed for camera() and geneSetTest()


group.mean <- function(data, group){
  # a function to calculate group mean
  
  id1 <- which(group == 1)
  nrep1 <- length(id1)
  nrep2 <- ncol(data)-nrep1
  mean1 <- apply(data[, id1], 1, mean)
  mean1 <- matrix(rep(mean1, nrep1), ncol=nrep1, byrow=F) 
  mean2 <- apply(data[, -id1], 1, mean)
  mean2 <- matrix(rep(mean2, nrep2), ncol=nrep2, byrow=F) 
  
  data.mean <- data.frame(matrix(NA, nrow(data), ncol(data)))
  data.mean[, id1] <- mean1
  data.mean[, -id1] <- mean2
  colnames(data.mean) <- colnames(data)
  rownames(data.mean) <- rownames(data)
  data.mean
}



standardize.microarray <- function(data, trt){
  ## standardize the expression data, to make it have variance 1 in each of treatmeent/control group
  #  first remove the means within each group to get sd, and then original data by sd
  data1 <- as.matrix(data)
  
  std_by_row <- function(row_dat, trt){
    
    group1 <- trt == 1
    group2 <- trt == 0
    
    x1 <- row_dat[group1]
    x2 <- row_dat[group2]
    
    x1_center <- x1-mean(x1)
    vas1 <- sd(x1_center)
    x1_new <- x1/vas1
    
    x2_center <- x2-mean(x2)
    vas2 <- sd(x2_center)
    x2_new <- x2/vas2
    
    x_return <- rep(NA, length(row_dat))  
    x_return[group1] <- x1_new
    x_return[group2] <- x2_new
    return(x_return)
  }
  new_data <- data.frame(matrix(NA, nrow(data), ncol(data)))
  
  for ( i in 1:nrow(data1)){
    new_data[i, ] <- std_by_row(data1[i, ], trt)
  }
  
  rownames(new_data) <- rownames(data)
  colnames(new_data) <- colnames(data)
  
  return(new_data)
  
}




estimate.sigma <- function(microarray, trt){
  
  microarray <- as.matrix(microarray)
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
  samp_rho <- cor(t(resid_mat))                 
  
  factors <- unique(trt)
  m1 <- sum(trt == factors[1])
  m2 <- sum(trt == factors[2]) 
  sigma2 <- 1/m1 + 1/m2

  ndim <- nrow(microarray)
  
  size <- length(trt)
  idex1 <- which(trt == factors[1])[1]
  idex2 <- which(trt == factors[2])[1]
  t_val <- group_mean[, idex2] - group_mean[, idex1]            # control - treatment

  sigma0 <- var(t_val) - sigma2
  sigma0 <- max(0, sigma0)                                     # in case of negative sigma0
  sigma <-  sigma0*diag(ndim) + sigma2*samp_rho
  
  
  return(list(sigma = sigma, t_val = t_val))

}
  

#####  our test
MOM_test <- function(microarray, trt, go_term, standardize=T){
##  we need to standardize it here.

  if (standardize == T){             # do the standardization 
    microarray <- standardize.microarray(microarray, trt)
  }

   est_sigma <- estimate.sigma(microarray, trt)
   sigma <- est_sigma$sigma
   t_val <- est_sigma$t_val

   beta0 <- mean(t_val)
   ones <- rep(1, length(t_val))
   
   numer <- (crossprod(go_term, (t_val - beta0*ones)))^2
   go_bar <- rep(mean(go_term), length(go_term))
   denom2 <- t(go_term-go_bar) %*% sigma %*% (go_term- go_bar)

   test_stat <- drop(numer/denom2)
   pval <- 1 - pchisq(test_stat, 1)
   
   return(list(stat = test_stat, p = pval))
}
  






compare_test <- function(dat){
########  a function to incorporate all test procedures

    microarray <- dat$data
  trt <- dat$trt
  go_term <- dat$go_term
 
  pval1 <- MOM_test(microarray, trt, go_term)$p             # our test
  
  est_sigma <- estimate.sigma(microarray, trt)
  t_val <- est_sigma$t_val
  pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]            # linear regression
  

  design <- model.matrix(~trt) 
  fit <- lmFit(microarray, design)
  fit <- eBayes(fit)								                # Emperical Bayes t test	
  stat <- fit$t[, 2]            							      #	use the moderated t statistics to do enrichment
  alter <- "mixed"                                  # This is the default option, see below for explanation.
  index1 <- which(go_term==1)                       # which genes are in GOTERM
  
  # up: positive statistics
  # down: negative stats
  # either: the set is either up or down-regulated as a whole
  # mixed: tests whether the genes in the set tend to be DE without regard for direction.
  #        In this case, the test will be significant if the set contains mostly large statistics, negative or positive
  
  tes1 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= F)         # moderated t GeneSet
  tes2 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= T)         # moderated t GeneSet Rank
  tes3 <- camera(microarray, index1,  design)$PValue                                    # camera proedure
  tes4 <- camera(microarray, index1,  design, use.ranks= T)$PValue                      # camera rank 
  
  tes5 <- GSEA.SingleSet(dat$data, dat$trt, dat$go_term, nperm=1000)$p.vals
  
  
  return(c(pval1, pval2, tes1, tes2, tes3, tes4, tes5))
  
}

