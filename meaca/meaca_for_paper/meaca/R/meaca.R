



#' meaca for single gene set test.
#'
#' @title meaca-single. 
#' @param expression_data  the expressoin matrix.
#' @param trt  treatment labels.
#' @param go_term  an indicator vector. 1 for genes in the test, 0 otherwise.
#' @param standardize  whether the data should be standaridzed.
#' @return a list 
#' \item{set_name}{The name of the gene set}
#' \item{testSetCor}{Average correlation for genes in the test set}
#' \item{interCor}{Average correlation between genes in the test set and those not in the test set}
#' \item{backSetCor}{Average correlations for genes not in the test set.}
#' @export
#' @examples



#####  our test
meaca_single <- function(expression_data, trt, go_term, standardize=F){
  ##  we need to standardize it here.
  
  if (standardize == T){             # do the standardization 
    expression_data <- standardize_expression_data(expression_data, trt)
  }
  
  est_sigma <- estimate_sigma(expression_data, trt)
  sigma <- est_sigma$sigma
  t_val <- est_sigma$t_val
  
  beta0 <- mean(t_val)
  ones <- rep(1, length(t_val))
  
  numer <- (crossprod(go_term, (t_val - beta0*ones)))^2
  go_bar <- rep(mean(go_term), length(go_term))
  denom2 <- t(go_term-go_bar) %*% sigma %*% (go_term- go_bar)
  
  test_stat <- drop(numer/denom2)
  pval <- 1 - pchisq(test_stat, 1)	
  
  ###  March 16:  one sided test modification
  t_temp <- mean(t_val[go_term == 1]) - mean(t_val[go_term ==0])
  
  ## decide whether the set is up- or down-regulated as a whole
  if (t_temp > 0 ){
    status <- "up"
  } else {
    status  <- "down"
  }  
  
  test_stat2 <- sqrt(test_stat)	
  pval2 <- 2*min(1- pnorm(test_stat2), test_stat2)  # two sided p-value
  
  return(list(stat = test_stat, p1 = pval, status = status, p2 = pval2))
}







#' meaca for testing multiple gene sets.
#'
#' @title meaca-multiple. 
#' @param expression_data  the expressoin matrix.
#' @param trt  treatment labels.
#' @param geneset  gene sets to be tested, an object from \code{read_gene_set}.
#' @param standardize  whether the data should be standaridzed.
#' @param minSetSize  the minimum number of genes contained for a gene set to be considered.
#' @param fdr_method  which method is ued to adjust the p values. see arguments in function \code{p.adjust}.
#' @return a data frame
#' @export
#' @examples





meaca_multiple <- function(expression_data, trt, geneset, standardize = T, minSetSize = 5, fdr_method = "BH"){
  ## conduct enrichment analysis for a battery of gene sets and return the raw p values 
  ## as well as adjusted p values
  ##  input: 
  ##      microarray: the raw expression matrix
  ##             trt: the treatment lables
  ##         geneset: object from readGeneSet function
  ##          
  
  ############### for GSE64810 USE the following two lines. Otherwise comment it
  # for GSE64810
  #	expression_data<- expression_data[, -(1:2)]
  #	all_genes <-  expression[, 2]
  
  # for other typical data, where row names are genes	
  all_genes <- rownames(expression_data)	
  
  if (standardize == T){             # do the standardization 
    expression_data <- standardize_expression_data(expression_data, trt)
  }
  
  est_sigma <- estimate_sigma(expression_data, trt)
  sigma <- est_sigma$sigma
  t_val <- est_sigma$t_val
  
  
  beta0 <- mean(t_val)
  ones <- rep(1, length(t_val))
  
  keep_term <- which(geneset$size >= minSetSize)
  print(paste("number of sets to be tested:", length(keep_term)))
  
  #    	pval <- c()
  pval <- pval2 <- status <- c() 
  set_size <- c()      
  set_name <- c()
  
  cat("\n this might take a long time...  \n")    
  cat("\n testing gene set: \n")
  
  for ( i in 1:length(keep_term)){
    #  for ( i in 1:50){
    gset1 <- geneset$gene_set[[keep_term[i]]][-(1:2)]
    cat('\r', i)      
    set_name[i] <- geneset$gene_set[[keep_term[i]]][1]
    go_term <- ifelse(all_genes %in% gset1, 1, 0)
    set_size[i] <- sum(go_term) 
    
    numer <- (crossprod(go_term, (t_val - beta0*ones)))^2
    go_bar <- rep(mean(go_term), length(go_term))
    denom2 <- t(go_term-go_bar) %*% sigma %*% (go_term- go_bar)
    
    test_stat <- drop(numer/denom2)
    pval[i] <- 1 - pchisq(test_stat, 1)
    
    
    ###  March 16:  one sided test modification
    t_temp <- mean(t_val[go_term == 1]) - mean(t_val[go_term ==0])
    if(!is.finite(t_temp)) {    # it's possible that there's no single gene matched between the data and the test set
      pval2[i] <- NA
      status[i] <- "Unknown"
    }  else {
      
      ## decide whether the set is up- or down-regulated as a whole
      if (t_temp > 0 ){ status[i] <- "up"
      } else { status[i]  <- "down"}  
      
      test_stat2 <- sqrt(test_stat)	
      pval2[i] <- 2*min(1- pnorm(test_stat2), test_stat2)  # two sided p-value
    }
    
    #  print(c (set_name[i], set_size[i], status[i], pval[i], pval2[i]))
  }  
  
  cat("\n")
  p1_fdr <- p.adjust(pval, method = fdr_method)
  
  p2_fdr <- p.adjust(pval, method = fdr_method)
  
  results <- data.frame(set_name = set_name, set_size = set_size, status = status,
                        p1=pval,p1_fdr = p1_fdr, p2 = pval2, p2_fdr = p2_fdr)
  
  ## calculate the between gene correlations
  btw <-btw_gene_corr(expression_data, trt, geneset, minSetSize=2, standardize = T)   
  result_comb <- merge(results, btw, by = "set_name")	
  
  return(result_comb)
}





#' Standardize the expression data.
#'
#' @title standardize expression data, with method described in the paper. 
#' @param expression_data  the expression matrix. 
#' @param trt  sample labels. 0 for control and 1 for treatment
#' @return a matrix of the same dimension as input data.
#' @export


standardize_expression_data <- function(expression_data, trt){
  ## standardize the expression data, to make it have variance 1 in each of treatmeent/control group
  #  first remove the means within each group to get sd, and then original data by sd
  
  group1 <- trt == 1
  group2 <- trt == 0
  col_index <- 1:ncol(expression_data)
  
  set1 <- expression_data[, group1]
  set2 <- expression_data[, group2]
  
  s1 <- apply(set1, 1, sd)
  s2 <- apply(set2, 1, sd)
  
  ### modified on Feb 18, 2016, use pooled variacne to do standardization....
  
  n1 <- length(group1)
  n2 <- length(group2)
  
  poolsd <-  sqrt( ( s1^2*(n1-1) + s2^2*(n2-1) ) / (n1 + n2-2) )  # calculate the pooled standard deviation  
  poolsd[poolsd ==0] <- 1           ### if the sd is 0, replace with 1. 
  
  new_data <- data.frame(matrix(NA, nrow(expression_data), ncol(expression_data)))
  
  new_data <- expression_data/poolsd
  
  rownames(new_data) <- rownames(expression_data); 
  colnames(new_data) <- colnames(expression_data)
  
  return(new_data)
}





trt_mean <- function(data, trt){
  # a function to calculate treatment means
  # return a matrix
  id1 <- which(trt == 1)
  nrep1 <- length(id1)
  nrep2 <- ncol(data)-nrep1
  mean1 <- apply(data[, id1], 1, mean)
  mean1 <- matrix(rep(mean1, nrep1), ncol=nrep1, byrow=F) 
  mean2 <- apply(data[, -id1], 1, mean)
  mean2 <- matrix(rep(mean2, nrep2), ncol=nrep2, byrow=F) 
  
  data_mean <- data.frame(matrix(NA, nrow(data), ncol(data)))
  data_mean[, id1] <- mean1
  data_mean[, -id1] <- mean2
  colnames(data_mean) <- colnames(data)
  rownames(data_mean) <- rownames(data)
  data_mean
}





#' Estimate sample covariance and calculate the gene-level statistics
#'
#' @title Estimate sample covariance. 
#' @param expression_data  the expression matrix. 
#' @param trt  sample labels. 0 for control and 1 for treatment
#' @return a list 
#' \item{sigma}{a covariance matrix }
#' \item{t_val}{a vector of gene level test statistics}
#' @export
# #' @examples


estimate_sigma <- function(expression_data, trt){
  
  expression_data <- as.matrix(expression_data)
  group_mean <- as.matrix(trt_mean(expression_data, trt))     # calculate correlation matrix
  resid_mat <- expression_data - group_mean                     # the trt effects are removed from matrix
  # samp_rho <- cov(t(resid_mat))
  samp_rho <- cor(t(resid_mat))      # because of the way we do standardization           
  
  # two lines correction for inflated type I error.	  
  #   nsamp <- ncol(expression_data)
  #  samp_rho <- cov(t(resid_mat)) * (nsamp-1)/(nsamp-2)
  
  factors <- unique(trt)
  m1 <- sum(trt == factors[1])
  m2 <- sum(trt == factors[2]) 
  sigma2 <- 1/m1 + 1/m2
  
  ndim <- nrow(expression_data)
  
  size <- length(trt)
  idex1 <- which(trt == factors[1])[1]
  idex2 <- which(trt == factors[2])[1]
  t_val <- group_mean[, idex2] - group_mean[, idex1]            # control - treatment
  
  sigma0 <- var(t_val) - sigma2
  sigma0 <- max(0, sigma0)                                     # in case of negative sigma0
  sigma <-  sigma0*diag(ndim) + sigma2*samp_rho
  
  
  return(list(sigma = sigma, t_val = t_val))
  
}



#' Average correlations for genes
#'
#' @title Estimate sample correlation. 
#' @param expression_data  the expressoin matrix.
#' @param trt  treatment labels
#' @param geneset  an object from \code{read_gene_set()}
#' @param standardize  whether the data should be standaridzed
#' @param minSetSize the minimum number of genes contained for a gene set to be considered.
#' @return a list 
#' \item{set_name}{The name of the gene set}
#' \item{testSetCor}{Average correlation for genes in the test set}
#' \item{interCor}{Average correlation between genes in the test set and those not in the test set}
#' \item{backSetCor}{Average correlations for genes not in the test set.}
#' @export
#' @examples


btw_gene_corr <- function(expression_data, trt, geneset, standardize = T, minSetSize = 5){
  ## calculate the mean correlation for testset genes, the background set genes and the inter-set correlations
  
  # for GSE64810
  # 	microarray <- expression_data[, -(1:2)]
  #  	all_genes <-  expression_data[, 2]
  
  # for other typical data, where row names are genes
  all_genes <- rownames(expression_data)
  
  if (standardize == T){             # do the standardization
    expression_data <- standardize_expression_data(expression_data, trt)
  }
  
  #	## calcuate the sample correlations
  expression_data <- as.matrix(expression_data)
  group_mean <- as.matrix(trt_mean(expression_data, trt))     # calculate correlation matrix
  resid_mat <- expression_data - group_mean                     # the trt effects are removed from matrix

#  samp_rho <- cor(t(resid_mat))
# modified on March 18, 2016
  samp_rho <- cov(t(resid_mat))
  
  keep_term <- which(geneset$size >= minSetSize)
  
  
  set_name <- set_size <- c()
  testSet_cor <- backSet_cor <- interCor <- sumTestCor <- rep(NA, length(keep_term))
  
  for ( i in 1:length(keep_term)){
    gset1 <- geneset$gene_set[[keep_term[i]]][-(1:2)]
    set_name[i] <- geneset$gene_set[[keep_term[i]]][1]
    go_term <- ifelse(all_genes %in% gset1, 1, 0)
    set_size[i] <- sum(go_term)
    
    go_ind <- which(go_term==1)    # which are go_terms
    sumTestCor[i] <- sum(samp_rho[go_ind, go_ind])
    testSet_cor[i] <- (sumTestCor[i] - set_size[i])/( set_size[i]*(set_size[i]-1) )
    interCor[i] <- mean(samp_rho[-go_ind, go_ind])
    
    # print(c(set_size[i], sumTestCor[i], interCor[i], testSet_cor[i])
  }
  
  # the correlations for the background genes
  sumAll <- sum(samp_rho); ndim <- nrow(samp_rho)
  n2 <- ndim - set_size
  backSet_cor <- ( sumAll - (interCor*n2*2 + sumTestCor) - n2 ) / (n2^2 - n2)
  
  cor_matrix <- data.frame(set_name = set_name, testSetCor = testSet_cor, interCor= interCor, backSetCor = backSet_cor)
  
  return(cor_matrix)
}








