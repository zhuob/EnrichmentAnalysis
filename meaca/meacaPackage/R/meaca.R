#' meaca for single gene set test.
#'
#' @title Testing the enrichment status of a single pre-defined gene set. 
#' @param expression_data  the expressoin matrix.
#' @param trt  treatment labels, should be 0s and 1s.
#' @param go_term  an indicator vector with value 1 for genes in the test, 0 otherwise.
#' @param standardize  whether the data should be standardized. Set to be TRUE for real data analysis.
#' @return a list 
#' \item{stat}{test statistic for meaca}
#' \item{p1}{chi square test p value}
#' \item{status}{"Up" for up-regulated gene sets and "Down" for down-regulated gene sets}
#' \item{p2}{a two-sided p value for conducting t-test from the test statstic \code{stat}}
#' @export
#' @examples  
#' m <- 100   # number of rows (genes)
#' n <- 50   # number of columns (samples)
#' y <- matrix(rnorm(5000), 100, 50)   # expression  matrix
#' trt <- rep(c(0, 1), each = n/2)      # treatment labels
#' go_term <- rep(0, m)
#' go_term[sample(1:m, 20)] <- 1     # the rows in the test set are labeled as 1
#' result <- meaca_single(y, trt, go_term)




#####  our test
meaca_single <- function(expression_data, trt, go_term, standardize=T){
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
    status <- "Up"
  } else {
    status  <- "Down"
  }  
  
  test_stat2 <- sqrt(test_stat)	
  pval2 <- 2*min(1- pnorm(test_stat2), test_stat2)  # two sided p-value
  
  return(list(stat = test_stat, p1 = pval, status = status, p2 = pval2))
}







#' meaca for testing multiple gene sets.
#'
#' @title meaca-multiple. 
#' @param expression_data  the expressoin matrix.
#' @param trt  treatment labels, 1 for treatment 0 for control.
#' @param geneset  gene sets to be tested, an object from \code{read_gene_set}.
#' @param standardize  whether the data should be standaridzed.
#' @param minSetSize  the minimum number of genes contained in a gene set for it to be considered as a test set.
#' @param fdr_method  which method is ued to adjust the p values. Options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", see arguments in function \code{p.adjust} for more details.
#' @return a data frame with the following columns
#' \item{set_name}{the name of the gene set }
#' \item{set_sizel}{the size of the gene set}
#' \item{status}{"Up" for up-regulated gene sets and "Down" for down-regulated gene sets} 
#' \item{p1}{chi square test p value}
#' \item{p1_fdr}{adjusted p1, using \code{p.adjust} with appropriate method}
#' \item{p2}{a two-sided t-test p value equivalent to the chi square test}
#' \item{p2_fdr}{adjusted p2 }
#' @export
#' @examples 
#' m <- 100   # number of rows (genes)
#' n <- 50    # number of columns (samples)
#' y <- matrix(rnorm(5000), 100, 50)              # expression  matrix
#' rownames(y) <- paste("gene", 1:m, sep = "")  
#' colnames(y) <- paste("sample", 1:n, sep = "")
#' trt <- rep(c(0, 1), each = n/2)                # treatment labels
#' ## create the gene set format, desired for this function.
#' total <- 15                     # number of gene sets in the database
#' set.seed(100)
#' size <- sample(1:30, replace=F, total)
#' go_term <- rep(0, m)
#' gs_line <- list()
#' set_name <- paste("set", 1:length(size), sep = "")
#' for (i in 1:total) {
#'     set.seed(i)
#'     gs_line[[i]] <- rownames(y)[sample(1:m, size[i])] 
#' }
#' gs <- list(total=total, size=size,set_name = set_name, gene_set=gs_line)
#' # run the multiple gene set test
#' result <- meaca_multiple(y, trt, gs)






meaca_multiple <- function(expression_data, trt, geneset, standardize = T, minSetSize = 5, fdr_method = "BH"){
  ## conduct enrichment analysis for a battery of gene sets and return the raw p values 
  ## as well as adjusted p values
  ##  input: 
  ##      microarray: the raw expression matrix
  ##             trt: the treatment lables
  ##         geneset: object from readGeneSet function
  ##          
  
 
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
  set_name <- geneset$set_name[keep_term]
  
  cat("\n this might take a long time...  \n")    
  cat("\n testing gene set: \n")
  
  
  for ( i in 1:length(keep_term)){
    #  for ( i in 1:50){
    gset1 <- geneset$gene_set[[keep_term[i]]]
    cat('\r', i)      
   
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
      if (t_temp > 0 ){ status[i] <- "Up"
      } else { status[i]  <- "Down"}  
      
      test_stat2 <- sqrt(test_stat)	
      pval2[i] <- 2*min(1- pnorm(test_stat2), test_stat2)  # two sided p-value
    }
    
    #  print(c (set_name[i], set_size[i], status[i], pval[i], pval2[i]))
  }  
  
  cat("\n")
  p1_fdr <- p.adjust(pval, method = fdr_method)
  
  p2_fdr <- p.adjust(pval, method = fdr_method)
  
  results <- data.frame(set_name = set_name, 
                        set_size = set_size, status = status,
                        p1=pval,p1_fdr = p1_fdr, p2 = pval2, p2_fdr = p2_fdr)
  
  ## calculate the between gene correlations
  btw <-btw_gene_corr(expression_data, trt, geneset, minSetSize=minSetSize, standardize = T)   
  result_comb <- merge(results, btw, by = "set_name")	
  
  return(result_comb)
}





#' Standardize the expression data.
#'
#' @title standardize expression data, with method described in the paper. 
#' @param expression_data  the expression matrix. 
#' @param trt  sample labels. 0 for control and 1 for treatment
#' @return the standardized expression, a matrix of the same dimension as input expression data.
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





#' Estimate test-statistic covariance
#'
#' @title Estimate test-statistic covariance. 
#' @param expression_data  the expression matrix. 
#' @param trt  sample labels. 0 for control and 1 for treatment
#' @return a list 
#' \item{sigma}{estimated covariance matrix for the gene-level test statistics. }
#' \item{t_val}{the gene level statistics, one for each gene}
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
#' @return a data frame with columns: 
#' \item{set_name}{The name of the gene set}
#' \item{testSetCor}{Average correlation for genes in the test set}
#' \item{interCor}{Average correlation between genes in the test set and those not in the test set}
#' \item{backSetCor}{Average correlations for genes not in the test set.}
#' @export
# #' @examples


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
  
  
  set_name <- geneset$set_name[keep_term]
  set_size <- c()
  testSet_cor <- backSet_cor <- interCor <- sumTestCor <- rep(NA, length(keep_term))
  
  for ( i in 1:length(keep_term)){
    gset1 <- geneset$gene_set[[keep_term[i]]]
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








