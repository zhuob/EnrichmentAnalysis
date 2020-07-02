#' meaca for single gene set test.
#'
#' @title meaca-single. 
#' @param expression_data  the expressoin matrix.
#' @param trt  treatment indicators, 1 for treatment, 0 for control group.
#' @param go_term  an indicator vector. 1 for genes in the test set, 0 otherwise.
#' @param standardize  whether the data should be standaridzed.
#' @return a list 
#' \item{stat}{the test statistic}
#' \item{p1}{chi-square test p value}
#' \item{status}{"up" or "down", the direction of differential expression}
#' \item{p2}{two-sided test p-value using normal distribution}
#' @export
#' @examples
#' t1 <- simulate_expression_data(size = 50, n_gene = 500, n_test = 100, 
#'                                prop = c(0.1, 0.1), de_mu = 2, de_sd = 1, 
#'                                rho1 = 0.1, rho2 = 0.05, rho3 = -0.05, 
#'                                data_gen_method = "chol", seed = 123)
#' meaca_single(t1$data, trt = t1$trt, go_term = t1$go_term)

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
#' @param min_set_size  the minimum number of genes contained for a gene set to be
#'   considered.
#' @param fdr_method  which method is ued to adjust the p values. see arguments
#'   in function \code{p.adjust}.
#' @return a data frame
#' @export
# #' @examples

meaca_multiple <- function(expression_data, trt, geneset, standardize = T, 
                           min_set_size = 5, fdr_method = "BH"){
  ## conduct enrichment analysis for a battery of gene sets and return the raw p
  ## values as well as adjusted p values
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
  
  keep_term <- which(geneset$size >= min_set_size)
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
    if(!is.finite(t_temp)) {    
      # it's possible that there's no single gene matched between the data and
      # the test set
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
  btw <-btw_gene_corr_multiple(expression_data, trt, geneset, 
                               min_set_size = min_set_size, standardize = T)   
  result_comb <- merge(results, btw, by = "set_name")	
  
  return(result_comb)
}




