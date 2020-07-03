#' implement the Sigpathay. 
#' 
#' @keywords internal
#' 
#' @title Sigpathway
#' @param num_gene  number of genes to be simulated.
#' @return p   the p-values.
#' @export
#' @examples

sig_path <- function(index, statistics, nsim = 9999){
  
  ssel <- statistics[index]
  ssel <- ssel[!is.na(ssel)]
  nsel <- length(ssel)
  if(nsel==0) return(1)
  stat <- statistics[!is.na(statistics)]
  n_gene <- length(statistics)
  msel1 <- mean(ssel)
  msel2 <- mean(statistics[-index])
  msel <- msel1  - msel2              ## use the difference as the test statistic
  
  ntail <- 1                          ## permuation test
  permu_sel <- c()
  for (i in 1:nsim) {
    samp_id <- sample(1: n_gene, nsel)
    permu_stat1 <- mean(stat[samp_id])
    permu_stat2 <- mean(stat[-samp_id])
    permu_sel[i] <- permu_stat1 - permu_stat2  # calculate the difference
    
    if (abs(permu_sel[i]) >= abs(msel)){
      ntail <- ntail + 1
    }
  }
  p_value <- ntail/(nsim+1)
  as.vector(p_value)
}



#' multiple gene set enrichemnt analysis
#' 
#' 
#' @keywords internal
#' @title Camera multiple
#' @param expression_data the expression data matrix
#' @param trt the treatment lables for each column of \code{expression_data}
#' @param  geneset  a list containing all the gene sets to be tested for
#'   enrichment status
#' @param use.rank  If \code{TRUE} then it corresponds to "Camara-Rank",
#'   otherwise "Camara"
#' @return a list that has exactly the same elements as that in \code{camara} of
#'   "limma" package (see user manual for more details)
# #' @export
Camera_multiple <- function(expression_data, trt, geneset, use.rank = F){
  ##  Use Camera procedure to do a battery of gene set test. 
  
  # for GSE64810 data
#  microarray <- expression[, -(1:2)]
#  all_genes <-  expression[, 2]
  
  # for other typical data, where row names are genes
          
  all_genes <- rownames(expression_data)
  
  gset1 <- list()
  set.name <- c()
  for ( i in 1:length(geneset$size)){
    gset1[[i]] <- geneset$gene_set[[i]][-(1:2)]
    set.name[i] <- geneset$gene_set[[i]][1]
  }
  names(gset1) <- set.name
  # it contains multiple lists
  c2.indices <- limma::ids2indices(gset1, all_genes)      
  
  design <- model.matrix(~trt)
  Results <- limma::camera(expression_data, c2.indices,  design, 
                           use.ranks = use.rank, sort = F)
  return(Results)
  
}




#' multiple gene set enrichemnt analysis 
#' 
#' 
#' @keywords internal
#' @title  MRGSE multiple
#' @param expression the expression data matrix
#' @param trt the treatment lables for each column of \code{expression_data}
#' @param  geneset  a list containing all the gene sets to be tested for enrichment status
#' @param use.rank  If \code{TRUE} then it corresponds to "MRGSE", otherwise "Camara"
#' @return a matrix
# #' @export

MRGSE_multiple <- function(expression, trt, geneset, use.rank = TRUE){
  
  
  # for GSE64810 data
  #        microarray <- expression[, -(1:2)]
  #        all_genes <-  expression[, 2]
  
  # for other typical data, where row names are genes
  microarray <- expression
  all_genes <- rownames(microarray)
  
  gset1 <- list()
  set.name <- c()
  for ( i in 1:length(geneset$size)){
    gset1[[i]] <- geneset$gene_set[[i]][-(1:2)]
    set.name[i] <- geneset$gene_set[[i]][1]
  }
  names(gset1) <- set.name
  # it contains multiple lists
  c2.indices <- limma::ids2indices(gset1, all_genes)	 
  
  design <- model.matrix(~trt)
  fit <- limma::lmFit(microarray, design)
  # Emperical Bayes t test    
  fit <- eBayes(fit)                                                                              
  # use the moderated t statistics to do
  stat <- fit$t[, 2]                                                                    
  alter <- "either"
  
  n_genes <- nrow(microarray)
  set.name <- names(c2.indices)
  set.size <- p.MRGSE <- rep(0, length(set.name))
  for ( i in 1:length(c2.indices)) {
    index1 <- c2.indices[[i]]
    # print(i)
    set.size[i] <- length(index1)
    p.MRGSE[i] <- limma::geneSetTest(index = index1, stat, alternative = alter, ranks.only= use.rank)
  }
  
  result <- data.frame(set.name, size = set.size, p.MRGSE)
  return(result)
}




#' The PLAGE method 
#' 
#' @param expression_data the expression matrix
#' @param trt treatment label, 1 for treatment and 0 for control
#' @param go_term gene go term indicator, 1 for genes in the test set, 0 otherwise
#' @param nperm   Number of permutation to run (to get the p value by permuted results t statistics [see equation 4 of the paper]  )
#' @param seed for reproducibility purpose
#' @return pval 
#' @export
#'   

plage <- function(expression_data, trt, go_term, nperm = 999, seed = 125){
  
  rownames(expression_data) <- paste("Gene", 1:nrow(expression_data), sep = "")
  gs <- list( set1 =  paste("Gene", which(go_term==1), sep = ""))
  
  ## calculate the activity level for each sample
  es_plage <- plage_score(expression_data, gs)
  # get the t-statistic for each treatment
  es_a <- es_plage[trt== 1]; 
  es_b <- es_plage[trt==0] 
  na <- length(es_a); nb <- length(es_b); n_total <- na + nb
  t_obs <- abs(mean(es_a)- mean(es_b))/sqrt(var(es_a)/na + var(es_b)/nb)
  
  t_perm <- rep(NA, nperm)
  # do the permutations
  for (k in 1:nperm){
    na_perm <- sample(1:n_total, na)
    nb_perm <- setdiff(1:n_total, na_perm)
    expression_data_perm <- expression_data[,c(na_perm, nb_perm)]
    es_plage_perm <- plage_score(expression_data_perm, gs)
    es_a_perm <- es_plage_perm[trt== 1]; 
    es_b_perm <- es_plage_perm[trt==0] 
    t_perm[k] <- abs(mean(es_a_perm)- mean(es_b_perm))/sqrt(var(es_a_perm)/na + var(es_b_perm)/nb)
  }
  
  # the pvalue is (b + 1)/(nperm + 1)
  b <- sum(t_perm >= t_obs)
  pval <- (b + 1)/(nperm + 1)
  return(pval)
}


###    this code is extracted from GSVA package 
plage_score <- function(X, gs){
  geneSets <- gs$set1
  p <- nrow(X)
  n <- ncol(X)
  
  Z <- t(apply(X, 1, function(x) (x-mean(x))/sd(x)))
  
  rightsingularsvdvectorgset <- function(gSetIdx, Z) {
    s <- svd(Z[gSetIdx, ])
    s$v[, 1]
  }
  
  es<-  rightsingularsvdvectorgset(geneSets, Z)
  names(es) <- colnames(X)
  
  return(es)
}


#' the over-representation method (see http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079217 ). 
#' 
#' @param expression_data  the expressoin matrix.
#' @param trt  treatment indicators, 1 for treatment, 0 for control group.
#' @param go_term  an indicator vector. 1 for genes in the test set, 0 otherwise.
#' @param method The method used to adjust for nominal p values if necessary.
#'   \code{BH} by default
#' @param thresh the threshold described in the paper, 200 by default.
#' @return pval  the enrichment p value.
#' 
ora <- function(expression_data, trt, go_term, method = "BH", thresh = 0.01){
  # the implementation of ORA method is described in this paper
  # http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079217
  
  design <- model.matrix(~trt) 
  fit <- limma::lmFit(expression_data, design)
  #stat <- fit$coefficients[, 2]
  # Emperical Bayes t test	
  fit <- limma::eBayes(fit)								                
  #	use the moderated t statistics to do enrichment
  stat <- fit$t[, 2]            							      
  
  fit_result <- limma::topTable(fit, number = length(fit$p.value), 
                         sort.by = "none", adjust.method = method)
  fold_change <- abs(fit_result$logFC)  ## the absolute fold change
  pvals <- fit_result$P.Value
  adjust_p <- fit_result$adj.P.Val
  
  thresh1 <- thresh
  if (thresh < 1 ) {
    thresh1 <- ceiling(thresh*length(adjust_p))
  }
  
  # procedure 1
  if (sum(adjust_p < 0.1)> thresh1) {
    DE <- adjust_p < 0.1
    cont <- ftable(data.frame(go_term, DE))
  }
  else {
    if (sum(pvals < 0.05 & fold_change > 1.5) > thresh1){ ## procedure 2
      DE <- pvals < 0.05 & fold_change > 1.5
      cont <- ftable(data.frame(go_term, DE))
    } else {    # procedure 3
      cutoff <- quantile(pvals, 0.01)
      DE <- pvals <= cutoff
      cont <- ftable(data.frame(go_term, DE))
    }
  }
  
  k <- colSums(cont)[2]  # the number of DE genes
  ros <- rowSums(cont)
  m <- ros[2]  # the number of GO term genes
  n <- ros[1]  # the number of background genes
  q <- cont[2, 2] # the number of DE genes from the GO term 
  
  # the probability of getting more than q DE genes
  # p_enrich <- phyper(q = q-1, m = m, n = n, k = k, lower.tail =  F) 
  ## BZ made change on June 11, 2017
  p_enrich <- phyper(q = q, m = m, n = n, k = k, lower.tail =  F)  # P(X > q)
  
  return(p_enrich)
}






