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





#' Compare meaca to other existing methods. 
#' 
#' 
#' @keywords internal
#' @title Compare meaca to existing methods
#' @param dat  the simulated data consisting of three lists:
#'   \item{data}{the expression matrix}
#'   \item{trt}{what treatment each column of \code{data} belongs to}
#'   \item{go_term}{whether each row of \code{data} belongs to this go term}
#' @return a matrix of p values 
# #' @export
# #' @examples

compare_test <- function(dat){
  ########  a function to incorporate all test procedures
  
  expression_data <- dat$data
  trt <- dat$trt
  go_term <- dat$go_term
  
  #  	pval1 <- MOM_test(microarray, trt, go_term, standardize=F)$p             
  # our test, NO standardization for the simulation
  
  ## modified on March 16  
  MEQ <- meaca_single(expression_data, trt, go_term, standardize = F)             
  # our test, NO standardization for the simulation
  #	print(MEQ)
  pval1 <- MEQ$p1							# chi-square test
  pval1_2 <- MEQ$p2							# normal reference distribution
  
  est_sigma <- estimate_sigma(expression_data, trt)
  t_val <- est_sigma$t_val
  pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]            # linear regression
  
  
  design <- model.matrix(~trt) 
  fit <- limma::lmFit(expression_data, design)
  #stat <- fit$coefficients[, 2]
  
  # Emperical Bayes t test	
  fit <- limma::eBayes(fit)								                
  #	use the moderated t statistics to do enrichment
  stat <- fit$t[, 2]            							      
  
# modified on March 18, 2016
  #  alter <- "mixed"                                  
  # This is the default option, see below for explanation.
  # this is comparable to what we are testing
  alternative <- "either"
  # which genes are in GOTERM
  index1 <- which(go_term==1)                       
  
  # up: positive statistics 
  # down: negative stats either: the set is either up or
  # down-regulated as a whole mixed: tests whether the genes in the set tend to
  # be DE without regard for direction. In this case, the test will be
  # significant if the set contains mostly large statistics, negative or
  # positive
  
   # MRGSE
  tes1 <- limma::geneSetTest(index = index1, stat, alternative = alternative, ranks.only= T)    
  # sigPathway methods
  tes2 <- sig_path(index = index1, stat)
  # camera proedure
  tes3 <- limma::camera(expression_data, index1,  design)$PValue 
  # camera rank 
  tes4 <- limma::camera(expression_data, index1,  design, use.ranks= T)$PValue                      
  
  tes5 <- GSEA.SingleSet(dat$data, dat$trt, dat$go_term, nperm=1000)$p.vals		# GSEA
  
  
  ## qusage <Yaari, 2013>
  geneSets <- list()
  geneSets[[paste("Set",1)]] <- which(go_term == 1)
  labels <- rep(NA, length(trt))
  labels[trt == 1] <- "B"; labels[trt==0] <- "A"
  # calculate the probability for all genes
  qsarray <- qusage::qusage(expression_data, labels, contrast = "B-A" , geneSets)  			
  # competitive test
  tes6 <- qusage::pdf.pVal(qsarray, alternative = "two.sided", selfContained = F) 		
  
  p_vec <- c(pval1, pval1_2, pval2, tes1, tes2, tes3, tes4, tes5, tes6)
  names(p_vec) <- c("meaca", "meaca_n", "LM", "MRGSE", "SigPathway", 
                    "CAMERA_ModT", "CAMERA_Rank", "GSEA", "QuSAGE")
  return(p_vec)
  
  
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

