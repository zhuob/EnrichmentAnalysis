#' Compare meaca to other existing methods. 
#' 
#' 
#' @keywords internal
#' @title Compare meaca to existing methods
#' @param dat  the simulated data consisting of three lists:
#'   \item{data}{the expression matrix}
#'   \item{trt}{what treatment each column of \code{data} belongs to}
#'   \item{go_term}{whether each row of \code{data} belongs to this go term}
#'  @param seed simulation seed for reproducibility
#' @return a matrix of p values 
# #' @export
# #' @examples

compare_test <- function(dat, seed){
  
  set.seed(seed)
  ########  a function to incorporate all test procedures
  expression_data <- dat$data
  trt <- dat$trt
  go_term <- dat$go_term
  
  # our test, NO standardization for the simulation
  #  	pval1 <- MOM_test(microarray, trt, go_term, standardize=F)$p             
  
  ## modified on March 16  
  # our test, NO standardization for the simulation
  MEQ <- meaca::meaca_single(expression_data = expression_data, trt = trt, 
                             go_term = go_term, standardize = T)             
  #	print(MEQ)
  p_meaca <- MEQ$p1							# chi-square test
  pval1_2 <- MEQ$p2							# normal reference distribution
  enrich_status <- MEQ$status
  # calculate the three rho's
  cor <- meaca::btw_gene_corr(expression_data = expression_data, trt = trt, 
                              go_term = go_term, standardize = T)
  
  # est_sigma <- estimate_sigma(expression_data, trt)
  # t_val <- est_sigma$t_val
  # pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]            # linear regression
  
  design <- model.matrix(~trt) 
  fit <- limma::lmFit(expression_data, design)
  #stat <- fit$coefficients[, 2]
  # Emperical Bayes t test	
  fit <- limma::eBayes(fit)								                
  #	use the moderated t statistics to do enrichment
  stat <- fit$t[, 2]            							      

  # modified on March 18, 2016
  # This is the default option, see below for explanation.
  #  alter <- "mixed"                                  
  # this is comparable to what we are testing
  alternative <- "either"                            
  # which genes are in GOTERM
  index1 <- which(go_term==1)                       
  
  # up: positive statistics
  # down: negative stats
  # either: the set is either up or down-regulated as a whole
  # mixed: tests whether the genes in the set tend to be DE without regard for direction.
  #        In this case, the test will be significant if the set contains mostly large statistics, negative or positive
  
  # MRGSE
  p_mrgse <- limma::geneSetTest(index = index1, stat, alternative = alternative, ranks.only= T)         
  # sigPathway methods
  p_sigpath <- sig_path(index = index1, stat, nsim = 999)                                                
  # camera proedure
  p_camera <- limma::camera(y = expression_data, index = index1,  
                            design = design, 
                            allow.neg.cor = TRUE, inter.gene.cor = NA)$PValue                                    
  # camera rank 
  p_camera_R <- limma::camera(y = expression_data, index = index1,  
                              design = design, use.ranks= T, 
                              allow.neg.cor = TRUE, inter.gene.cor = NA)$PValue                      
  
  # GSEA
  p_gsea <- GSEA.SingleSet(data = expression_data, trt = trt, go_term = go_term, nperm = 1000)$p.vals		
  
  ## qusage <Yaari, 2013>
  geneSets <- list()
  geneSets[[paste("Set",1)]] <- which(go_term == 1)
  labels <- rep(NA, length(trt))
  labels[trt == 1] <- "B"; labels[trt==0] <- "A"
  # calculate the probability for all genes
  qsarray <- qusage::qusage(expression_data, labels, contrast = "B-A" , geneSets)  			
  # competitive test
  p_qusage <- qusage::pdf.pVal(qsarray, alternative = "two.sided", selfContained = F) 		
  
  
  ## Modified on May 17, 2017
  ## the GSVA method
  # gs <- list(set1 = rownames(expression_data)[geneSets$`Set 1`])
  # es_gsva <- gsva(expression_data, gs, method = "gsva",
  #                     verbose=FALSE, parallel.sz=1)$es.obs
  # p_gsva <- t.test(es_gsva[trt==1], es_gsva[trt==0])$p.value
  
  # over-representation method
  p_ora <- ora(expression_data = dat$data, trt = dat$trt, go_term = dat$go_term, 
               method = "BH", thresh = 0.01) 
  
  p_vec <- c(p_meaca, pval1_2, p_mrgse, p_sigpath, p_camera, p_camera_R,
             p_gsea, p_qusage, p_ora, cor)
  names(p_vec)[1:9] <- c("meaca", "meaca_n", "MRGSE", "SigPathway", 
                          "CAMERA_ModT", "CAMERA_Rank", "GSEA", "QuSAGE", "p_ora")
  p_vec <- data.frame(p_vec)
  p_vec$status <- enrich_status
  return(p_vec)
  
}
