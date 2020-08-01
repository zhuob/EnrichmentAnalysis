#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## this file is to run for real-data based type 1 error simulation
##   Method:  use real data, extract test statistics used in different methods
##            use real data to derive between-gene correlation
##            but the gene level statistics for the background gene set is sampled
##            from the test set statitics (of course, with replacement)
##            Therefore, we create a hypothetical null, based on which we 
##            evaluate the model performance on type 1 error
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 

create_bootstrap_data <- function(expression_data, go_term, trt, seed = 123, 
                                  raw_data = FALSE){
  
  expression_data <- as.matrix(expression_data)
  group_mean <- as.matrix(meaca:::trt_mean(expression_data, trt))
  resid_mat <- expression_data - group_mean
  
  test_genes <- which(go_term == 1)
  back_genes <- which(go_term == 0)
  # t_test_set <- t_val0[test_genes]
  n_back_set <- sum(go_term == 0)
  n_test_set <- sum(go_term == 1)
  
  ## bootstrap beta0, beta1
  
  set.seed(seed)
  resample_test_genes <- base::sample(1:length(go_term), size = n_test_set, replace = TRUE)
  resample_back_genes <- base::sample(1:length(go_term), size = n_back_set, replace = TRUE)
  group_mean_hat <- group_mean # initiate the matrix
  if(!raw_data){
    group_mean_hat[test_genes, ] <- group_mean[resample_test_genes, ]
    group_mean_hat[back_genes, ] <- group_mean[resample_back_genes, ]
  }

  # bootstrap residual independently from beta0 and beta1
  set.seed(seed + 1e5)
  resample_test_genes2 <- base::sample(test_genes, size = n_test_set, replace = TRUE)
  resample_back_genes2 <- base::sample(back_genes, size = n_back_set, replace = TRUE)
  resid_mat_hat <- resid_mat
  if(!raw_data){
    resid_mat_hat[test_genes, ] <- resid_mat[resample_test_genes2, ]
    resid_mat_hat[back_genes, ] <- resid_mat[resample_back_genes2, ]
  }
  expression_data_hat <- group_mean_hat + resid_mat_hat
  
  return(expression_data_hat)
  
}



compare_test_new <- function(dat, seed, nsim = 1e3, ncore = 4, 
                             package_used = c("tidyverse"), verbose_show = FALSE){
  
  ########  a function to incorporate all test procedures
  expression_data <- dat$data
  trt <- dat$trt
  go_term <- dat$go_term

  library(foreach)
  cl <- parallel::makeCluster(ncore, type = "FORK")
  doParallel::registerDoParallel(cl)
  set.seed(seed)
  subseed <- sample(1:1e6, nsim)
  ## add a standardization step 
  expression_data <- meaca::standardize_expression_data(expression_data, trt)
  
  pvals <- foreach(kk = 1:nsim, .combine = bind_rows, 
                   .packages = package_used, 
                   .verbose = verbose_show) %dopar% {
    
    expression_data_hat <- create_bootstrap_data(expression_data = expression_data, 
                                                 go_term = go_term, trt = trt, 
                                                 seed = subseed[kk], raw_data = FALSE)
    
    expression_data_hat <- meaca::standardize_expression_data(expression_data_hat, trt)
    
    temp <- meaca::meaca_single(expression_data = expression_data_hat, 
                                trt = trt, go_term = go_term, 
                                standardize = FALSE)
    p_meaca <- temp$p1
    
    design <- model.matrix(~trt) 
    fit <- limma::lmFit(expression_data_hat, design)
    # Emperical Bayes t test	
    fit <- limma::eBayes(fit)								                
    #	use the moderated t statistics to do enrichment
    stat <- fit$t[, 2]            							      

    alternative <- "either"                            
    # which genes are in GOTERM
    index1 <- which(go_term==1)                       
    
    # MRGSE
    p_mrgse <- limma::geneSetTest(index = index1, stat, alternative = alternative, ranks.only= T)         
    # sigPathway methods
    p_sigpath <- sig_path(index = index1, stat, nsim = 999)                                                
    # camera proedure
    p_camera <- limma::camera(y = expression_data_hat, index = index1,  
                              design = design, 
                              allow.neg.cor = TRUE, inter.gene.cor = NA)$PValue                                    
    # camera rank 
    p_camera_R <- limma::camera(y = expression_data_hat, index = index1,  
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
    qsarray <- qusage::qusage(expression_data_hat, labels, contrast = "B-A" , geneSets)  			
    # competitive test
    p_qusage <- qusage::pdf.pVal(qsarray, alternative = "two.sided", selfContained = F) 		
    
    p_plage <- plage(expression_data = expression_data_hat, trt = trt, 
                     go_term = go_term, seed = kk, nperm = 999)
    # modified on May 14, 2017
    # the over-representation method
    p_ora <- ora(expression_data = expression_data_hat, trt = trt, 
                 go_term = go_term, method = "BH", thresh = 0.01) 
    
    temp1 <- tibble::tibble(p_meaca = p_meaca, 
                            p_mrgse = p_mrgse, 
                            p_sigpath = p_sigpath,
                            p_camera = p_camera, 
                            p_camera_R = p_camera_R,
                            p_gsea  = p_gsea,
                            p_qusage = p_qusage, 
                            p_plage = p_plage, 
                            p_ora = p_ora)
    
      return(temp1)
    }
  parallel::stopCluster(cl)
  
  return(pvals)
}

#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                    RUN Simulation
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

df1 <- prep_padog_data("GSE8762")
expression_data <- df1$data; trt <- df1$trt; go_term <- df1$go_term

result <- compare_test_new(dat = df1, seed = 1234, nsim = 1000, ncore = 4)

write_csv(result, "real-data/padog-package/padog-real-data-type1error-simulation-all-genes-data-corr-test-bootstrap-separately-v3.csv")

