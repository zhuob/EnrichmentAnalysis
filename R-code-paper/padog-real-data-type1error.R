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
  resample_test_genes2 <- base::sample(1:length(go_term), size = n_test_set, replace = TRUE)
  resample_back_genes2 <- base::sample(1:length(go_term), size = n_back_set, replace = TRUE)
  resid_mat_hat <- resid_mat
  if(!raw_data){
    resid_mat_hat[test_genes, ] <- resid_mat[resample_test_genes2, ]
    resid_mat_hat[back_genes, ] <- resid_mat[resample_back_genes2, ]
  }
  expression_data_hat <- group_mean_hat + resid_mat_hat
  
  return(expression_data_hat)
  
}

create_bootstrap_data <- function(expression_data, go_term, trt, seed = 123, 
                                  raw_data = FALSE){
  
  test_genes <- which(go_term == 1)
  back_genes <- which(go_term == 0)
  # t_test_set <- t_val0[test_genes]
  n_back_set <- sum(go_term == 0)
  n_test_set <- sum(go_term == 1)
  
  set.seed(seed)
  resample_test_genes <- base::sample(1:length(go_term), size = n_test_set, replace = TRUE)
  resample_back_genes <- base::sample(1:length(go_term), size = n_back_set, replace = TRUE)
  expression_data_hat <- expression_data # initiate the matrix
  if(!raw_data){
    expression_data_hat[test_genes, ] <- expression_data[resample_test_genes, ]
    expression_data_hat[back_genes, ] <- expression_data[resample_back_genes, ]
  }
  
  return(expression_data_hat)
} 


## about 2 + 3 minutes to run 100 reps
alpha_meaca <- function(expression_data, trt, go_term, standardize = TRUE, 
                        nrep = 1e3, sim_seed = 123){
  if (standardize == T) {
    expression_data <- meaca::standardize_expression_data(expression_data, trt)
  }
  
  # perform resample of gene level statistics 
  p_alpha <- NULL
  set.seed(sim_seed)
  system.time(
  for(kk in 1:nrep){
    cat("\r", kk)
    expression_data_hat <- create_bootstrap_data(expression_data = expression_data, 
                                                 go_term = go_term, trt = trt, 
                                                 seed = kk, raw_data = FALSE)
    temp <- meaca::meaca_single(expression_data = expression_data_hat, 
                                trt = trt, go_term = go_term, 
                                standardize = FALSE)
    p_alpha <- c(p_alpha, temp$p1)
  }
  )
  return(p_alpha)
}

alpha_mrgse_sigpath <- function(expression_data, trt, go_term,  
                        nrep = 1e3, sim_seed = 123){
  
  design <- model.matrix(~trt) 
  fit <- limma::lmFit(expression_data, design)
  #stat <- fit$coefficients[, 2]
  # Emperical Bayes t test	
  fit <- limma::eBayes(fit)								                
  #	use the moderated t statistics to do enrichment
  stat <- fit$t[, 2]            							      
  alternative <- "either"                            
  # which genes are in GOTERM
  test_genes <- which(go_term == 1)
  back_genes <- which(go_term == 0)
  n_back_set <- sum(go_term == 0)
  n_test_set <- sum(go_term == 1)
  
  df_pval <- matrix(nrow = nrep, ncol = 2)
  set.seed(sim_seed)
  for(i in 1:nrep){
    stat_permute <- stat
    resample_test_genes <- base::sample(1:length(go_term), size = n_test_set, replace = TRUE)
    resample_back_genes <- base::sample(1:length(go_term), size = n_back_set, replace = TRUE)
    stat_permute[test_genes] <- stat[resample_test_genes]
    stat_permute[back_genes] <- stat[resample_back_genes]
    # MRGSE
    p1 <- limma::geneSetTest(index = test_genes, stat_permute, 
                             alternative = alternative, ranks.only= T)         
    # sigPathway methods
    p2 <- sig_path(index = test_genes, stat_permute, nsim = 999)                                                
    df_pval[i, ] <- c(p1, p2)
  }
  names(df_pval) <- c("p_mrgse", "p_sigpath")
  df_pval <- data.frame(df_pval)
  return(df_pval)
}

# 4 minutes to run 100 simulations
# system.time(r2 <- alpha_mrgse_sigpath(expression_data = expression_data, trt = trt, 
#                                       go_term = go_term, nrep = 100, sim_seed = 123))

## ora 

alpha_ora <- function(expression_data, trt, go_term, method = "BH", 
                      sim_seed = 123, thresh = 0.01, nrep = 1e3){
  
  design <- model.matrix(~trt) 
  fit <- limma::lmFit(expression_data, design)
  #stat <- fit$coefficients[, 2]
  # Emperical Bayes t test	
  fit <- limma::eBayes(fit)								                
  #	use the moderated t statistics to do enrichment
  stat <- fit$t[, 2]            							      
  
  fit_result0 <- limma::topTable(fit, number = length(fit$p.value), 
                                sort.by = "none", adjust.method = method)
  
  # identify which genes are from which group
  test_genes <- which(go_term == 1)
  back_genes <- which(go_term == 0)
  n_back_set <- sum(go_term == 0)
  n_test_set <- sum(go_term == 1)
  
  p_ora <- NULL
  set.seed(sim_seed)
  for(kk in 1:nrep){
    fit_result <- fit_result0
    
    resample_test_genes <- base::sample(1:length(go_term), size = n_test_set, replace = TRUE)
    resample_back_genes <- base::sample(1:length(go_term), size = n_back_set, replace = TRUE)
    fit_result[test_genes, ] <- slice(fit_result0, resample_test_genes)
    fit_result[back_genes, ] <- slice(fit_result0, resample_back_genes)
    
    fit_result <- fit_result %>% mutate(adj.P.Val = p.adjust(adj.P.Val, method = method))
    
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
    
    p_ora <- c(p_ora, p_enrich)
  }
  
  return(p_ora)

}

###  
# system.time(r3 <- alpha_ora(expression_data = expression_data, trt = trt, 
#                             go_term = go_term, nrep = 100, sim_seed = 123, 
#                             method = "BH", thresh = 0.01))

# The GSEA is too complicated, let's wait first
# p_gsea <- GSEA.SingleSet(data = expression_data, trt = trt, go_term = go_term, nperm = 1000)$p.vals		

alpha_simu <- function(expression_data, trt, go_term, standardize = TRUE,
                       nrep = 1e3, sim_seed = 123, method = "BH", thresh = 0.01){
  
  m1 <- alpha_meaca(expression_data, trt, go_term, standardize = standardize, 
                    nrep = nrep, sim_seed = sim_seed)
  m2 <- alpha_mrgse_sigpath(expression_data, trt, go_term, 
                            nrep = nrep, sim_seed = sim_seed)
  m3 <- alpha_ora(expression_data = expression_data, trt = trt, 
                  go_term = go_term, nrep = nrep, sim_seed = sim_seed,
                  method = method, thresh = thresh)

  result <- bind_cols(m1, m2, m3)
  names(result) <- c("p_meaca", "p_mrgse", "p_sigpath", "p_ora")
  return(result)
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
  pvals <- foreach(kk = 1:nsim, .combine = bind_rows, 
                   .packages = package_used, 
                   .verbose = verbose_show) %dopar% {
    
    expression_data_hat <- create_bootstrap_data(expression_data = expression_data, 
                                                 go_term = go_term, trt = trt, 
                                                 seed = kk, raw_data = FALSE)
    
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

result <- compare_test_new(dat = df1, seed = 2, nsim = 1000, ncore = 23)

system.time(result <- alpha_simu(expression_data = expression_data, trt = trt, 
                                 go_term = go_term, standardize = TRUE, 
                                 nrep = 1e3, sim_seed = 123, 
                                 method = "BH", thresh = 0.01))


write_csv(result, "padog-real-data-type1error-simulation-all-genes-v2.csv")

