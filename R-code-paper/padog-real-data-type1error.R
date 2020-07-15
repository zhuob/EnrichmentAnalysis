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

## about 2 + 3 minutes to run 100 reps
alpha_meaca <- function(expression_data, trt, go_term, standardize = TRUE, 
                        nrep = 1e3, sim_seed = 123){
  if (standardize == T) {
    expression_data <- meaca::standardize_expression_data(expression_data, trt)
  }
  est_sigma <- meaca:::estimate_sigma(expression_data, trt)
  sigma <- est_sigma$sigma
  t_val0 <- est_sigma$t_val
  
  test_genes <- which(go_term == 1)
  t_test_set <- t_val0[test_genes]
  n_back_set <- sum(go_term == 0)
  
  
  # perform resample of gene level statistics 
  p_alpha <- NULL
  set.seed(sim_seed)
  system.time(
  for(kk in 1:nrep){
    t_val <- t_val0
    t_val[-test_genes] <- base::sample(t_test_set, size = n_back_set, replace = TRUE)
    beta0 <- mean(t_val)
    ones <- rep(1, length(t_val))
    numer <- (crossprod(go_term, (t_val - beta0 * ones)))^2
    go_bar <- rep(mean(go_term), length(go_term))
    denom2 <- t(go_term - go_bar) %*% sigma %*% (go_term - go_bar)
    test_stat <- drop(numer/denom2)
    pval <- 1 - pchisq(test_stat, 1)
    t_temp <- mean(t_val[go_term == 1]) - mean(t_val[go_term == 0])
    if (t_temp > 0) {
      status <- "up"
    }
    else {
      status <- "down"
    }
    test_stat2 <- sqrt(test_stat)
    pval2 <- 2 * min(1 - pnorm(test_stat2), test_stat2)
    p_alpha <- c(p_alpha, pval)
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
  index1 <- which(go_term==1)
  n_back_set <- sum(go_term == 0)

  df_pval <- matrix(nrow = nrep, ncol = 2)
  set.seed(sim_seed)
  for(i in 1:nrep){
    stat_permute <- stat
    stat_permute[-index1] <- base::sample(stat[index1], size = n_back_set, replace = TRUE)
    # MRGSE
    p1 <- limma::geneSetTest(index = index1, stat_permute, 
                             alternative = alternative, ranks.only= T)         
    # sigPathway methods
    p2 <- sig_path(index = index1, stat_permute, nsim = 999)                                                
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
  
  test_genes <- which(go_term == 1)
  n_back_set <- sum(go_term == 0)
  
  p_ora <- NULL
  set.seed(sim_seed)
  for(kk in 1:nrep){
    index0 <- base::sample(test_genes, size = n_back_set, replace = TRUE)
    fit_result <- fit_result0
    fit_result[-test_genes, ] <- slice(fit_result0, index0)
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
                       nrep = 1e3, sim_seed = 123, method = "BH", thresh = 0.1){
  
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



#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                    RUN Simulation
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

df1 <- prep_padog_data("GSE8762")
expression_data <- df1$data; trt <- df1$trt; go_term <- df1$go_term

system.time(result <- alpha_simu(expression_data = expression_data, trt = trt, 
                                 go_term = go_term, standardize = TRUE, 
                                 nrep = 1e3, sim_seed = 123, 
                                 method = "BH", thresh = 0.1))


write_csv(result, "padog-real-data-type1erro-simulation.csv")

